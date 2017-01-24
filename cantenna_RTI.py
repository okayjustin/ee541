"""
Produces an RTI (range x time intensity) image of the cantenna recording.
Also applies a simple two-pulse canceller to filter out clutter and CFAR
normalization to improve visual detection. See inside the function for
for additional parameters.

Pythonized code for Raytheon Low Cost Radar Initiative
Adapted from:
MIT Short Program 2013 Laptop Radar Course
(c) 2014 Massachusetts Institute of Technology
Version 1 12/6/2013
Version 2 8/11/2015
    Updated print statement to function for Python 3.4 compatibility.
    Corrected inverted y-axis in plot.
Version 3 6/29/2016
    Updated style to PEP8 standards.
    Added channel check.
"""

import common as com
import constants as C
import scipy.signal as sig
import os.path
from matplotlib import pyplot as plt
from numpy import (append, arange, argwhere, array, delete, diff, exp, floor,
 log10, median, nonzero, ones, sign, sin, tile, zeros)
from scipy import fftpack
from time import time


# Inputs: 
#    wavFile  The .wav filename as string.
def cantenna_rti_v6(wavFile):
    # Use default wave file if wavFile does not exist.
    if not os.path.isfile(wavFile):
        print("File %s not found." % (wavFile))
        wavFile = 'data/running_outside_20ms.wav'
        print("Using %s" % wavFile)
        
    [fs, y] = com.read_wavefile(wavFile)
    [data, beat] = com.chan_check(y, 'RTI')
    # Samples per pulse.
    #nP = round(C.TP * fs)   
    # Transmit bandwdith.
    bw      = C.FSTOP - C.FSTART
    # Range resolution.cantenna_dop_v4(wavFile)
    delta_r = C.C / (2 * bw)

    # The trigger signal is in the first channel.    
    tr_sig = -1*(beat)
    # The mixer output is in the second channel.
    s      = -1*(data)
    
    print("Preprocess Sync...")
    trig = com.schmitt_trigger(tr_sig, C.SCHMITT_THRESH)

    # Parse the trigger signal (look for threshold crossings).
    print("Parsing the recording...")
    riseFallGroup = com.process_sync(trig)
    pulseStarts   = riseFallGroup[0]

    # Refine using measured parameters: the pulse width.
    nP = round(min(diff(pulseStarts)) / 2)
    tp = float(nP) /float(fs)
    print("Measured pulsewidth of %f ms." % (tp * 1000.0))

    # Pre-compute some windows and other vectors.
    # Number of output range samples.
    nRange     = int(floor(C.OSAMP_RNG * nP / 2))
    # Range of each bin.
    dataRange  = arange(0.0, nRange, 1.0) * (delta_r / C.OSAMP_RNG)
    # Keep this for plotting purposes.
    rangeMax   = max(dataRange)
    dataRgIdx  = dataRange > C.MAX_RANGE
    dataRange  = delete(dataRange, argwhere(dataRgIdx))
    # Number of range bins to keep.
    nRangeKeep = len(dataRange)           
    # The window applied to reduce range sidelobes.
    rngWin     = com.hann_window(1, nP + 1, nP + 1)  
    # The window applied to the padded data.
    padWin     = (sin(array(range(1, C.NUM_PAD + 1), 'float')
                     / (C.NUM_PAD + 1) * C.PI/2)**2)
    # The window applied to the trigger data.
    trgWin     = com.hann_window(1, C.NUM_PAD*2 + 2, C.NUM_PAD*2 + 2)

    # Get data parameters, e.g. number of pulses.
    nSamples   = len(s)

    pulseStarts = pulseStarts[pulseStarts + nP + C.NUM_PAD <= nSamples]
    numPulses   = pulseStarts.size
    print("Found %d pulses." % (numPulses))

    # Process pulses into a data matrix.
    sif         = zeros((nRangeKeep, numPulses), dtype=complex)
    print("Processing pulse data...")

    for pIdx in range(numPulses):
        # Bandlimited interpolate the trigger signal.
        tmp        = (trig[pulseStarts[pIdx] + 
                     arange(-C.NUM_PAD,C.NUM_PAD+1)] * trgWin)
        interpTmp  = com.fft_interp(tmp, C.OSAMP_TRG)
        interpTmp  = (interpTmp[(C.NUM_PAD*C.OSAMP_TRG )
                      + arange(-C.OSAMP_TRG, C.OSAMP_TRG+1)])
        interpOffs = (arange(-C.OSAMP_TRG, C.OSAMP_TRG + 1, dtype=float)
                      / C.OSAMP_TRG)
        myIdx      = nonzero(diff(sign(interpTmp)) == 2)
        myIdx      = myIdx[0] + 1
        myIdx      = append(myIdx -1, myIdx)
        if 0 == len(myIdx):
            print("Trigger edge not found... skipping pulse %d", pIdx)
            continue

        tmp2 = interpTmp[myIdx]
        # Linear interpolate to find the zero crossing.
        fracOffset = (-(interpOffs[myIdx[1]] - tmp2[1] / (tmp2[1]-tmp2[0])
                    /  C.OSAMP_TRG))
        # Time-align the data to the trigger event (the zero crossing).
        cInds            = pulseStarts[pIdx] + arange(-C.NUM_PAD,(nP +C.NUM_PAD))
        tmp              = s[cInds.astype(int)]
        tmp              = tmp.astype(float)
        tmp[0:C.NUM_PAD] = tmp[0:C.NUM_PAD] * padWin
        tmp[arange(len(tmp) - 1, len(tmp) - C.NUM_PAD-1, -1)] =  \
            tmp[arange(len(tmp) - 1, len(tmp) - C.NUM_PAD - 1, -1)] * padWin
        # Time delay applied in the frequency domain below.
        tmp          = fftpack.fft(tmp)
        tmp          = (tmp * exp( -1j*(arange(0, nP + 2 * C.NUM_PAD))
            / (nP + 2 * C.NUM_PAD) * 2 * C.PI * fracOffset))

        tmp          = fftpack.ifft(tmp)
        # Compute & scale range data from the time-aligned mixer output.
        tmp          = (fftpack.ifft(tmp[C.NUM_PAD + arange(0, nP, dtype=int)]
                       * rngWin, 2*nRange))
        sif[0:,pIdx] = tmp[0:nRangeKeep]

    sif = sif.T

    # Display the RTI.
    plt.figure()
    plt.imshow(20*log10(abs(sif)), cmap='jet', \
        extent=[0, rangeMax, 2*numPulses*tp, 0], aspect='auto')
    plt.title('RTI without clutter rejection')
    plt.xlabel('Range (m)')
    plt.ylabel('Time (s)')
    plt.colorbar()


    # Apply the N-pulse canceller.
    mti_filter         = -ones((C.NPULSE_CANCEL, 1)) / C.NPULSE_CANCEL
    midIdx             = int(floor((C.NPULSE_CANCEL + 1) / 2))
    mti_filter[midIdx] = mti_filter[midIdx] + 1

    # Temp used since 'same' option behavior differs between Python and Matlab.
    # Code below is equivalent to conv(sif, mit_filter, 'same') in Matlab.
    
    temp = sig.convolve(sif, mti_filter)
    sif  = temp[range(1, temp.shape[0]), :]
    # Display the MTI results.
    plt.figure()
    plt.imshow(20*log10(abs(sif)), cmap='jet', \
        extent=[0, rangeMax, 2*numPulses*tp, 0], aspect='auto')
    plt.title('RTI with MTI clutter rejection')    
    plt.xlabel('Range (m)')    
    plt.ylabel('Time (s)')
    plt.colorbar()

    # Apply the median CFAR normalization. 
    sif_dB = 20*log10(abs(sif))
    med    = tile(median(sif_dB,0), numPulses).reshape(numPulses, len(dataRange))
    sif_dB = sif_dB - med

    med    = None
    med    = tile(median(sif_dB,1), sif.shape[1]).reshape(len(dataRange), numPulses)
    med    = med.T
    sif_dB = sif_dB - med
    
    # Plot the CFAR normalized results.
    plt.figure()
    plt.imshow(sif_dB, cmap='jet', extent=[dataRange[0], 
        rangeMax, 2*numPulses*tp, 0], aspect='auto', vmin=-3, vmax=37)
    plt.title('RTI with MTI+CFAR')
    plt.xlabel('Range (m)')
    plt.ylabel('Time (s)')
    plt.colorbar()

# If specifying path to a data file, the file separator is a forward slash, 
# even on Windows. For example, see commented line below.
# wavFile = 'Z:/rasPI/data/running_outside_20ms.wav'
wavFile = 'running_outside_20ms.wav'
# Start timer.
startTime = time()
cantenna_rti_v6(wavFile)
stopTime  = time()
print("Total time was", stopTime - startTime , "seconds")



