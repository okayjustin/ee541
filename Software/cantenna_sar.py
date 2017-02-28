"""
Produces a SAR image of the cantenna recording. This processing assumes the
data was collected by switching the radar on/off at several positions.
See inside the function for additional parameters.

Pythonized code for Raytheon Low Cost Radar Initiative
Adapted from:
MIT Short Program 2013 Laptop Radar Course
(c) 2014 Massachusetts Institute of Technology
Version 1 1/3/2014
Version 2 8/11/2015
    Updated print statement to function for Python 3.4 compatibility
Version 3 7/15/2016 
    Updated format for PEP8 standard. Add chan_check call.
    
"""


import common as com
import constants as C
import gc
import matplotlib.pyplot as plt
import os.path
from numpy import (append, arange, array, diff, digitize, exp, fft, flipud, 
floor, isnan, linspace, log10, negative, ones, sign, sin, sqrt, tan, where, 
zeros)
 
def cantenna_sar(wavFile):
    if not os.path.isfile(wavFile):
       wavFile = 'data/towardswarehouse.wav'
    [fs, y] = com.read_wavefile(wavFile)
    # Confirm which channel has the data.
    [data, beat] = com.chan_check(y, 'SAR')

    # Transmit bandwidth in Hertz.
    bw      = C.FSTOP - C.FSTART 
    # Range resolution in meters.
    delta_r = C.C / (2.0 * bw)
    # Center frequency in Hertz.
    fc      = (C.FSTART + C.FSTOP) / 2.0 
    # Image pixels in x (azimuth) dimension.
    imgAz   = linspace(-C.SCENE_SZ_AZ / 2, C.SCENE_SZ_AZ / 2, 641)
    # Image pixels in y-dim
    imgRg   = linspace(0, C.SCENE_SZ_RG, 481)
  
    # The trigger signal is in the first channel.  
    tr_sig  = negative(beat)
    # The mixer output is in the second channel.
    s       = negative(data)

    print("Preprocess Sync...")
    trig    = com.schmitt_trigger(tr_sig, C.SCHMITT_THRESH)

    gc.collect()
    # Parse the trigger signal (look for threshold crossings).
    print("Parsing the recording...")
    #com.process_sync(trig) inputs trig variable which has data captured
    #captured from the schmitt trigger
    #riseFallGroup is a 2-D array with [rise, fall, groupNum] all representing
    #arrays
    riseFallGroup = com.process_sync(trig)
    pulseStarts   = riseFallGroup[0]
    breakNumber   = riseFallGroup[2]
     
    # Rework of Matlab's histc function.    
    currentBreak = breakNumber[len(breakNumber) - 1]
    bins         = digitize(breakNumber, range(currentBreak + 1), right=False)
    plsPerBreak  = zeros(bins[len(bins) - 1], dtype=int)
    for i in range(1, bins[len(bins) - 1] + 1):
        plsPerBreak[i-1] = sum(i == bins)
    
    numBreaks = currentBreak
    numPulses = pulseStarts.size
    print("Found %d breaks and %d pulses." % (numBreaks, numPulses))
    
    # Refine using measured parameters: the pulsewidth.
    np = round(max(diff(pulseStarts)) / 2)
    for bIdx in range(2, numBreaks + 1):
        myPulses = where(breakNumber == bIdx)[0]
        # function diff() produces out[n] = a[n+1] - a[n]
        # tempMin is storing the average min time
        # error: TypeError: only length-1 arrays can be converted to Python scalars
        tempMin  = round(diff(pulseStarts[myPulses[1:len(myPulses) - 1]]) / 2).min()
        np       = array([tempMin, np]).min()

    tp = np / fs
    print("Measured pulsewidth of %2.3f ms." % (tp*1e3))
    
    minpulsesPerBreak = plsPerBreak[1:len(plsPerBreak) - 1].min()
    print("Minimum of %d pulses between breaks" % minpulsesPerBreak)

    # Compute aperture, range, and cross-range dimensions.
    # Cross range position of radar on aperture L in meters.
    xa = (range(1, numBreaks) - (numBreaks - 2.0) / 2.0) * C.DELTA_AZ
    # Range position of radar, assumes aperture is along a straight line.
    ya = 0 * xa
    # Radar altitude, assumes aperture is at 0 m elevation relative to scene.
    za = 0 * xa

    # Make figure for output display.
    myImg = zeros([imgRg.size, imgAz.size], dtype='complex')

    # Pre-compute some windows and other vectors.
    # Ranges of each bin.    
    dataRange  = arange(C.OSAMP_FACT * np) * (delta_r / C.OSAMP_FACT) 
    # Scale factor applied to the data as a fn of range.
    rangeScale = dataRange**(3.0 / 2.0)
    # Apply window to reduce range sidelobes
    rngWin     = com.hann_window(1, np + 1, np + 1)
    # Apply window to the padded data.
    padWin     = (sin( arange(1, C.NUM_PAD + 1, dtype='float')
        / (C.NUM_PAD + 1) * C.PI / 2) **2)
    # Apply window to the trigger data.
    trgWin     = com.hann_window(1, 2 * C.NUM_PAD + 2, 2 * C.NUM_PAD + 2)
    slowWin    = ones(numBreaks - 1)
    
    if C.USE_TAPER:
        slowWin = com.hann_window(1, numBreaks - 1, numBreaks)
    
    # The (residual) carrier phase of each range bin.
    carrierPhase = exp( -1j * (4 * C.PI * fc / C.C * delta_r) * arange(np))

    if C.ANGLE_LIMIT < 90:
        for xIdx in range(imgAz.size):
            # determine which pixels are outside the angle limit of the image
            clipPixels = ( (abs(imgRg) < abs(imgAz[xIdx] 
               * tan(C.PI / 2 - C.ANGLE_LIMIT * (C.PI / 180)))))
            # Set out-of-bounds pixels to NaN.
            myImg[clipPixels, xIdx] = float('nan')

    # loop over break & process pulses.
    print("Processing position data...\n")
    for bIdx in range(2, numBreaks + 1):
        posIdx = bIdx - 2
        # All pulses for this position.
        myPulses = where(breakNumber == bIdx)[0]
        
        # Compute the zero-doppler mixer output for this position.
        tmpRP = zeros(np)
        for pIdx in range(1, minpulsesPerBreak-1):
            # Bandlimited interpolate the trigger signal.
            tmp        = ( trig[pulseStarts[myPulses[pIdx]]
                         + arange(-C.NUM_PAD, C.NUM_PAD + 1)]
                         * trgWin)
            interpTmp  = com.fft_interp(tmp, int(C.OSAMP_FACT))
            interpTmp  = (interpTmp[int((C.NUM_PAD*C.OSAMP_FACT))
                         + arange(-C.OSAMP_FACT, C.OSAMP_FACT + 1, dtype=int)])
            interpOffs = arange(-C.OSAMP_FACT, C.OSAMP_FACT + 1) / C.OSAMP_FACT
            myIdx      = where(2 == diff(sign(interpTmp)))[0] + 1
            if not myIdx.nonzero():
                print("Trigger edge not found... skipping pulse %d " % pIdx)
                continue
            tmp2       = interpTmp[myIdx + arange(-1, 1) ]
            
            # Linear interpolation to find zero crossings.
            fracOffset = (-(interpOffs[myIdx] - tmp2[1] / (tmp2[1]-tmp2[0])
                          / C.OSAMP_FACT))
            
            # Time-align the data to the trigger the zero crossing.
            cInds            = (pulseStarts[myPulses[pIdx]]
                               + arange(-C.NUM_PAD, np + C.NUM_PAD))
            tmp              = s[cInds.astype(int)]
            tmp              = tmp.astype(float)
            tmp[0:C.NUM_PAD] = tmp[0:C.NUM_PAD] * padWin
       
            tmp[len(tmp)-1:len(tmp) - 1-C.NUM_PAD:-1] = \
                tmp[(len(tmp) - 1):len(tmp) - 1-C.NUM_PAD:-1] * padWin
  
            # Time delay applied in the frequency domain below.
            tmp = fft.fft(tmp)
            tmp = (tmp * exp(-1j * arange(0, np + 2 * C.NUM_PAD)
                  / (np + 2 * C.NUM_PAD) * 2 * C.PI * fracOffset))
            tmp = fft.ifft(tmp)

            # Incorporate this pulse in the average (computing the zero-doppler mixer output).
            tmpRP = tmpRP + tmp[C.NUM_PAD + arange(0, np, dtype=int)]

        # Compute the range profile from the mixer output.
        # Apply fast & slow-time windows, then ifft.
        tmpRP = fft.ifft(tmpRP * (rngWin * slowWin[posIdx]))
        # Baseband remove carrier phase, interpolate up, and scale signal vs range.
        tmpRP = com.fft_interp(tmpRP * carrierPhase, int(C.OSAMP_FACT)) * rangeScale         
        # Compute the first difference in range (used for linear interpolation).
        diffRP = diff(append(tmpRP, 0), axis=0)

        # Incorporate this position into the image via backprojection.
        for xIdx in range(0, imgAz.size):
            # Compute the range to each image pixel & the matched filter terms.
            rangeVec = (sqrt(((imgAz[xIdx] - xa[posIdx])**2 + za[posIdx]**2)
                       + (imgRg - ya[posIdx])**2))
            matchVec = exp(1j * (4 * C.PI * fc / C.C) * rangeVec)
            
            # Compute integer and fractional range indices to each image pixel
            # for linear interpolation.
            rangeInds = rangeVec * (C.OSAMP_FACT / delta_r) + 1
            rangeLo   = floor(rangeInds)
            rangeLo   = rangeLo.astype(int)
            rangeFrac = rangeInds - rangeLo
            
            # Perform linear interpolation and apply the matched filter terms.
            # (Backprojection integral)
            myImg[:, xIdx] = (myImg[:, xIdx] + (tmpRP[rangeLo] 
                             + diffRP[rangeLo] * rangeFrac) * matchVec)
                             
                             
    # Make figure for output display.
    img_dB = 20*log10(abs(myImg))
    # Convert NaN to 0 for plotting purposes.
    img_dB[isnan(img_dB)] = 0.0
    img_dB = flipud(img_dB)
    plt.figure()
    plt.imshow(img_dB, cmap='jet', extent=(imgAz[0], imgAz[imgAz.size-1],
               imgRg[0], imgRg[imgRg.size - 1]), aspect='auto')

    plt.xlabel('X (meters)')
    plt.ylabel('Y (meters)')
    dynamicRng = array([-50.0, 10.0])
    plt.clim(img_dB.max() + dynamicRng)
    plt.colorbar()
    plt.show()   


wavFile = 'towardswarehouse.wav'
cantenna_sar(wavFile)
print("Processing of %s completed." % wavFile)
