"""
Produces a DTI (Doppler x time intensity) image of the
cantenna recording. Assumes the data is collected with
the radar in a continuous wave (CW) mode. See inside
the function for more parameters.
   wavFile = the filename of the .WAV file to process
Pythonized code for Raytheon Low Cost Radar Initiative
Adapted from:
MIT Short Program 2013 Laptop Radar Course
(c) 2014 Massachusetts Institute of Technology
Version 1 12/6/2013
Version 2 9/7/2014
    Added import of os.path
Version 3 8/11/2015
    Updated print statement to function for Python 3.4 compatibility
Version 4 6/29/2016
    Updated style to PEP8 standards.
"""

import common as com
import constants as C
import matplotlib.pyplot as plt
from numpy import (arange, array, floor, fft, flipud, log10, mean, round,
  where, zeros)
import os.path
import time

def cantenna_dop_v4(wavFile):
    # Use default wave file if wavFile does not exist.
    if not os.path.exists(wavFile):
        wavFile = 'data/dop_Stocker_LaBrea_SW_11Jan2014_bigRig.wav'
    [fs, y] = com.read_wavefile( wavFile)

    # Derived parameters.
    # Samples per pulse.
    ns = round(C.CPI * fs)

    # Note following negation was included in the ML implementation and is
    # retained here as a comment. Inverting the y vector does not impact DTI.
    # Reverse inverted input.
    x =  -1*y[:,0]

    # Grab an integer number of overlapped frames.
    m = floor(len(x) / ns * C.OLAP_FACT) - C.OLAP_FACT + 1
    m = m.astype(int)
    # Compute axes parameters for the plot.
    # The Doppler data are oversampled by OSAMP_DOP
    # Frequency in Hertz.
    delta_f = arange( C.OSAMP_DOP*ns/2-1 ) / (C.OSAMP_DOP*ns) * fs
    # Doppler -> speed.
    speed   = delta_f * C.LAMBDA / 2
    # Collection time (sec)
    time    = arange(1, m) * C.CPI / C.OLAP_FACT

    # Limit the speed axis to a reasonable range.
    speed = speed[where(speed <= C.MAX_SPEED)]
    speedLen = len(speed)
    # Compute Doppler window.
    dopWin = com.hann_window(0, ns, ns-1)

    # Compute the Doppler vs. time plot.
    print("Processing...")
    dti = zeros([speedLen, m])
    for mIdx in  range(m):
        # Compute indices.
        xInds       = arange(int(ns)) + int(mIdx*floor(ns/C.OLAP_FACT))
        # Apply Doppler window.
        tmp         = x[xInds] * dopWin
        # Remove DC component if it exists.
        tmp         = tmp - mean(tmp)
        # Compute oversampled Doppler spectrum.
        tmp         = fft.fft(tmp, int(C.OSAMP_DOP*ns))
        # Grab result in dB.
        dti[:,mIdx] = (abs(tmp[0:speedLen])) # 20*log10(abs(tmp[0:speedLen]))
    # Transpose data.
    dti = dti.T
    dti = flipud(dti)

    # Make Doppler vs. time plot.
    plt.figure()
    plt.imshow(dti, extent=[speed.min(), speed.max(), time.min(), \
        time.max()], cmap='jet', aspect='auto')
    vAx = plt.gca()
    vAx.invert_yaxis()
    dynamicRng = array([-60.0, 0.0])
    plt.clim(dti.max() + dynamicRng)
    plt.xlabel('Speed (m/sec)')
    plt.ylabel('Time (sec)')
    plt.colorbar()


# If specifying path to a data file, file separator is forward slash, even on
# Windows. For example, see commented line below.
wavFile = 'Z:/rasPI/data/dop_volume_up_Stocker_LaBrea_SW_11Jan2014_bigRig.wav'
#wavFile = 'dop_volume_up_Stocker_LaBrea_SW_11Jan2014_bigRig.wav'
# Start timer.
startTime = time.time()
cantenna_dop_v4(wavFile)
stopTime  = time.time()
print("Total time was ", stopTime - startTime, " seconds.")

