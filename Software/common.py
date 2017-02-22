"""
Common function module used by cantenna RTI and Doppler files
Pythonized code for Raytheon Low Cost Radar Initiative
Adapted from:
MIT Short Program 2013 Laptop Radar Course
(c) 2014 Massachusetts Institute of Technology
Version 1 6/1/2014
Version 2 7.1.2016
    Add chan_check function.
"""

import constants as C
from numpy import (append, arange, argwhere, array, concatenate, convolve, cos,
                   cumsum, delete, diff, fft, floor, fromstring, hstack, isnan,
                   nonzero, ones, setdiff1d, sin, tile, zeros)
import sys
sys.path.append('.')
import wave


def chan_check(raw2d, caller):
    """
    Kits may have audio cable wired differently. Confirm which channel has the 
    data. Both channels have variable signal levels due to tuning, but beat 
    frequency channel should have more values closer to zero in SAR mode.
    If RTI, beat channel should have no zero crossings. Doppler mode uses only
    one channel.
    Inputs: Raw data from .wav file
    Outputs: data, beat
    """
    ch0 = raw2d[:,0]
    ch1 = raw2d[:,1]
    if 'SAR' == caller:
        lowVals0 = argwhere(abs(ch0) < 20).size
        lowVals1 = argwhere(abs(ch1) < 20).size
        if lowVals0 < lowVals1:
            data = ch0
            beat = ch1
        else:
            data = ch1
            beat = ch0 
    elif 'RTI' == caller:
        numZeros0 = argwhere(ch0 == 0).size
        numZeros1 = argwhere(ch1 == 0).size
        if numZeros0 < numZeros1:
            beat = ch0
            data = ch1
        else:
            data = ch0
            beat = ch1
    return(data, beat)

    
def fft_interp(x, M):
    """
    Perform approximate bandlimited interpolation of x by a factor of M.
    Inputs: x, M
    """
    L       = 4
    mArray  = arange(-L*M, L*M+1, dtype=float)
    winInds = mArray/M * C.PI

    # Get the ideal anti aliasing filter's impulse response of length 2*M + 1.
    winInds[L *M] = 1
    myWin         = sin(winInds) / winInds
    myWin[L * M]  = 1

    # Use the window method; apply a hann window.
    myWin = myWin * hann_window(1, 2*L*M + 2, 2*L*M + 2)

    # Insert zeros in data and apply antialias filter via FFT.
    nFFT = x.size * M
    fftx = fft.fft(x)

    fftx = tile(fftx, M)
    # Reality check.
    if 'cfloat' != x.dtype:
        invY = fft.ifft(fft.fft(myWin, nFFT) * fftx)
        invY = invY.real
    else:
        invY = fft.ifft(fft.fft(myWin, nFFT) * fftx)

    y = concatenate( (invY[(L*M):], invY[0:(L*M)]), axis=0)
    return y


def hann_window(start, stop, denominator):
    """
    Create a Hann (cosine squared) window.
    Doppler and RTI ML code used distinct implementations of the Hann function. The 
    three input arguments allow each version to be run as they originally appeared. 
    Inputs: 
      start             Start of window interval
      stop              End of window interval
      denominator       Full width at half maximum
    Outputs: Hann window as float type.
    """
    nArray = arange(start, stop, 1, dtype=float)
    return .5 + .5 * cos(C.TWOPI * (nArray/denominator - .5))


def nearest_interp(xi, x, y):
    """ Perform nearest neighbor interpolation/extrapolation.
    Inputs: 
        xi  Vector to interpolate/extrapolate
        x   Domain vector
        y   Range vector
    Outputs:
       Range vector associated with xi
     """
    idx = abs(x - xi[:,None])
    return y[idx.argmin(axis=1)]


def schmitt_trigger(x, threshold):
    """
    Schmitt trigger implementation.
    Inputs:
        x    Vector of any numeric type
    Outputs:
        cleantrig Denoise vector based on C.SCHMITT_TRESH threshold
    """
    cleanTrig = array(x, copy=True)
    gtTrig    = (threshold < x)
    ltTrig    = (-threshold> x)
    setFlag   = nonzero(gtTrig + ltTrig)[0]

    cleanTrig[nonzero(gtTrig)] =  1
    cleanTrig[nonzero(ltTrig)] = -1
    lockedIdx                  = zeros(len(x))
    if 0 == len(setFlag):
        return x
    lockedIdx[setFlag[0]:]     = 1
    locked                     = nonzero(lockedIdx)[0]
    unlocked                   = (0 == lockedIdx).nonzero()
    cleanTrig[unlocked]        = 0

    valsToChange = setdiff1d(locked, setFlag)
    for ii in valsToChange:
        cleanTrig[ii] = cleanTrig[ii -1]
    return cleanTrig


def process_sync(x, *args):
    # Identify starts and stops of pulses and pulse groups.
    if len(args) < 1:
        print("process_sync using default parameters.")
        samplesPerSec      = C.RATE
        heightThresh       = .5
        # Minimum pulse duration.
        tMinPulse          = 15e-3
        # Maximum pulse duration.
        tMaxPulse          = 25e-3
        # Minimum duration of gap between pulse groups.
        tMinBreak          = 0.5
        # Minimum number of pulses per group.
        nMinPulsesPerGroup = 25
        # Drop this many pulses at the beginning of a group.
        nDiscardFirst      = 3
        # Drop this many pulses at the end of a group.
        nDiscardLast       = 3
    else:
    # Get parameters from arguments.
        param = args[0]
        samplesPerSec      = param.samples_per_sec
        heightThresh       = param.height_thresh
        tMinPulse          = param.t_min_pulse
        tMaxPulse          = param.t_max_pulse
        tMinBreak          = param.t_min_break
        nMinPulsesPerGroup = param.nMinPulsesPerGroup
        nDiscardFirst      = param.nDiscardFirst
        nDiscardLast       = param.nDiscardLast
    
    pulsewidthThresh = floor(tMinPulse*samplesPerSec)
    breakwidthThresh = floor(tMinBreak*samplesPerSec)
    
    # Find when signal has been above threshold for some duration.
    x1 = x >= heightThresh
    x1 = x1.astype(int)
    x2 = convolve(x1, ones(pulsewidthThresh), 'same')

    x3 = (x2 >= pulsewidthThresh).astype(int)
    x4 = append([0], diff(x3))
    
    # Initial rise and fall estimates.
    rise0 = nonzero(x4 == 1)[0]
    fall0 = nonzero(x4 == -1)[0]

    # Zero crossings of original signal.
    x5 = x > 0
    
    # This assumes if the first sample is high, it's a rise, and if the
    # last sample is high, then the following sample is a fall
    x6 = diff(concatenate(([0], x5, [0]), axis=0))

    # Candidate rise and fall times (zero crossings).
    riseCandidates = nonzero(x6 == 1)[0]
    fallCandidates = nonzero(x6 == -1)[0]
    
    assert(all(fallCandidates > riseCandidates))
    
    # For each rise and fall estimate, find nearest candidate.
    rise = nearest_interp(rise0, riseCandidates, riseCandidates)
    fall = nearest_interp(fall0, fallCandidates, fallCandidates)

    assert(~any(isnan(rise) | isnan(fall)))
    
    # This can happen with small positive signals because rise0 is based on
    # (non-zero) threshold crossings, and riseCandiates is based on zero
    # crossings
    indDiscard = argwhere(rise > rise0) | argwhere(fall < fall0)
    indDiscard = nonzero(indDiscard)[0]

    rise       = delete(rise, indDiscard)
    fall       = delete(fall, indDiscard)
    
    assert(all(fall>rise))
    
    if C.VERBOSE:
        print("Found %d raw pulses." % rise.size)
    
    # Filter on pulsewidth.
    delta      = fall - rise
    indDiscard = (delta < tMinPulse*samplesPerSec) | (delta > tMaxPulse*samplesPerSec)
    indDiscard = nonzero(indDiscard)[0]
    rise       = delete(rise, indDiscard)
    fall       = delete(fall, indDiscard)
    if C.VERBOSE:
        print("Discarded %d pulses based on pulsewidth." %  nonzero(indDiscard)[0].size)
    
    # Find group for each pulse.
    if rise.size > 1:
        diff_temp         = (diff(rise) >= breakwidthThresh).astype(int)
        isFirstOfGroup    = hstack([1, diff_temp])
    else: 
        isFirstOfGroup    = ones(rise.size)

    # Filter on number of pulses per group here.
    firstPulseInd = nonzero(isFirstOfGroup)[0]
    
    if firstPulseInd.size > 1:
        lastPulseInd  = append(firstPulseInd[1:] - 1, isFirstOfGroup.size)
    else:
        lastPulseInd = isFirstOfGroup.size
    
    pulsesPerGroup  = 1 + lastPulseInd - firstPulseInd
    indDiscardGroup = pulsesPerGroup < nMinPulsesPerGroup
    indDiscard      = zeros(rise.size, 'int')
    #TODO Inefficient, update later
    for k in range(pulsesPerGroup.size):
        if indDiscardGroup[k]:
            indDiscard[firstPulseInd[k]:lastPulseInd[k]+1] = True

    indDiscard     = nonzero(indDiscard)[0]
    rise           = delete(rise, indDiscard)
    fall           = delete(fall, indDiscard)
    isFirstOfGroup = delete(isFirstOfGroup, indDiscard)
    
    if C.VERBOSE:
        print("Discarded %d pulses based on pulses per group" % indDiscard.size)
    
    # Find which group each pulse belongs.
    groupNum = cumsum(isFirstOfGroup)

    # Toss first and last n pulses in each group.
    temp       = ones(nDiscardFirst + nDiscardLast)
   
    indDiscard = convolve(isFirstOfGroup, temp)
    # First nDiscardLast pts.
    indDiscard = delete(indDiscard, range(0, nDiscardLast))
    # Last nDiscardFirst-1 pts.
    indDiscard = delete(indDiscard, range(indDiscard.size-nDiscardFirst, indDiscard.size-1))
   
    indDiscard = indDiscard.nonzero()[0]
    rise       = delete(rise, indDiscard)
    fall       = delete(fall, indDiscard)
    groupNum   = delete(groupNum, indDiscard)
    
    if C.VERBOSE:
        print("Discarded  %d pulses at beginnings and ends of groups." % indDiscard.size)
    
    if C.VERBOSE:
        print("Final answer: %d pulses in  %d groups." % (rise.size, groupNum[groupNum.size-1]))

    return rise, fall, groupNum


def read_wavefile(waveFile):
    """Read a .wav file. Assumes int16 encoding.
    Inputs:
       waveFile, the wave file name as string.
    Outputs:
       the frame rate and decoded channels.
    """
    print("Loading WAV file...")
    
    #read and store waveFile onto fp
    fp         = wave.open(waveFile, "rb")
    
    #returns and stores the number of audio frames found in waveFile
    numFrames = fp._nframes
    
    #returns and stores the number of audio channels in waveFile
    numChanns = fp._nchannels

    # Get both channels interleaved as strings.
    # readframes(numFrames) returns numFrames as a string of bytes
    # fromstring returns a 1-D string array of bytes, interpreted as 16bit
    # integers => (int16)
    raw = fromstring(fp.readframes(numFrames), 'int16')
    fp.close()
    
    
    # Convert strings to Integers
    # zeros(params) takes in the array of numFrames,numChanns
    # and converts the return data type as 32 bit size integer array
    y   = zeros([numFrames,numChanns], dtype = 'int32')
    
    # for loop that iterates the length of numChanns
    # ii is a temp data type to store individual values from numChanns
    for ii in range(numChanns):
        # creates an array(y) that goes through every ii data value
        # raw[ii::2] is an array starting at the value ii and increases
        # ii by 2
        y[:,ii] = raw[ii::2]
        
    # returns an array of frame rates and y values
    return (fp.getframerate(), y)
