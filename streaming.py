"""
 Class module used by cantenna RTI and Doppler streaming mode.
Pythonized code for Raytheon Low Cost Radar Initiative
Adapted from:
MIT Short Program 2013 Laptop Radar Course
(c) 2014 Massachusetts Institute of Technology
 Python version 1, 4/12/2014
"""

import constants as C
import common as common
from scipy import fftpack
import scipy.signal as sig
from numpy import append, arange, argwhere, array, ceil, concatenate, \
    convolve, cumsum, delete, diff, exp, floor, fromstring, hstack, isnan, linspace, log10, \
    median, nonzero, ones, shape, sign, sin, squeeze, tile, zeros
import pyaudio as pa
import pylab as pl
import warnings as warn

class dopProcess:
    """ Class supporting functions for doppler processing
    """
    def __init__(self):
        self.lastProcessedInd = 0
        self.freq             = []
        self.fThresh          = 0.0

    def get_thresh(self):
        if max(self.freq) < 350:
            self.fThresh = max(self.freq + 4)
        else:
            self.fThresh = 350


    def process_doppler(self, x):
        self.freq = abs(fftpack.fftshift(fftpack.fft(x, C.NFFT)))

    def update(self, parLen):
        self.lastProcessedInd = parLen

    def update_freq(self, parPosFreqIdx):
        self.freq                         = self.freq[parPosFreqIdx]
        self.freq[self.freq<self.fThresh] = 0

class RTIProcess:
    """
    Class for RTI streaming
    """
    def __init__(self):
        self.clean         = 0
        self.data          = 0
        self.dataRange     = 0
        self.fall          = 0
        self.groupNum      = 0
        self.lastProcessed = 0
        self.maxRange      = 100.0
        self.mti           = []
        self.mtidB         = []
        self.mtiCFARdB     = []
        self.nRange        = 0
        self.numPulses     = 0
        self.nP            = float('Inf')
        self.nRange        = 0
        self.nRangeKeep    = 0
        self.padWin        = 0
        self.pulseCount    = 0
        self.rise          = 0
        self.rngWin        = 0
        self.sif           = []#zeros(shape=[C.PULSE_BUFFER, 0], dtype='float32')
        self.tP            = 0
        self.trgWin        = 0
        self.trig          = 0
        self.h             = 0

        self.heightThresh       = .5
        # Minimum number of pulses per group.
        self.nMinPulsesPerGroup = 1
        # Drop this many pulses at the beginning of a group.
        self.nDiscardFirst      = 0
        # Drop this many pulses at the end of a group.
        self.nDiscardLast       = 0

        # Minimum duration of gap between pulse groups.
        tMinBreak          = 0.5
        # Minimum pulse duration.
        tMinPulse          = 15e-3
        # Maximum pulse duration.
        tMaxPulse          = 25e-3
        self.pulsewidthLThresh  = floor(tMinPulse*C.RATE)
        self.pulsewidthUThresh  = floor(tMaxPulse*C.RATE)
        self.breakwidthThresh   = floor(tMinBreak*C.RATE)

        # Initizlize mti filter.
        self.filter         = -ones((C.NPULSE_CANCEL, 1)) / C.NPULSE_CANCEL
        midIdx              = int(round((C.NPULSE_CANCEL) / 2))
        self.filter[midIdx] = self.filter[midIdx] + 1

    def mti_CFAR(self):
        if 0 == len(self.mti):
            #mti = zeros([self.numPulses, self.nRangeKeep], dtype='complex')
            self.mti = zeros([C.PULSE_BUFFER, self.nRangeKeep], dtype='complex')
        if 0 == len(self.mtidB):
            self.mtidB = zeros([C.PULSE_BUFFER, self.nRangeKeep], dtype='float')
        if 0 == len(self.mtiCFARdB):
            self.mtiCFARdB = zeros([C.PULSE_BUFFER, self.nRangeKeep], dtype='float')

         # Apply N pulse canceller.
        if self.pulseCount < C.PULSE_BUFFER:
            # 'temp' used to replicate 'same' option behavior in Matlab.
            # Code below is equivalent to conv(sif, mit_filter, 'same').
            temp = convolve(self.sif[range(self.pulseCount),:], self.filter)
            #print "any zeros in temp?= "
            #print (0 == temp).all()
            #print "size of new temp:"
            #print temp[range(1, temp.shape[0]),:].shape
            self.mti[range(self.pulseCount),:] = temp[range(1, temp.shape[0]), :]
            #print "any zeros in mti segment?"
            #print (0 == self.mti[range(self.pulseCount),:]).all()
            self.mtidB[range(self.pulseCount),:] = 20*log10(abs(self.mti[range(self.pulseCount),:]))
            #print "mtidB isinf="
            #print isinf(self.mtidB[range(self.pulseCount),:]).any()
            # Over time.
            print("pulse count: %d" % self.pulseCount)
            med       = tile(median(self.mtidB[0:self.pulseCount,:],0), \
                self.pulseCount).reshape(self.pulseCount, self.nRangeKeep)
            #print "1) size and sum of med:"
            #print med.shape
            #print sum(med)
            #print "1) size and sum of mtidB:"
            #print self.mtidB[0:self.pulseCount].shape
            #print sum(self.mtidB[0:self.pulseCount])
            #print isinf(med).any()
            #print isinf(self.mtidB[0:self.pulseCount]).any()
            test = self.mtidB[0:self.pulseCount] - med

           # print "size of test:"
            #print test.shape
            #print "2) Size of mtidB"
            #print self.mtidB[self.pulseCount:,:].shape
            self.mtiCFARdB = vstack([test,      self.mtidB[self.pulseCount:,:]])
            #print "3) Size of mtidB"
            #print self.mtidB[self.pulseCount:,:].shape
            # Over range.
            med       = transpose(tile(median( \
                self.mtiCFARdB[0:self.pulseCount, :], 1), [self.nRangeKeep, 1]))
            #print "2) size of med:"
            #print med.shape
            self.mtiCFARdB = vstack((self.mtiCFARdB[0:self.pulseCount] - med, \
                self.mtidB[self.pulseCount:,:]))
           # print "4) size of mtiCFARdB:"
          #  print self.mtiCFARdB.shape

        else:
            print("in else clause")
            # temp used since 'same' option behavior differs between Python and Matlab.
            # Code below is equivalent to conv(sif, mit_filter, 'same') in Matlab.
            temp       = convolve(self.sif, self.filter)
            self.mti   = temp[range(1,temp.shape[0]), :]
            self.mtidB = 20*log10(abs(self.mti))

            #med        = tile(median(self.mtidB, 0), self.numPulses).reshape(self.numPulses, self.nRangeKeep)
            med        = tile(median(self.mtidB, 0), self.mtidB.shape[0])
            #print "size of med"
            #print med.shape
            #print "size of mtiCFARdB:"
            #print self.mtiCFARdB.shape
            self.mtiCFARdB  = self.mtidB - med
            #med        = tile(median(self.mtiCFARdB, 1), self.numPulses).reshape(self.numPulses, self.nRangeKeep)
            med        = tile(median(self.mtiCFARdB, 1), self.mtiCFARdB.shape[0])
            self.mtiCFARdB  = self.mtiCFARdB - med



    def nearest_interp(self, xi, x, y):
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

    def process_sync(self):
        # Identify starts and stops of pulses and pulse groups.
        # Find when signal has been above threshold for some duration.
        x1 = self.clean >= self.heightThresh
        x1 = x1.astype(int)
        pl.figure()
        pl.plot(self.clean)
        #print "x1 max = %f " % x1.max()
        pl.figure()
        pl.plot(x1)
        x2 = convolve(x1, ones(floor(self.pulsewidthLThresh)), 'same')
        #print "x2 max = %f " % x2.max()
        pl.figure()
        pl.plot(x2)
        #print self.pulsewidthLThresh
        x3 = (x2 >= floor(self.pulsewidthLThresh).astype(int))
        #print "x3 max = %f" % x3.max()
        pl.plot(x3)
        x4 = append([0], diff(x3))

        # Initial rise and fall estimates.
        rise0 = nonzero(x4 == 1)[0]
        fall0 = nonzero(x4 == -1)[0]

        # Zero crossings of original signal.
        x5 = self.clean > 0

        # This assumes if the first sample is high, it's a rise, and if the
        # last sample is high, then the following sample is a fall
        x6 = diff(concatenate(([0], x5, [0]), axis=0))

        # Candidate rise and fall times (zero crossings).
        riseCandidates = nonzero(x6 == 1)[0]
        fallCandidates = nonzero(x6 == -1)[0]

        assert(all(fallCandidates > riseCandidates))

        #print x5.sum()

        # For each rise and fall estimate, find nearest candidate.
        rise = self.nearest_interp(rise0, riseCandidates, riseCandidates)
        fall = self.nearest_interp(fall0, fallCandidates, fallCandidates)
        #print rise
        #print fall
        assert(~any(isnan(rise)) | ~any(isnan(fall)))

        # This can happen with small positive signals because rise0 is based on
        # (non-zero) threshold crossings, and riseCandidates is based on zero
        # crossings
        indDiscard = []
        #print squeeze(argwhere(rise>rise0))
        #print squeeze(argwhere(fall < fall0))
        if 0 < len(squeeze(argwhere(rise > rise0))):
            if 0 < len(squeeze(argwhere(fall < fall0))):
                indDiscard = concatenate(squeeze(argwhere(rise > rise0)), squeeze(argwhere(fall < fall0)))
            else:
                indDiscard = squeeze(argwhere(rise > rise0))
        elif 0 < len(squeeze(argwhere(fall < fall0))):
             indDiscard = squeeze(argwhere(fall < fall0))

        #indDiscard = argwhere(rise > rise0) | argwhere(fall < fall0)
        #indDiscard = nonzero(indDiscard)[0]
        #print "indDiscard:"
     #   print indDiscard.size
        rise       = delete(rise, indDiscard)
        fall       = delete(fall, indDiscard)

        assert(all(fall>rise))

        if C.VERBOSE:
            print('Found %d raw pulses.' % rise.size)

        # Filter on pulsewidth.
        delta      = fall - rise
        indDiscard = (delta < self.pulsewidthLThresh) | (delta > self.pulsewidthUThresh)
        indDiscard = nonzero(indDiscard)[0]
        rise       = delete(rise, indDiscard)
        fall       = delete(fall, indDiscard)
        if C.VERBOSE:
            print('Discarded  %d  pulses based on pulsewidth.' % nonzero(indDiscard)[0].size)

        # Find group for each pulse.
        if rise.size > 1:
            diff_temp         = (diff(rise) >= floor(self.breakwidthThresh)).astype(int)
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
        indDiscardGroup = pulsesPerGroup < self.nMinPulsesPerGroup
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
            print('Discarded %d pulses based on pulses per group' % indDiscard.size)

        # Find which group each pulse belongs.
        groupNum = cumsum(isFirstOfGroup)
        #print "groupNun:\n"
        #print groupNum
        # Toss first and last n pulses in each group.
        #print "temp:\n"
        temp = ones(self.nDiscardFirst + self.nDiscardLast)
        #print temp
        if 0 == len(temp):
            temp = array([0])

        indDiscard = convolve(isFirstOfGroup, temp)
        # First nDiscardLast pts.
        indDiscard = delete(indDiscard, range(0, self.nDiscardLast))

        # Last nDiscardFirst-1 pts.
        indDiscard = delete(indDiscard, range(indDiscard.size - self.nDiscardFirst, indDiscard.size))

        indDiscard = indDiscard.nonzero()[0]

        self.rise       = delete(rise, indDiscard)
        self.fall       = delete(fall, indDiscard)
        self.groupNum   = delete(groupNum, indDiscard)

        if C.VERBOSE:
            print('Discarded  %d pulses at beginnings and ends of groups.' % indDiscard.size)
            if 0 < len(self.rise):
                print('Final answer: %d pulses in  %d groups.' %  (rise.size, self.groupNum[self.groupNum.size-1]))

    def set_data(self, parTrig, parData):
        self.trig = -parTrig
        self.data = -parData

    def set_last_processed(self, x):
        self.lastProcessed = x

    def update_params(self):
        # Refine using measured parameters: the pulse width.
        if 0 < len(self.rise):
            nPTmp = round(min(diff(self.rise)) / 2)
        else:
            nPTmp = 0
        tPTmp = nPTmp / float(C.RATE)
        if C.VERBOSE:
            print('Measured pulsewidth of %f ms.' % (tPTmp * 1000.0))
        if nPTmp < self.nP:
            self.nP = nPTmp
            self.tP = tPTmp
            # Number of ouput range samples
            self.nRange     = int(floor(C.OSAMP_RNG * self.nP / 2))
            # Range of each bin.
            self.dataRange  = arange(0.0, self.nRange, 1.0) * (C.DELTA_RG / C.OSAMP_RNG)
            # Keep this for plotting purposes.
            dataRgIdx       = self.dataRange > self.maxRange
            self.dataRange  = delete(self.dataRange, argwhere(dataRgIdx))
            # Number of range bins to keep.
            self.nRangeKeep = len(self.dataRange)
            # The window applied to reduce range sidelobes.
            self.rngWin     = common.hann_window(1, self.nP + 1, self.nP + 1)
            # The window applied to the padded data.
            self.padWin     = sin(array(range(1, C.NUM_PAD + 1), 'float') / \
                (C.NUM_PAD + 1) * C.PI/2)**2
            # The window applied to the trigger data.
            self.trgWin = common.hann_window(1, C.NUM_PAD*2 + 2, C.NUM_PAD*2 + 2)

        nSamples         = len(self.data)
        self.pulseStarts = self.rise[self.rise + self.nP + C.NUM_PAD <= nSamples]
        self.pulseStarts = self.pulseStarts[self.pulseStarts - C.NUM_PAD > 0]
        self.numPulses   = self.pulseStarts.size
        if C.VERBOSE:
            print('Found %d pulses.' % self.numPulses)
        if 0 == self.numPulses:
            warn.warn('No pulses were found. Check that radar is on and sync pulse is enabled.')
        self.pulseCount  = self.pulseCount + self.numPulses

    def update_sif(self):
        if 0 < self.nRangeKeep:
            if 0 == len(self.sif):
                self.sif = zeros(shape=[C.PULSE_BUFFER, self.nRangeKeep], dtype='complex')
            else:
                # Trim column length to new minimum number of samples.
                self.sif = self.sif[:,range(self.nRangeKeep)]
            # Clear rows for new pulses.
            self.sif = vstack((zeros([self.numPulses,self.nRangeKeep]),  \
                    self.sif[range(self.sif.shape[0]-self.numPulses),:]))
        else:
            print('Error, no valid pulses.')
            return

        for pIdx in range(self.numPulses):
            # Bandlimited interpolate the trigger signal.
            tmp        = self.clean[self.pulseStarts[pIdx] + \
                arange(-C.NUM_PAD,C.NUM_PAD+1)] * self.trgWin
            interpTmp  = common.fft_interp(tmp, C.OSAMP_TRG)
            interpTmp  = interpTmp[(C.NUM_PAD*C.OSAMP_TRG ) + \
                arange(-C.OSAMP_TRG, C.OSAMP_TRG+1)]
            interpOffs = arange(-C.OSAMP_TRG, C.OSAMP_TRG + 1, dtype=float) /  \
                C.OSAMP_TRG
            myIdx      = int(squeeze(nonzero(diff(sign(interpTmp)) == 2))) + 1
            if 0 == myIdx:
                 print('Trigger edge not found... skipping pulse %d' % pIdx)
                 continue
            tmp2 = array([interpTmp[myIdx-1], interpTmp[myIdx]])
            # Linear interpolate to find the zero crossing.
            fracOffset = -(interpOffs[myIdx] - tmp2[1]/(tmp2[1]-tmp2[0]) / C.OSAMP_TRG)
            # Time-align the data to the trigger event (the zero crossing).
            cInds            = self.pulseStarts[pIdx] + arange(-C.NUM_PAD,(self.nP +C.NUM_PAD))
            tmp              = self.data[cInds.astype(int)]
            tmp              = tmp.astype(float)
            tmp[0:C.NUM_PAD] = tmp[0:C.NUM_PAD] * self.padWin
            tmp[arange(len(tmp) - 1, len(tmp) - C.NUM_PAD-1, -1)] =  \
                tmp[arange(len(tmp) - 1, len(tmp) - C.NUM_PAD - 1, -1)] * self.padWin
            # Time delay applied in the frequency domain below.
            tmp          = fftpack.fft(tmp)
            #print fracOffset
            #print sum(exp(-1j*(arange(0, self.nP + 2 * C.NUM_PAD)) / (self.nP + 2 * C.NUM_PAD) * 2 * C.PI * fracOffset))
            tmp          = tmp * exp( -1j*(arange(0, self.nP + 2 * C.NUM_PAD)) / \
                (self.nP + 2 * C.NUM_PAD) * 2 * C.PI * fracOffset)
            tmp          = fftpack.ifft(tmp)
            # Compute & scale range data from the time-aligned mixer output.
            tmp          = fftpack.ifft(tmp[C.NUM_PAD + arange(0, self.nP, dtype=int)] * self.rngWin, 2*self.nRange)
            self.sif[pIdx,:] = tmp[0:self.nRangeKeep]
        print("any zeros in sif?")
        print((self.sif[range(self.numPulses),:]==0).all()) #DEBUG

class TIPlot:
    """ Class supporting graphic display for doppler and RTI streaming.
    Attributes:
        ti: 2D array of floats containing the doppler returns.
        posFreqInd: Index of positive doppler frequency bins.

    """
    def __init__(self, parType, parXlbl, parYlbl):
        """ Initializes TIPlot class with elements required for plotting.
            parType: String
                'Doppler' or 'RTI'
            parXlbl: String
                The plot's xlabel.
            parYlbl: String
                The plot's ylabel.
        """
        if 'DOPPLER' == parType.upper():
            self.posFreqInd = arange(ceil(C.NFFT/2), C.NFFT, dtype='int32')
            self.ti         = zeros(shape=[C.NUM_ROWS, len(self.posFreqInd)],
                                       dtype='float32')
            F               = linspace(0, C.RATE/2, len(self.posFreqInd))
            self.velAxisMPS = F*(C.C/C.FC)/2.0
            self.velAxisMPH = self.velAxisMPS*.447

            self.fig, self.ax = pl.subplots(1,1)

            self.image = self.ax.imshow(20*(self.ti), #20*log10(self.ti),
               cmap='jet',
               extent=[0,self.velAxisMPS[len(self.velAxisMPS)-1], C.NUM_ROWS,0],
               aspect='auto',
               vmin=0,
               vmax=.1)

            pl.gca().invert_yaxis()
            self.fig.canvas.draw()
            pl.xlim(1, 20)
            pl.xlabel(parXlbl)
            pl.ylabel(parYlbl)
            self.fig.colorbar(self.image)
            pl.show()
            if True == C.USE_MPH:
                self.ax.set_title('Current Speed: %2.2f MPH' % 0.0)
            else:
                self.ax.set_title('Current Speed: %2.2f m/s' % 0.0)
        else: # RTI
            self.fig, self.ax = pl.subplots(1,1)
            self.image        = self.ax.imshow([[1],[1]],
                                           cmap='jet',
                                           aspect='auto',
                                           extent=[0.5, 1.5, 0.5, 1.5],
                                           vmin=0,
                                           vmax=37)
            self.ax.set_title("RTI with MTI+CFAR")
            pl.xlabel(parXlbl)
            pl.ylabel(parYlbl)
            self.fig.colorbar(self.image)
            pl.show()

    def RTI_plot(self, parTitle, parXlabel, parYlabel):
        """ Set image parameters for RTI plot
        """
#        h = imagesc(nan);
#        set(gca,'YDir','normal');
#        colormap(jet(256));
        #caxis([0 40]-3);
        self.fig, self.ax = pl.subplots(1,1)
        self.image        = self.ax.imshow(20*log10(self.ti),
                                           cmap='jet',
                                           aspect='auto')
        self.ax.set_title(parTitle)
        pl.xlabel(parXlabel)
        pl.ylabel(parYlabel)
        self.fig.colorbar(self.image)

    def doppler_plot(self, parXlabel, parYlabel):
        """ Set image parameters for Doppler plot
            fig: Figure handle object.
            ax: Axes object of fig.
            F:
            velAxis: Velocity axis for fig
        """
        F            = linspace(0, C.RATE/2, len(self.posFreqInd))
        self.velAxis = F*(C.C/C.FC)/2
        if True == C.USE_MPH:
            self.velAxis = self.velAxis*.447
            self.ax.set_title('Current Speed: %2.2f MPH' % 0.0)
        else:
            self.ax.set_title('Current Speed: %2.2f m/s' % 0.0)
        self.fig, self.ax = pl.subplots(1,1)
        self.image = self.ax.imshow(20*log10(self.ti),
               cmap='jet',
               extent=[0,self.velAxis[len(self.velAxis)-1], C.NUM_ROWS,0],
               aspect='auto',
               vmin=0,
               vmax=0)
        pl.gca().invert_yaxis()
        self.fig.canvas.draw()
        pl.xlim(1, 20)
        pl.xlabel(parXlabel)
        pl.ylabel(parYlabel)
        self.fig.colorbar(self.image)

    def update_doppler(self, newData):
        """ Updates ti graph with new row of data.
        """
        maxInd = newData.argmax()
        if True == C.USE_MPH:
            self.ax.set_title('Current Speed: %2.2f MPH' % self.velAxisMPH[maxInd])
        else:
            self.ax.set_title('Current Speed: %2.2f m/s' % self.velAxisMPS[maxInd])

        self.ti[1:C.NUM_ROWS,:] = self.ti[0:C.NUM_ROWS-1,:]
        newData                     = 20*log10(newData)
        newData[newData < 0]        = 0
        self.ti[0,:]                = newData
        self.image.set_data(self.ti)
        self.image.set_clim(vmin=0, vmax=self.ti.max())

    def update_rti(self, objR):
        self.yData = array([range(0, objR.sif.shape[0])], dtype='float')*objR.tP*2
        self.xData = array([range(0, objR.dataRange.shape[0])], dtype='float')

     #   pl.xlim(objR.dataRange[0], objR.dataRange.max())
       # pl.ylim(0, 20)
        self.image.set_data(objR.mtiCFARdB)
        self.image.set_extent([0.0, 100, \
        0.0, 20])
       # self.image.set_aspect = 'auto'
        #self.image.set_extent=(0,100, 20,0)
        #self.ax.relim()
        #self.ax.autoscale_view(True,True,True)
        #self.image.set(h,'XData',xData);
        #set(h,'YData',yData)
#        axis(get(h,'Parent'),[min(xData) max(xData) min(yData) max(yData)]);
#

class Stream:
    """ Class that manages pyaudio input stream from cantenna.
        Input recording device is set to iMic if present.

        Attributes:
             p: Pyaudio object. Opens and closes recording stream.
             lastProcessed: Integer. Index holding last processed data in audio
             stream.
             recordLen: Float. Duration in seconds for data capture.
             data: formatted audio data
             raw: Raw input from stream.
             sync: sync data for DTI mode.
             devInd: Integer. Index identifying input audio device.
    """
    def __init__(self, parRecordLen=0, parFormat=0):
        """
		    Initializes Stream object.
        """
        self.p             = pa.PyAudio()
        self.lastProcessed = 0
        self.recordLen     = parRecordLen
        self.data          = zeros(shape([int(C.RATE * self.recordLen), 2])
                               , dtype='float32')
        self.raw           = self.data
        self.sync          = self.data
        self.devInd        = []
        self.get_device()
        if 'float32' == parFormat:
            self.FORMAT = pa.paFloat32
        elif 'int16' == parFormat:
            self.FORMAT = self.p.paInt16

    def __del__(self):
        self.p.terminate()

    def format_raw(self):
        """
            Formats raw stream data to floats.
        """
        self.raw  = fromstring(self.raw, 'float32')
        self.raw  = self.raw.reshape( int(C.RATE * self.recordLen ),2)
        self.data = self.raw[:, 1]
        # Sync pulse not relevant for Doppler mode
        self.sync = self.raw[:, 0]

    def get_data(self):
        """
            Reads from stream. Note Matlab's audio port input stream seems to
            be three orders of magintude greater than than of PyAudio. Does not
            appear to adversely affect results; Schmitt Trigger simply becomes
            OBE.
        """
        self.raw = self.stream.read(int(C.RATE * self.recordLen))


    def get_device(self):
        """
            Sets recording device to iMic if present.
        """
        x           = self.p.get_default_input_device_info()
        self.devInd = x['index']
        for ii in range(self.p.get_device_count()):
            x = self.p.get_device_info_by_index(ii)
            if 'iMic ' in x['name']:
                self.devInd = ii
                print('Audio device set to %s' % x['name'])
                break

    def open_stream(self):
        """
            Opens stream and saves raw data.
        """
        self.stream = self.p.open(format=self.FORMAT,
                                 channels=C.NUM_CHANNELS,
                                 rate=C.RATE,
                                 input=True,
                                 frames_per_buffer=C.CHUNK)

    def set_raw_num_samp(self):
        self.numRawSamp = len(self.raw)
        print('last processed: %f' % self.lastProcessed)
        self.raw = self.raw[self.lastProcessed:]
        self.lastProcessed = self.numRawSamp

    def stop_stream(self):
        self.stream.stop_stream()
        self.stream.close()

