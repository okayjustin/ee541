#Python code to process the waveFile generated in Record_Data.py
import common as com
import constants as C
import scipy.signal as sig
import os.path
from matplotlib import pyplot as plt
from numpy import (append, arange, argwhere, array, delete, diff, exp, floor, fft
 log10, median, nonzero, ones, sign, sin, tile, zeros)
from scipy import fftpack
from time import time
import wave as w

#modify Fields according to what we need
CHUNK = 1024  #Number of Samples taken at individual times
FORMAT = pyaudio.paInt16  #data type extracted from pyaudio
CHANNELS = 1  #Numer of channels we are recording (1?)
RATE = 44100  #Sampling Rate 
RECORD_TIME = 20  #to modify for different record lengths
RPM = 5  #to modify according to actual speed of rotation

#short time fourier transform function
def stft(x, fftsize=CHUNK, overlap=4):   
    hop = fftsize / overlap
    w = scipy.hanning(fftsize+1)[:-1]      # better reconstruction with this trick +1)[:-1]  
    return np.array([np.fft.rfft(w*x[i:i+fftsize]) for i in range(0, len(x)-fftsize, hop)])
    

# another function for calculating short time Fourier transform
# method returns stft data 
# data = a numpy array containing the signal to be processed
# fs = a scalar which is the sampling frequency of the data
def stft2(data, fs, fft_size, overlap_fac):  
	hop_size = np.int32(np.floor(fft_size * (1-overlap_fac)))
	pad_end_size = fft_size          # the last segment can overlap the end of the data array by no more than one window size
	total_segments = np.int32(np.ceil(len(data) / np.float32(hop_size)))
	t_max = len(data) / np.float32(fs)
 
	window = np.hanning(fft_size)  # our half cosine window
	inner_pad = np.zeros(fft_size) # the zeros which will be used to double each segment size
 
	proc = np.concatenate((data, np.zeros(pad_end_size)))              # the data to process
	result = np.empty((total_segments, fft_size), dtype=np.float32)    # space to hold the result
 
	for i in xrange(total_segments):                      # for each segment
    	current_hop = hop_size * i                        # figure out the current segment offset
    	segment = proc[current_hop:current_hop+fft_size]  # get the current segment
    	windowed = segment * window                       # multiply by the half cosine function
    	padded = np.append(windowed, inner_pad)           # add 0s to double the length of the data
    	spectrum = np.fft.fft(padded) / fft_size          # take the Fourier Transform and scale by the number of samples
    	autopower = np.abs(spectrum * np.conj(spectrum))  # find the autopower spectrum
    	result[i, :] = autopower[:fft_size]               # append to the results array
 
	result = 20*np.log10(result)          # scale to db
	result = np.clip(result, -40, 200)    # clip values
	
	#create frequency list which correspond to bins in result
	freqs = np.arange(fft_size / np.float32(fft_size * 2) * fs
	
	return result, freqs

def polar_Image(wavFile):

	# Display if file is not found
   if not os.path.isfile(wavFile):
   	print("File %s not found." % (wavFile))
   
   #data contains the data of the wave, (magnitudes and frequencies
   [fs, data] = com.read_wavefile(wavFile)
   
   #make the list of data into an array
   data = numpy.array(data)
   
   #get the number of elements in the array
   npDLength = np.size(data, 0)
   
   #calculate the degree/element
   #this angle rate corresponds to which angle we will sample and do
   #fft on to extract frequencies and magnitudes
   DEG_PER_ELEM = (RPM/60 * RECORD_TIME * 360 / npDLength
   
   #creates sampling interval given CHUNK
   #each element in the array contains 1024 samples
   #k = np.aragne(npDLength)
   #T = npDLength/RATE
   
   #will use to plot frequency (frequency is proportional to distance)
   #frq = k/T  #double sided 
   #frq = frq[range(n/2)] #one sided frequency
   
   #stores Fourier data
   #Ydata = np.fft.fft(npData)/npDLength ##fft computing and normalizing
   #Ydata = Ydata[range(npDLength/2)]
   
   #get the Fourier frequencies
   #freqs = np.fft.fftfreq(npDLength)
   
   #use for loop to construct angles list
   #based on the DEG_PER_ELEM and data length
   #must have angles length = data length
   angles = []
   
   for i in range(1, npDLenght):
   	angles.append(i * DEG_PER_ELEM)
   
   #get the Fourier data 
   [result , freqs] = stft2(data, RATE, CHUNK, 0.5)
   
   
   
   
   
   
   
   
   	
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
    