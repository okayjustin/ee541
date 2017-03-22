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

def polar_Image(wavFile):

	# Display if file is not found
   if not os.path.isfile(wavFile):
   	print("File %s not found." % (wavFile))
   
   #data contains the data of the wave, (magnitudes and frequencies
   [fs, data] = com.read_wavefile(wavFile)
   
   #make the list of data into an array
   npData = numpy.array(data)
   
   #get the number of elements in the array
   npDLength = np.size(npData, 0)
   
   #calculate the degree/element
   #this angle rate corresponds to which angle we will sample and do
   #fft on to extract frequencies and magnitudes
   DEG_PER_ELEM = (RPM/60 * RECORD_TIME * 360 / npDLength
   
   #creates sampling interval given CHUNK
   #each element in the array contains 1024 samples
   k = np.aragne(npDLength)
   T = npDLength/CHUNK
   
   #will use to plot frequency (frequency is proportional to distance)
   frq = k/T  #double sided 
   frq = frq[range(n/2)] #one sided frequency
   
   #stores Fourier data
   Ydata = np.fft.fft(npData)/npDLength ##fft computing and normalizing
   Ydata = Ydata[range(npDLength/2)]
   
   #get the Fourier frequencies
   freqs = np.fft.fftfreq(npDLength)
   
   #use for loop to construct angles list
   #based on the DEG_PER_ELEM and data length
   #must have angles length = data length
   angles = []
   
   for i in range(1, npDLenght):
   	angles.append(i * DEG_PER_ELEM)
   
   
   
   
   	
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
    