#Python code to process the waveFile generated in Record_Data.py
import common as com
import constants as C
import scipy.signal as sig
import os.path
from matplotlib import pyplot as plt
from numpy import (append, arange, argwhere, array, delete, diff, exp, floor,
 log10, median, nonzero, ones, sign, sin, tile, zeros)
from scipy import fftpack
from time import time

#modify Fields according to what we need
CHUNK = 1024  #Number of Samples taken at individual times
FORMAT = pyaudio.paInt16  #data type extracted from pyaudio
CHANNELS = 1  #Numer of channels we are recording (1)
RATE = 44100  #Sampling Rate 
RECORD_TIME = 20  #to modify for different record lengths
RPM = 5  #to modify according to actual speed of rotation

def polar_Image(waveFile):

	# Display if file is not found
   if not os.path.isfile(wavFile):
   	print("File %s not found." % (wavFile))
   
   #take i
   
    