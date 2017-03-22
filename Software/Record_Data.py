#Python code for recording 

import pyaudio
import wave
import serial
import numpy as np

#modify Fields according to what we need
CHUNK = 1024  #Number of Samples taken at individual times
FORMAT = pyaudio.paInt16  #data type extracted from pyaudio
CHANNELS = 1  #Numer of channels we are recording (1)
RATE = 44100  #Sampling Rate 
RECORD_TIME = 20  #to modify for different record lengths
WAVE_OUTPUT_FILENAME = "output.wav"

#creates a variable to do pyaudio
p = pyaudio.PyAudio()

stream = p.open(format=FORMAT,
                channels=CHANNELS,
                rate=RATE,
                input=True,
                frames_per_buffer=CHUNK)
                
ser = serial.Serial('/dev/tty.usbmodem1421', 9600);
print("* recording")

frames = []
 
#keep waiting until a 'tick'
while (ser != 'a'):
   pass
   
print("* recording")
for i in range(0, int(RATE / CHUNK * RECORD_TIME)):
   data = stream.read(CHUNK)
   frames.append(data)
   
print("finished recording")
   
#close the stream to stop recording
stream.stop_stream()
stream.close()
p.terminate()

#create the wave 
wf = wave.open(WAVE_OUTPUT_FILENAME, 'wb')
wf.setnchannels(CHANNELS)
wf.setsampwidth(p.get_sample_size(FORMAT))
wf.setframerate(RATE)
wf.writeframes(b''.join(frames))
wf.close()