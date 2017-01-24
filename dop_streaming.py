"""
dop_streaming.py
Requires common, constants, matplotlib, numpy, pyaudio, and streaming.
Pythonized code for Raytheon Low Cost Radar Initiative
Adapted from:
MIT Short Program 2013 Laptop Radar Course
(c) 2014 Massachusetts Institute of Technology
Python version 1, 6/1/2014

"""

import streaming as st
import matplotlib.animation as animation

def timer_callback(i,
                  objStrm,
                  objP,
                  objTI):
    objStrm.get_data()
    objStrm.format_raw()
    objP.process_doppler(objStrm.data)
    objP.update_freq(objTI.posFreqInd)
    objTI.update_doppler(objP.freq)

# Autocalibration
print("Calibrating...")
recordLen = 0.2
format    = 'float32'

S = st.Stream(recordLen, format)
P = st.dopProcess()
D = st.TIPlot('Doppler', 'Speed', 'Time')

S.open_stream()
S.get_data()
S.format_raw()
P.process_doppler(S.data)
P.update_freq(D.posFreqInd)

P.get_thresh()

print("Threshold: %f" % P.fThresh)
# End calibration

anim = animation.FuncAnimation(D.fig,
                               timer_callback,
                               fargs=(S, P, D),
                               interval=1)
#Why does blit-True cause errors?
