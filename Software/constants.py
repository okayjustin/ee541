"""
List of constants used for 
cantenna_RTI.py
cantenna_dop.py
cantenna_sar.py
dop_streaming.py
RTI_streaming.py
constants.py
Version 1 6/1/2014
Version 2 6/29/2016
Version 3 7/15/2016
"""

AUDIODEVNAME   = 'iMic'
C              = 299.0e6 # Speed of light, m/s
CH_CHK_SAMPS   = 100000
CPI            = 0.50 # Coherent processing interval
FC             = 2590.0e6 # (Hz) Center frequency (connect VCO Vtune to +5)
FALSE          = 0
FSTART         = 2400e6 # (Hz) LFM start frequency
FSTOP          = 2480e6 # (Hz) LFM stop frequency
LAMBDA         = C/FC  # Wavelength
MAX_RANGE      = 300.0
MAX_SPEED      = 30.0 # Maximum speed to display m/s
NFFT           = 20000
NPULSE_CANCEL  = 2 # Number of pulses to use for canceller
NUM_BITS       = 16
NUM_PAD        = 64 # number of samples to pad for bandlimited interpolation & shifting
OLAP_FACT      = 8 # (unitless) number of overlapped pulse windows (1 for no overlap)
OSAMP_DOP      = 4 # (unitless) oversample factor for Doppler axis
OSAMP_TRG      = 16 # oversampling factor applied when interpolating the trigger signal
OSAMP_RNG      = 2 # oversampling factor applied when interpolating the range data
PI             = 3.14159265358979
SCHMITT_THRESH = 1000 #Schmitt trigger threshold
TP             = 10e-3 # minimum pulse length
TWOPI          = 2*PI
USE_MPH        = False
VERBOSE        = True

# Constants for RTI
BW = FSTOP - FSTART
DELTA_RG = C/(2*BW)

# Constants for SAR processing
TRP            = .1
DELTA_AZ       = -2*0.3048/12 # 2 inch antenna spacing
OSAMP_FACT     = 16.0
SCENE_SZ_AZ    = 140
SCENE_SZ_RG    = SCENE_SZ_AZ*3/4
ANGLE_LIMIT    = 45.0 # Angular limit in degrees
USE_TAPER      = False

# Constants for streaming
CHUNK          = 1024
DATA_FORMAT    = 'float32'
FCSTREAM       = 2.4e9
NUM_CHANNELS   = 2
NUM_ROWS       = 200
PULSE_BUFFER   = 500
RATE           = 44100
RECORD_LEN     = 0.2
STREAM_SCHMITT = 10000
#RECORD_SECONDS = 4.0
TIMER_PERIOD   = .2
#WAVE_OUTPUT_FILENAME = "output.wav"


