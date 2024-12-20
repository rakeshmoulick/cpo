import numpy as np
import matplotlib.pyplot as plt # type: ignore
from numpy.fft import fft, ifft
import matplotlib as mp
# plt.style.use('seaborn-poster')

# sampling rate
sr = 2000
# sampling interval
ts = 1.0/sr # time step or time interval
t = np.arange(0,1,ts) # time array having 2000 or sr numbers of data points
# print(len(t))
# -------------------------------------------------------------------
freq = 1
x = 3*np.sin(2*np.pi*freq*t) # signal_1

freq = 4
x += np.sin(2*np.pi*freq*t) # signal_1 + signal_2

freq = 7   
x += 0.5* np.sin(2*np.pi*freq*t) # signal_1 + signal_2 + signal_3
# --------------------------------------------------------------------
X = fft(x)
N = len(X)
n = np.arange(N)
T = N/sr
freq = n/T 
# --------------------------------------------------------------------
figsize = np.array([300, 60]) #Figure size in mm (FOR SINGLE FIGURE)
figsize1 = np.array([100, 60])
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ----------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
ax[0].plot(t, x, 'r')
ax[0].set_xlabel('Time (s)')
ax[0].set_ylabel('Amplitude')
ax[0].grid(True)

ax[1].stem(freq, np.abs(X)/(sr/2), 'b', markerfmt=" ", basefmt="-b")
# To get the actual amplitudes of the spectrum, we have to normalize the output of (fft) by N/2, the number of samples
ax[1].set_xlabel('Freq (Hz)')
ax[1].set_ylabel('FFT Amplitude |X(freq)|')
ax[1].set_xlim(0, 10)
ax[1].grid(True)

ax[2].plot(t, ifft(X), 'r')
ax[2].set_xlabel('Time (s)')
ax[2].set_ylabel('Amplitude')
ax[2].grid(True)

plt.show()