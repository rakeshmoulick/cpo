import numpy as np
from numpy import diff
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt # type: ignore
import matplotlib as mp
from mpl_toolkits import mplot3d
from scipy.constants import value as constants
from scipy.constants import pi
import os.path
from os.path import join as pjoin
from celluloid import Camera # type: ignore
import sys

t = np.linspace(0, 1000, 1000)
A = 0.01
x = np.linspace(0, 2*pi, 100)

nd = np.zeros([len(t), len(x)])
# print(nd)

for i in range(len(t)):
    for j in range(len(x)):
        nd[i,j] = A*np.cos(x[j])*np.cos(t[i]) + (A**2)*np.cos(2*x[j])*(0.25)*t[i]*np.sin(t[i])
        # nd[i,j] = A*np.cos(x[j])*np.cos(t[i])
        # nd[i,j] = (A**2)*np.cos(2*x[j])*(0.25)*t[i]*np.sin(t[i])

# print(nd[1, :])
# -------------------------------------------------------------------------------- 
Time = np.array([])
Amplitude = np.array([])
m = 0
MaxStp = len(t)
# --------------------------------------------------------------------------------
figsize = np.array([150, 100]) #Figure size in mm (FOR SINGLE FIGURE)
figsize1 = np.array([300, 100])
dpi = 600                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# --------------------------------------------------------------------------------
fig, ax = plt.subplots(1,2,figsize=figsize1/25.4,constrained_layout=True,dpi=ppi)
while m<MaxStp:    
    # -----------------------------------------------------------------------------
    wpet = t[m]
    y = nd[m,:]     # array of dn for ts time step. 
    sr = len(y)      # number of samples in the signal
    F = fft(y)       # FFT of the signal
    N = len(F)       # length of the FFT array F
    n = np.arange(N) # Create an array of N length
    T = N/sr         # Divide the length of the FFT array of F by the no. of samples of original signal
    k = n/T          # Find the inverse space/frequency 

    # Since the amplitudes of FFT will be symmetric about the whole range of k, 
    # take just half of the array of k and F for plotting
    halflen = np.array(F.shape, dtype=int)//2
    # print(halflen, type(halflen))
    k = k[:halflen[0]]
    F = F[:halflen[0]]    
    
    amp = F[1]/(sr/2) # normalized amplitude of the first modal value (k=1) not for k=0. This is the dominant mode at t=0    
    Time = np.append(Time, wpet)
    Amplitude = np.append(Amplitude, amp)
    # ------------------------------------------------------------------------
    ax[0].cla()
    ax[0].plot(x,nd[m,:],color='black')
    ax[0].set_ylim(-0.025, 0.025)
    ax[0].set_xlabel('x')
    ax[0].set_ylabel('$\delta n_{d}$')
    ax[0].grid(True)

    ax[1].cla()
    ax[1].stem(k, np.abs(F)/(sr/2), 'b', markerfmt=" ", basefmt="-b")
    ax[1].set_ylim(0, 0.05)
    ax[1].set_xlabel('$k$')
    ax[1].set_ylabel('$Amplitude$')
    ax[1].set_title("$\omega_{pe}t$ : %.2f"%(wpet))
    ax[1].grid(True)

    # ------------------------------------------------------------------------
    plt.pause(0.5)
    m += 10 
# ----------------------------------------------------------------------   
fig, ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
ax.plot(Time, np.abs(np.real(Amplitude)), color='black', linestyle='--')
ax.set_xlabel('$\omega_{pe}t$')
ax.set_ylabel('$\delta n_{d}$')
plt.show()