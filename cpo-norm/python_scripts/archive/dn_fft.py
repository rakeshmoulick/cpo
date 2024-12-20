
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
# ---------------------------------------------
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
eV = e
# ---------------------------------------------
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS = 30000 
write_interval = 5 
DT_coeff = 0.01
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
n0 = 1E10
mi = me
#-------------------------------------------------------
# SIM Vars
NC = 512
Time = 0
#-----------------------------------------------------
ni0 = n0 
ne0 = n0

file_name = 'processed_results_all.npz'
path = '../data/files/'

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATA_TS = int(NUM_TS/write_interval) + 1
L = 2*pi/NC
wp = np.sqrt(n0*e**2/(eps0*me)) # Total Plasma Frequency
DT = DT_coeff*(1.0/wp)
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    data = np.load(pjoin(path,file_name))
    x = data['x']
    ndi = data['ndi']
    nde = data['nde']
    rho = data['rho']
    phi = data['phi']   
    EF = data['EF']
else:
    print('No data')
    exit()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = x.reshape(DATA_TS,(NC+1))
xx = x[0,:] # SPATIAL GRID: Select only a column, all the columns have the same value
dx = xx[1]-xx[0] # Normalized StepSize

# create array of the w_pe*t
ts = np.arange(0, NUM_TS+write_interval, write_interval)*DT_coeff

# Reshape the data 
# -----------------------------------------------------------------------------
EF = EF.reshape(DATA_TS,(NC+1))
#rho = rho.reshape(DATA_TS,(NC+1))
phi = phi.reshape(DATA_TS,(NC+1))
nde = nde.reshape(DATA_TS,(NC+1))
ndi = ndi.reshape(DATA_TS,(NC+1))
# -----------------------------------------------------------------------------
dn = ndi - nde 
# print("The shape of delta_n (reduced delta_n) is: ", dn.shape)
# -----------------------------------------------------------------------------
#                               FFT
# -----------------------------------------------------------------------------
ts = 1000        # time step whose dn signal is required
wpet = ts*write_interval*DT_coeff 

y = dn[ts,:]     # array of dn for ts time step. 
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
# -----------------------------------------------------------------------------
figsize = np.array([300, 100]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 600                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ----------------------------------------------------------------------
fig, ax = plt.subplots(1, 2, figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
ax[0].plot(xx, y, 'r')
ax[0].set_xlabel('Position')
ax[0].set_ylabel('Amplitude')
ax[0].grid(True)



ax[1].stem(k, np.abs(F)/(sr/2), 'b', markerfmt=" ", basefmt="-b")
# To get the actual amplitudes of the spectrum, we have to normalize the output of (fft) by N/2, the number of samples
ax[1].set_xlabel('k')
ax[1].set_ylabel('FFT Amplitude')
ax[1].set_title('Modes at $\omega_{pe}t=%.2f$'%(wpet))
# ax[1].set_xlim(0, 10)
ax[1].grid(True)
plt.show()

