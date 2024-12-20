
import numpy as np
from numpy import diff
from numpy.fft import fft, ifft
import matplotlib.pyplot as plt # type: ignore
import matplotlib as mp
from mpl_toolkits import mplot3d
from scipy.constants import value as constants
from scipy.constants import pi
import scipy.special as sp
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
m = 0 # Time data index
# ----------------------------------------------------------------------------
Time = np.array([])
Amplitude = np.array([])
# ---------------------------------------------------------------------------- 
MaxStp = DATA_TS
while m<MaxStp:
    wpet = m*write_interval*DT_coeff    
    # -----------------------------------------------------------------------------
    y = dn[m,:]     # array of dn for ts time step. 
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
    #amp = F[1]
    Time = np.append(Time, wpet)
    Amplitude = np.append(Amplitude, amp)
    m += 1 
# ----------------------------------------------------------------------------
#                           Analytical Bessel Calculations
# -----------------------------------------------------------------------------
t = np.linspace(0, np.max(Time), 10000)
A = 0.05
x = (A**2*t**3)/(24*np.sqrt(2))
y0 = sp.jv(0, x)
y1 = sp.jv(1, x)
y = (A/2)*np.sqrt(y0**2 + y1**2)*2*np.pi # In the main expression amplitude is A/2, but I don't think it is there
# -----------------------------------------------------------------------------
figsize = np.array([150, 100]) #Figure size in mm (FOR SINGLE FIGURE)
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
fig, ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
ax.plot(Time, (np.pi)*np.abs(np.real(Amplitude)), color='black', linestyle=':')
ax.plot(t, y, color='red', linestyle='-')
ax.set_xlim(0, 300)
ax.set_ylim(0, 0.16)
ax.set_xlabel('$\omega_{pe}t$')
ax.set_ylabel('$\delta n_{d}$')
plt.savefig(pjoin('../','fluctuation.png'),dpi=dpi)
plt.show()


