"""
dispersion_npz_EF.py
-------------------------------------------------------
OBJ: To plot the dispersion relation 
-------------------------------------------------------
RUN AS: python3 dispersion_npz_EF.py ../input.ini
--------------------------------------------------------
DEVELOPER: DR. RAKESH MOULICK
"""
import numpy as np
import matplotlib.pyplot as plt # type: ignore
import matplotlib as mp
from scipy.constants import value as constants
import os.path
from os.path import join as pjoin
import sys
import configparser 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
config_path = sys.argv[1]
# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)
#------------------------------------------------------------------------------
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
# ---------------IMPORT REQUIRED DATA FROM INPUT.INP-----------
NUM_TS = config.getint('time','NUM_TS')
DT_coeff = config.getfloat('time','DT_coeff')
write_interval = config.getint('diagnostics','write_interval')
n0 = config.getfloat('population','density')
mi = config.getfloat('population','massI')
mi = mi*AMU
NC = config.getint('grid','NC')
nParticlesE = config.getint('population','nParticlesE')
nParticlesI = config.getint('population','nParticlesI')
# ------------------------------------------------------

file_name = 'processed_results_all.npz'
path = '../data/files/'

#-----------------------------------------------------
ni0 = n0 
ne0 = n0

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATA_TS = int(NUM_TS/write_interval) + 1
wpe = np.sqrt(ne0*e**2/(eps0*me)) # Electron Plasma Frequency
Lx = 2*np.pi
L  = Lx/NC
W = wpe
DT = DT_coeff*(1.0/W)
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
# ----------------------------------------------------------------------------
xl = NC*dx
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
EF = EF.reshape(DATA_TS,(NC+1))

print("The shape of EF is: ", EF.shape)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
wpet_1 = 0
wpet_2 = NUM_TS*DT_coeff
y1 = wpet_1/(DT_coeff*write_interval)
y2 = wpet_2/(DT_coeff*write_interval)
E = EF[int(y1):int(y2),:]

print("The shape of E (reduced EF) is: ", E.shape)
#+++++++++++++++++++++++++ FFT of E-field data ++++++++++++++++++++++++++++++++
F = np.fft.fftn(E, norm='ortho')
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define Omega and K
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
NUM_TS1 = wpet_1/DT_coeff
NUM_TS2 = wpet_2/DT_coeff
actual_sim_time = (NUM_TS2 - NUM_TS1)*(DT) # Previously it was (NUM_TS*DT)
# -----------------------------------------------------------------------------
omega = 2*np.pi*np.arange(NUM_TS)/(actual_sim_time) #(DT*NUM_TS) #(actual_sim_time) # Unnormalized Omega
k     = 2*np.pi*np.arange(NC+1)/(NC*dx*L)       # un-normalized k
print('Length of k: ',len(k))
print('Max of k: ',np.max(k))

Omega, K = np.meshgrid(omega, k, indexing='ij')
print('Shape of Omega: ',Omega.shape)
#------------------------------------------------------------------------------
halflen = np.array(F.shape, dtype=int)//2
Omega = Omega[:halflen[0],:halflen[1]]
K = K[:halflen[0],:halflen[1]]
F = F[:halflen[0],:halflen[1]]
Omega /=W # Normalized Omega
K = K*L  # Normalized K
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Z = np.log(np.abs(F))
#print(F)
#+++++++++++++++++++++ Raw Analytic Part Calculation ++++++++++++++++++++++++++++++

raw_analytic = True
if raw_analytic:    
    # wla = np.sqrt((wpe**2*(1 + 3*k**2*L**2)));
    wla = np.sqrt(wpe**2)*np.ones(len(k));    
    ep = wla/W    # Electron Plasma Wave (Langmuir Wave)  
   
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 1200                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
extent = [np.min(K), np.max(K), np.min(Omega), np.max(Omega)]
ax.imshow(Z, extent=extent, interpolation='nearest', aspect='auto', origin='lower')

# plt.pcolor(K, Omega, Z,cmap='rainbow',shading='auto',vmin=np.min(Z),vmax=np.max(Z))
# plt.pcolor(K, Omega, Z,cmap='rainbow',shading='auto',vmin=-5,vmax=8)
# print(np.min(Z), np.max(Z))

# Note: Changing vmin & vmax will change the depth of the plots.
# cbar = plt.colorbar()
# cbar.set_label('$\zeta$')

if raw_analytic:
    plt.plot(k*L, ep, color='w', linestyle='--', lw = 1.0, alpha = 0.4, label='$EPW$')

ax.set_xlabel('$\kappa \; [kL]$')
ax.set_ylabel('$\Omega \; [\omega/\omega_{pe}]$')

ax.set_xlim([0, 3.0])
ax.set_ylim([0, 3.0])
leg = ax.legend(loc='upper right',framealpha=0.5)
plt.savefig(pjoin(pjoin(path,'dispersion.png')),dpi=dpi)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
plt.show()
