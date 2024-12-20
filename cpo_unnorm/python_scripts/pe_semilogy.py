"""
pe_semilogy.py
-------------------------------------------------------
OBJ: To observe energy conservation by ploting 
the electric potential energy versus the kinetic energy 
-------------------------------------------------------
RUN AS: python3 pe_semilogy.py ../input.ini
--------------------------------------------------------
DEVELOPER: DR. RAKESH MOULICK
"""
import numpy as np
import matplotlib.pyplot as plt # type: ignore
import matplotlib as mp
from scipy.constants import value as constants
from scipy.constants import pi
import os.path
from os.path import join as pjoin
import scipy.integrate as intg
import sys
from scipy.signal import find_peaks
import configparser 
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Constants
eps0 = constants('electric constant')
kb = constants('Boltzmann constant')
me = constants('electron mass')
AMU = constants('atomic mass constant')
e = constants('elementary charge')
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
config_path = sys.argv[1]
# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)
# ---------------IMPORT REQUIRED DATA FROM INPUT.INP-----------
NUM_TS = config.getint('time','NUM_TS')
DT_coeff = config.getfloat('diagnostics','DT_coeff')
write_interval = config.getint('diagnostics','write_interval')
n0 = config.getfloat('population','density')
mi = config.getfloat('population','massI')
mi = mi*AMU
NC = config.getint('grid','NC')
nParticlesE = config.getint('population','nParticlesE')
nParticlesI = config.getint('population','nParticlesI')
# ------------------------------------------------------
ni0 = n0
ne0 = n0
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATA_TS = int(NUM_TS/write_interval) + 1

Lx = 2*pi
L = Lx/NC

wpi = np.sqrt(ni0*e**2/(eps0*mi)) # Ion Plasma Frequency
wpe = np.sqrt(ne0*e**2/(eps0*me)) # Electron Plasma Frequency

tcf = 1E9  # Time Conversion Factor to convert from seconds to nano-seconds, microseconds etc.

print("wpe:%g [rad/sec]"%(wpe))
print("Time Period:%g [ns]"%(2*pi*tcf/wpe))
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
file_name = 'processed_results_all.npz'
path = '../data/files/'    
#++++++++++++++++++++ UNPACK +++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    data = np.load(pjoin(path,file_name))
    x = data['x']
    nde = data['nde']
    ndi = data['ndi']
    phi = data['phi']
    EF = data['EF']
else:
    print('No data')
    exit()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
x = x.reshape(DATA_TS,(NC+1))
x = x[0,:] # SPATIAL GRID: Select only a column, all the columns have the same value
dx = x[1]-x[0] # Normalized StepSize
xl = NC*dx
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
"""
In each row of EF, there are NC+1 = 1025 columns. Data is put in the first row first
and then it goes into the second row. Therefore, all the values of the EF for the
first time step is written in first row in 1025 columns. The data of next time step
is in the second row and so on. Similar is the case for writing and reshaping x.
w_pe*t = NUM_TS*DT_coeff
"""
EF = EF.reshape(DATA_TS,(NC+1))
phi = phi.reshape(DATA_TS,(NC+1))

DT = DT_coeff/wpe
# create array of the w_pe*t
ts = np.linspace(0, NUM_TS*DT, DATA_TS)

#print("The shape of EF is: ", EF.shape)

# ++++++++++++++++++++++++++++++  Method-2 : Integrate Electric Field ++++++++++++++++++++++++++++++++
# each row of EF correspond to a separate time step. Hence, no of rows implies no of timesteps
PE1 = np.zeros(EF.shape[0])
for i in range(len(PE1)):
    xdata = x # un-normalized x-data   
    ydata = (EF[i,:]**2)  # un-normalized EF-square data
    # integrate the square electric field over the volume to get the un-normalized electric potential energy
    E_int = intg.trapz(ydata,xdata)    

    # multiply 0.5*eps0 to the un-normalized electric potential energy
    PE1[i] = 0.5*eps0*(E_int)

# Normalize by the thermal energy 
# Tot_cold_energy = ((electron_spwt)*nParticlesE)*Te

# Normalize the electric potential energy by its maximum value
# PE1 /= np.max(PE1)
# +++++++++++++++  Assign Appropriate Potential Energy +++++++++++++++++++++++
Pot_En = PE1 # use Method-2 to calculate Pot_En
#print(len(Pot_En), len(ts))
#print(np.shape(Field_En))


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 600                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24.5 #Screen resolution

mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=10)
mp.rc('axes', labelsize=10)
mp.rc('xtick', labelsize=10)
mp.rc('ytick', labelsize=10)
mp.rc('legend', fontsize=10)
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot-1: Only Potential Semilog Plot
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# fig,ax = plt.subplots(figsize=figsize/25.4,constrained_layout=True,dpi=ppi)
# #ax.semilogy(ts,Pot_En)
# ax.plot(ts,Pot_En)
# ax.set_xlabel('$\omega_{pe}t$')
# ax.set_ylabel('Eelectrostatic Field Energy')
# # ax.set_ylim(1E-4,1E2)
# # ax.legend(loc='upper right', framealpha=0.5)
# ax.grid(True)
# # plt.savefig('../pe_semilogy.png',dpi=dpi)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot-2: Compare potential and kinetic energies
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
file_name = 'ke.txt'
#++++++++++++++++++++ UNPACK ++++++++++++++++++++++++++++++++++++++++++++++++++
if os.path.exists(pjoin(path,file_name)):
    t, kei, kee = np.loadtxt(pjoin(path,file_name),unpack=True)
else:
    print('No data')
    exit()
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
Kin_En = (kei + kee)
# Kin_En /= np.max(Kin_En)    
Tot_En = Pot_En + Kin_En
# ----------------------------------------------------------------------------
# Compare the times 'ts' (constructed in this file) and 't' (from simulation)
# They should be the same
file = open("test.txt","w")
for i in range(len(t)):
    file.write(str(ts[i])+'\t'+str(t[i])+'\n')
# ----------------------------------------------------------------------------
fig,ax = plt.subplots(figsize=figsize/20,constrained_layout=True,dpi=ppi)

# ax.plot(ts, Tot_En, 'g-', label='$TE$')
ax.plot(ts*tcf, Pot_En*1E3, 'r-', label='$PE$')
ax.plot(t*tcf, Kin_En*1E3, 'b-', label='$KE$')

peaks, _ = find_peaks(Pot_En)
# print('Peaks are at indices:', peaks)
# ax.plot(ts[peaks]*tcf, Pot_En[peaks]*1E3, 'gd', label='$PE\;Peaks$')
print("Time period of one oscillation: %f [ns]"%((ts[peaks[2]]-ts[peaks[0]])*tcf))

# ax.set_xlim(0, 80)
ax.set_xlabel('$Time \; [ns]$')
ax.set_ylabel('$Energy$ [mJ]')

ax.legend(loc='upper right', framealpha=0.5)
plt.savefig(pjoin(path, 'energy_conservation.jpg'),dpi=dpi)
plt.show()

