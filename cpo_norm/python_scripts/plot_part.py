"""
plot_part.py
-------------------------------------------------------
OBJ: To plot the particle phase space 
-------------------------------------------------------
RUN AS: python3 plot_part.py ../input.ini
--------------------------------------------------------
DEVELOPER: DR. RAKESH MOULICK
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
import os.path
from os.path import join as pjoin
from celluloid import Camera
import sys
import configparser 
# ------------------------------------------------
config_path = sys.argv[1]
# Read the configuration file
config = configparser.ConfigParser()
config.read(config_path)
# ------------------------------------------------
NUM_TS = config.getint('time','NUM_TS')
write_interval_phase = config.getint('diagnostics','write_interval_phase')
DT_coeff = config.getfloat('time','DT_coeff')
# ------------------------------------------------
figsize = np.array([80,80/1.618]) #Figure size in mm (FOR SINGLE FIGURE)
dpi = 600                        #Print resolution
ppi = np.sqrt(1920**2+1200**2)/24.5 #Screen resolution
# ------------------------------------------------
mp.rc('text', usetex=False)
mp.rc('font', family='sans-serif', size=10, serif='Computer Modern Roman')
mp.rc('axes', titlesize=15)
mp.rc('axes', labelsize=15)
mp.rc('xtick', labelsize=15)
mp.rc('ytick', labelsize=15)
mp.rc('legend', fontsize=15)
# ------------------------------------------------
fig = plt.figure(figsize=(8,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])
# ------------------------------------------------
camera = Camera(fig)
# ------------------------------------------------
k = 0 # file index
while k<=NUM_TS:
	wpet = k*DT_coeff
	#  Load files
	file_name1 = 'e%d.txt'%(k)
	file_name2 = 'i%d.txt'%(k)
	data_path = '../data/'
	
	# Load electron phase space
	if os.path.exists(pjoin(data_path,file_name1)):
		xe, ye = np.loadtxt(pjoin(data_path,file_name1),unpack=True)
	else:
		print('No Data')
		exit()
	# Load ion phase space
	if os.path.exists(pjoin(data_path,file_name2)):
		xi, yi = np.loadtxt(pjoin(data_path,file_name2),unpack=True)
	else:
		print('No Data')
		exit()
	# ---------------------------------------------------------------
	#plt.cla() # Remove this command and you see a cluster of all previous timesteps		
	ax.scatter(xe, ye, s=0.3, alpha=0.5, color='r')
	ax.scatter(xi, yi, s=0.3, alpha=0.5, color='b')	
	
	# txt = "wpet: " + str(wpet)
	txt = "$wpet: %.2f$"%(wpet)
	ax.text(200,-15,txt,ha='center',color='black', size=15)	
	ax.set_xlabel('$x$')
	ax.set_ylabel('$v$')
	ax.grid(True)

	#plt.pause(0.1)	
	camera.snap()	
	k = k + write_interval_phase
animation = camera.animate()
#animation.save('CPO.gif', writer = 'pillow', fps=5.0)
path = '../data/files/'
animation.save(pjoin(path,'CPO.mp4'), fps=5.0)	
plt.show()
