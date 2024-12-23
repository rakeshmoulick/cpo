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
dpi = 600                         #Print resolution
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
# ------------------------------------------------------
k = 0 # file index
camera = Camera(fig)
while k<=NUM_TS:
	wpet = k*DT_coeff
	file_name1 = 'e%d.txt'%(k)
	file_name2 = 'i%d.txt'%(k)
	path = '../data/'
	
	# Load electron phase space data
	if os.path.exists(pjoin(path,file_name1)):
		xe, ye = np.loadtxt(pjoin(path,file_name1),unpack=True)
	else:
		print('No Data')
		exit()
	# Load ion phase space data
	if os.path.exists(pjoin(path,file_name2)):
		xi, yi = np.loadtxt(pjoin(path,file_name2),unpack=True)
	else:
		print('No Data')
		exit()
	
	#plt.cla() # Remove this command and you see a cluster of all previous timesteps		
	plt.scatter(xe, ye, s=0.3, alpha=0.5, color='r')
	plt.scatter(xi, yi, s=0.3, alpha=0.5, color='b')	
	
	# txt = "ts:"+str(k)
	txt = "$wpet: %.2f$"%(wpet)
	plt.text(2, -1.5E7,txt,ha='center',color='black',size=15)	
	#plt.xlim([-0.0005, 0.0025])
	#plt.ylim([-3E6, 3E6])
	plt.xlabel('$x$')
	plt.ylabel('$v$')
	#plt.title(k)
	plt.grid(True)	
	#plt.pause(0.1)	
	camera.snap()	
	k = k + write_interval_phase
animation = camera.animate()
#animation.save('cpo.gif', writer = 'pillow', fps=5.0)
animation.save('cpo.mp4', fps=5.0)	
plt.show()
