
import numpy as np
import matplotlib.pyplot as plt
import os.path
from os.path import join as pjoin
from celluloid import Camera
k = 0
fig = plt.figure(figsize=(7,5))
left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
ax = fig.add_axes([left, bottom, width, height])
camera = Camera(fig)
while k<=50000:
	file_name1 = 'ec%d.txt'%(k)
	file_name2 = 'eh%d.txt'%(k)
	#file_name3 = 'i%d.txt'%(k)
	path = '../data/'
	
	# Load cold electron phase space
	if os.path.exists(pjoin(path,file_name1)):
		xc, yc = np.loadtxt(pjoin(path,file_name1),unpack=True)
	else:
		print('No Data')
		exit()
	# Load hot electrons phase space
	if os.path.exists(pjoin(path,file_name2)):
		xh, yh = np.loadtxt(pjoin(path,file_name2),unpack=True)
	else:
		print('No Data')
		exit()
	
	#plt.cla() # Remove this command and you see a cluster of all previous timesteps		
	plt.scatter(xh, yh, s=0.3, alpha=0.5, color='r')
	plt.scatter(xc, yc, s=0.3, alpha=0.5, color='b')	
	
	txt = "ts:"+str(k)
	plt.text(0.001,0.002,txt,ha='center',color='black')	
	#plt.xlim([-0.0005, 0.0025])
	#plt.ylim([-3E6, 3E6])
	plt.xlabel('x')
	plt.ylabel('v')
	#plt.title(k)
	plt.grid(True)	
	#plt.pause(0.1)	
	camera.snap()	
	k = k + 200
animation = camera.animate()
#animation.save('ebp.gif', writer = 'pillow', fps=5.0)
animation.save('ebp.mp4', fps=5.0)	
plt.show()
