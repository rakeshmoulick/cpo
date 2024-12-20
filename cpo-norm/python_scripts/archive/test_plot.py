import numpy as np 
import matplotlib.pyplot as plt
from scipy.constants import value as constants
from scipy.constants import pi
# -------------------------------------------------------------------------
# https://stackoverflow.com/questions/40642061/how-to-set-axis-ticks-in-multiples-of-pi-python-matplotlib
# -------------------------------------------------------------------------

x = np.linspace(0, 2*pi, 100)
A = 0.05
d = -A
n = 1 + d*np.cos(x)

fig, ax = plt.subplots()
ax.plot(x,n)
ax.axhline(1, color='grey', linestyle='--',linewidth=2.0)

# ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi / 2))
# ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi / 12))
# ax.xaxis.set_major_formatter(plt.FuncFormatter(multiple_formatter()))
ax.set_xlabel('x')
ax.set_ylabel('n')
ax.grid(True)
plt.show()