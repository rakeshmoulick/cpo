import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt # type: ignore

plt.matplotlib.rc('text', usetex = False)
plt.matplotlib.rc('grid', linestyle = 'dotted')
plt.matplotlib.rc('figure', figsize = (6.4, 4.8)) # (width,height) inches
t = np.linspace(0, 1000, 10000)

A = 0.01
x = (A**2*t**3)/(24*np.sqrt(2))

y0 = sp.jv(0, x)
y1 = sp.jv(1, x)
y = (A/2)*np.sqrt(y0**2 + y1**2)*(2*np.pi)

plt.plot(t, y)

plt.xlabel('$t$')
plt.ylabel('${J}_n(x)$')
plt.grid(True)
plt.savefig('../bessel.jpg')

plt.show()
# save the data for later use by pgfplots
# np.savetxt('example-04.txt',list(zip(x,sp.jv(0,x),sp.jv(1,x),sp.jv(2,x),sp.jv(3,x),sp.jv(4,x),sp.jv(5,x))),fmt="% .10e")