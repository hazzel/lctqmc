import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import numpy as np

def f1(k):
	return 1. + np.exp(1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a1 - a2))

def f3(k):
	return np.exp(-1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, 2.*a1-a2))

def f_pi_flux(k):
	return -1. + np.exp(1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, a1 + a2))
	#phi = np.pi/4.
	#return 1.*np.exp(-1.j*phi) + np.exp(1.j * np.dot(k, a1))*np.exp(1.j*phi) + np.exp(1.j * np.dot(k, a2))*np.exp(1.j*phi) + np.exp(1.j * np.dot(k, a1 + a2))*np.exp(-1.j*phi)

# BZ:
#  __
# |  |
# |__|

a = 1.
a1 = np.array([1., 1.])*a
a2 = np.array([1., -1.])*a
delta = [np.array([-1., 0.])*a, np.array([1., 0.])*a, np.array([0., -1.])*a, np.array([0., 1.])*a]

b1 = np.array([np.pi/a, np.pi/a])
b2 = np.array([np.pi/a, -np.pi/a])
K = np.array([0.5*np.pi, 0.5*np.pi])
Kp = np.array([0.5*np.pi, -0.5*np.pi])
M = np.array([np.pi, 0.])
Gamma = np.array([0., 0.])

corner = []
corner.append([-np.pi, -np.pi])
corner.append([-np.pi, np.pi])
corner.append([np.pi, np.pi])
corner.append([np.pi, -np.pi])
corner.append([-np.pi, -np.pi])
corner = np.array(corner)

'''
# BZ:
#  /\
# |  |
#  \/

a1 = np.array([3./2., np.sqrt(3.)/2.])
a2 = np.array([3./2., -np.sqrt(3.)/2.])
delta = np.array([1./2., np.sqrt(3.)/2.])

b1 = np.array([2.*np.pi/3., 2.*np.pi/np.sqrt(3.)])
b2 = np.array([2.*np.pi/3., -2.*np.pi/np.sqrt(3.)])
K = np.array([2.*np.pi/3., 2.*np.pi/3./np.sqrt(3.)])
Kp = np.array([2.*np.pi/3., -2.*np.pi/3./np.sqrt(3.)])
M = np.array([2.*np.pi/3., 0.])
Gamma = np.array([0., 0.])

x0 = 2.0/(3.*3.**0.5)*np.pi
y0 = 2.0/3.0*np.pi
corner = []
corner.append([-y0, -x0])		#down left
corner.append([-y0, x0])		#down right
corner.append([0., 2.*x0])		#right
corner.append([y0, x0])			#up right
corner.append([y0, -x0])		#up left
corner.append([0., -2.*x0])	#left
corner.append([-y0, -x0])		#down left
corner = np.array(corner)
'''

t1 = 1.
t3 = 0.

dx, dy = 0.01, 0.01
xmin, xmax = corner[0][0], corner[2][0]
ymin, ymax = corner[0][1], corner[2][1]
y, x = np.mgrid[slice(ymin, ymax + dy, dy), slice(xmin, xmax + dx, dx)]
z = np.zeros(x.shape)
for i in range(z.shape[0]):
	for j in range(z.shape[1]):
		k = np.array([x[i, j], y[i, j]])
		#z[i, j] = np.linalg.norm(t1*f1(k) + t3*f3(k))
		z[i, j] = np.linalg.norm(f_pi_flux(k))

# x and y are bounds, so z should be the value *inside* those bounds.
# Therefore, remove the last value from the z array.
z = z[:-1, :-1]
levels = MaxNLocator(nbins=100).tick_values(z.min(), z.max())

# pick the desired colormap, sensible levels, and define a normalization
# instance which takes data values and translates those into levels.
cmap = plt.get_cmap('PiYG')
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

fig, ax = plt.subplots()
cf = ax.contourf(x[:-1, :-1] + dx/2., y[:-1, :-1] + dy/2., z, levels=levels, cmap=cmap, norm=norm)
fig.colorbar(cf, ax=ax)

ax.plot(corner[:,0], corner[:,1], marker="o", lw=4.0, ms=12.0, c="k")

fig.tight_layout()
plt.show()
