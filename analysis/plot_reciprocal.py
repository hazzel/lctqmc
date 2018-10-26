import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
import sys
import io
import numpy as np
sys.path.append("/home/stephan/mc/qising-SSE")
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/qising-SSE")
from texify import *

#latexify()
fig, ax = plt.subplots()
ax.set_aspect("equal")

Lx = 8
Ly = 8


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

print(np.linalg.norm(b1) * 3./2. * 2.**0.5)

'''
# BZ:
#  __
# /  \
# \__/

b1 = (2.*np.pi/np.sqrt(3.), -2.*np.pi/3.)
b2 = (0., 4.*np.pi/3.)
K = (2.*np.pi/3./np.sqrt(3.), 2.*np.pi/3.)
Kp = (4.*np.pi/(3.*np.sqrt(3.)), 0.)
M = (0., 2.*np.pi/3.)
Gamma = (0., 0.)
'''

x = []
y = []
for i in range(Lx + 1):
	for j in range(Ly + 1):
		x.append(float(i) / Lx * b1[0] + float(j) / Ly * b2[0])
		y.append(float(i) / Lx * b1[1] + float(j) / Ly * b2[1])

ax.scatter(x, y, s=125, c="k")
ax.plot([K[0]], [K[1]], c="r", marker="o", markersize=15, fillstyle="top", label="K")
ax.plot([Kp[0]], [Kp[1]], c="g", marker="o", markersize=15, fillstyle="top", label="Kp")
ax.plot([M[0]], [M[1]], c="orange", marker="o", markersize=15, fillstyle="top", label="M")
ax.plot([Gamma[0]], [Gamma[1]], c="b", marker="o", markersize=15, fillstyle="top", label="Gamma")
ax.arrow(0, 0, b1[0], b1[1], head_width=0.25, head_length=0.25, fc='k', ec='k')
ax.arrow(0, 0, b2[0], b2[1], head_width=0.25, head_length=0.25, fc='k', ec='k')
ax.arrow(b1[0], b1[1], b2[0], b2[1], head_width=0., head_length=0., fc='k', ec='k')
ax.arrow(b2[0], b2[1], b1[0], b1[1], head_width=0., head_length=0., fc='k', ec='k')

print(4.*2.**0.5*Lx * 3./2*np.linalg.norm(K - (Lx//3+1)*2*b1/Lx - (Lx//3+1)*b2/Lx))

#for kx, ky in sorted(zip(x, y), key=lambda n: n[0]):
#	print(kx, ky)

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
