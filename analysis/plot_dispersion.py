import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
import glob
import sys
import io
import numpy as np
sys.path.append("/home/stephan/mc/qising-SSE")
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/qising-SSE")
from texify import *

def f1(k):
	return 1. + np.exp(1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a1 - a2))
def f3(k):
	'''
	return np.exp(-1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, 2.*a1-a2)) \
		+ np.exp(-1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a2-a1)) + np.exp(1.j * np.dot(k, a1+a2)) \
		+ np.exp(1.j * np.dot(k, 2.*a1)) + np.exp(1.j * np.dot(k, 2.*a1-2.*a2)) + np.exp(1.j * np.dot(k, a1-2.*a2))
	'''

	return np.exp(-1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, 2.*a1-a2))

def generate_path(t1, t3, P1, P2, N):
	x = []
	y = []
	#P1 -> P2
	for i in range(N):
		k = P1 + float(i) / N * (P2 - P1)
		x.append(float(i) / N * np.linalg.norm(P2 - P1))
		y.append(np.linalg.norm(t1*f1(k) + t3*f3(k)))
	return [np.array(x), np.array(y)]

'''
a1 = np.array([3./2., np.sqrt(3.)/2.])
a2 = np.array([3./2., -np.sqrt(3.)/2.])
b1 = np.array([2.*np.pi/3., 2.*np.pi/np.sqrt(3.)])
b2 = np.array([2.*np.pi/3., -2.*np.pi/np.sqrt(3.)])
K = np.array([2.*np.pi/3., 2.*np.pi/3./np.sqrt(3.)])
Kp = np.array([2.*np.pi/3., -2.*np.pi/3./np.sqrt(3.)])
M = np.array([2.*np.pi/3., 0.])
Gamma = np.array([0., 0.])
'''

#  __
# /  \
# \__/
a1 = np.array([1./2., -np.sqrt(3.)/2.])
a2 = np.array([1., 0.])
b1 = np.array([2.*np.pi, -2.*np.pi/3.**0.5])
b2 = np.array([0., 4.*np.pi/np.sqrt(3.)])
K = np.array([4.0/3. * np.pi, 0.])
Kp = np.array([2.0/3. * np.pi, 2.0/3.**0.5 * np.pi])
M = np.array([np.pi, np.pi / 3.**0.5])
Gamma = np.array([0., 0.])

t1 = 1.
t3 = 0.5
N = 10000
path_gamma_k = generate_path(t1, t3, Gamma, K, N)
path_k_m = generate_path(t1, t3, K, M, N)
path_m_gamma = generate_path(t1, t3, M, Gamma, N)

pos = [0., np.linalg.norm(Gamma-K), np.linalg.norm(Gamma-K)+np.linalg.norm(K-M), np.linalg.norm(Gamma-K)+np.linalg.norm(K-M)+np.linalg.norm(M-Gamma)]
#latexify()
fig, ax = plt.subplots()
#ax.set_aspect("equal")
ax.plot(path_gamma_k[0], path_gamma_k[1], c="b", lw=3.0)
ax.plot(path_gamma_k[0], -path_gamma_k[1], c="b", lw=3.0)
ax.plot(path_k_m[0]+pos[1], path_k_m[1], c="g", lw=3.0)
ax.plot(path_k_m[0]+pos[1], -path_k_m[1], c="g", lw=3.0)
ax.plot(path_m_gamma[0]+pos[2], path_m_gamma[1], c="orange", lw=3.0)
ax.plot(path_m_gamma[0]+pos[2], -path_m_gamma[1], c="orange", lw=3.0)

plt.xticks(pos, ['G', 'K', 'M', 'G'])
ax.set_ylabel(r"$\epsilon(k)$")

plt.tight_layout()
plt.show()
