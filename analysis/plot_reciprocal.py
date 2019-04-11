import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import glob
import sys
import io
import numpy as np
import scipy.optimize

def f1(k):
	return 1. + np.exp(1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a1 - a2))
def f3(k):
	return np.exp(-1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, 2.*a1-a2))
def f_pi_flux(k):
	return -1. + np.exp(1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a2)) + np.exp(1.j * np.dot(k, a1 + a2))
	#phi = np.pi/4.
	#return 1.*np.exp(-1.j*phi) + np.exp(1.j * np.dot(k, a1))*np.exp(1.j*phi) + np.exp(1.j * np.dot(k, a2))*np.exp(1.j*phi) + np.exp(1.j * np.dot(k, a1 + a2))*np.exp(-1.j*phi)
def FitFunction(x, a, b, c):
	return a*x*x + c

#latexify()
fig, ax = plt.subplots()
ax.set_aspect("equal")

Lx = 7
Ly = 7


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
#K = np.array([0.5*np.pi, 0.5*np.pi])
#Kp = np.array([0.5*np.pi, -0.5*np.pi])
K = np.array([np.pi, np.pi/2.])
Kp = np.array([np.pi, -np.pi/2.])
M = np.array([np.pi, 0.])
Gamma = np.array([0., 0.])

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
'''

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
d = []
for i in range(Lx + 1):
	for j in range(Ly + 1):
		x.append(float(i) / Lx * b1[0] + float(j) / Ly * b2[0])
		y.append(float(i) / Lx * b1[1] + float(j) / Ly * b2[1])
		q = np.array([x[-1], y[-1]])
		d.append(np.linalg.norm(q-K) / np.pi * Lx)
		#print(i*(Ly+1) + j)
		#print(q)
		#print("|K-q| = ", np.linalg.norm(q-K)*Lx)
		#print("E = ", np.linalg.norm(f_pi_flux(q)))
print("|K-q| = ", min(d))

ax.scatter(x, y, s=125, c="k")
ax.plot([K[0]], [K[1]], c="r", marker="o", markersize=15, fillstyle="top", label="X")
ax.plot([Kp[0]], [Kp[1]], c="g", marker="o", markersize=15, fillstyle="top", label="X'")
#ax.plot([M[0]], [M[1]], c="orange", marker="o", markersize=15, fillstyle="top", label="M")
ax.plot([Gamma[0]], [Gamma[1]], c="b", marker="o", markersize=15, fillstyle="top", label="$\Gamma$")
ax.arrow(0, 0, b1[0], b1[1], head_width=0.25, head_length=0.25, fc='k', ec='k')
ax.arrow(0, 0, b2[0], b2[1], head_width=0.25, head_length=0.25, fc='k', ec='k')
ax.arrow(b1[0], b1[1], b2[0], b2[1], head_width=0., head_length=0., fc='k', ec='k')
ax.arrow(b2[0], b2[1], b1[0], b1[1], head_width=0., head_length=0., fc='k', ec='k')

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.tight_layout(rect=[0, 0.03, 1, 0.95])


#######################

fig2, ax2 = plt.subplots()
ax2.set_aspect("equal")

t1 = 1.
t3 = 0.

q_abs = []
E_q = []
n_fit = Lx//8

for i in range(-Lx//2 + Lx%2, Lx//2 + 1):
	q = K + i*b1 / Lx
	q_abs.append(np.sign(i) * np.linalg.norm(K-q))
	#E_q.append(np.linalg.norm(t1*f1(q) + t3*f3(q)))
	E_q.append(np.linalg.norm(t1*f_pi_flux(q)))
ax2.plot(q_abs, E_q, c="k", marker="o", markersize=15, fillstyle="top", lw=3.)
#parameter, perr = scipy.optimize.curve_fit( FitFunction, q_abs[len(q_abs)//2-n_fit:len(q_abs)//2+n_fit], E_q[len(q_abs)//2-n_fit:len(q_abs)//2+n_fit], p0=[0.1, 0.1, 1.], method='trf')
#xnew = np.linspace(q_abs[len(q_abs)//2-n_fit], q_abs[len(q_abs)//2+n_fit], 10000)
#ax2.plot( xnew, FitFunction(xnew, *parameter), color="r", markersize=0.0, linewidth=3.0, ls="--")


plt.show()
