import re
import os
import glob
import sys
sys.path.append('/home/stephan/mc/ctqmc')
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/ctqmc")
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pylab
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/qising-SSE")
sys.path.append("/home/stephan/mc/qising-SSE")
from Fit import *
from texify import *
import scipy.integrate
import scipy.interpolate

#plt.rc('text', usetex=True)
#plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']

def closest_k_point(L):
	k = []
	for i in range(L):
		for j in range(L):
			k.append(float(i) * b1 / L + float(j) * b2 / L)
	q = k[0]
	d_q = np.linalg.norm(q - K)
	for i in range(len(k)):
		d_k = np.linalg.norm(k[i] - K)
		if d_k < d_q:
			q = k[i]
			d_q = d_k
	return [q, d_q]

filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/K_point/delta_kek.txt")[0]
with open(filename) as f:
	lines = (line for line in f if not line.startswith('L'))
	new_del = []
	for l in lines:
		new_del.append(','.join(l.split()))
	data = np.loadtxt(new_del, delimiter=',', skiprows=0)

cnt = 0
L_list = []
ax = plt.gca()

e_k = {"3" : "-1.73205", "4" : "-1", "5" : "-0.618034", "6" : "-1", "7" : "-0.554958", "8" : "-0.414214", "9" : "-0.68404", "10" : "-0.381966", "11" : "-0.309721", "13" : "0.29079", "14" : "0.24698"}
b1 = np.array([2.*np.pi/3., 2.*np.pi/np.sqrt(3.)])
b2 = np.array([2.*np.pi/3., -2.*np.pi/np.sqrt(3.)])
K = [2.*np.pi/3., 2.*np.pi/3./np.sqrt(3.)]

for i in range(len(data[:,0])):
	if len(L_list) == 0 or data[L_list[-1],0] != data[i,0]:
		L_list.append(i)
L_list.append(len(data[:,0]))

for i in range(len(L_list)-1):
	L = int(data[L_list[i], 0])
	
	q, d_q = closest_k_point(L)
	#data[L_list[i]:L_list[i+1],2] *= 3./2. * d_q / np.abs(float(e_k[str(L)]))
	#data[L_list[i]:L_list[i+1],3] *= 3./2. * d_q / np.abs(float(e_k[str(L)]))
	
	ax.plot( data[L_list[i]:L_list[i+1],1], data[L_list[i]:L_list[i+1],2], marker="o", color=color_cycle[cnt], markersize=10.0, linewidth=0.0, label="L="+str(int(L)))
	(_, caps, _) = ax.errorbar(data[L_list[i]:L_list[i+1],1], data[L_list[i]:L_list[i+1],2], yerr=data[L_list[i]:L_list[i+1],3], marker='None', capsize=8, linewidth=2.0, color=color_cycle[cnt])
	for cap in caps:
		cap.set_markeredgewidth(1.6)
		
	x = data[L_list[i]:L_list[i+1],1]
	y = data[L_list[i]:L_list[i+1],2]
	xnew = np.linspace(x.min(), x.max(), 300)
	f = scipy.interpolate.interp1d(x, y, kind=1)
	ax.plot( xnew, f(xnew), color=color_cycle[cnt], markersize=0.0, linewidth=2.0)
	cnt += 1

pylab.xlabel(r"$V$", fontsize=18)
pylab.ylabel(r"$\Delta_{kek}^{+} \sqrt{N}$", fontsize=18)
ax.axvline(1.355, color='k', linestyle='--')
plt.legend()
plt.tight_layout()
plt.show()
