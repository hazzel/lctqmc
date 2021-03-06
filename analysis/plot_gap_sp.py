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

def LinearFunction(x, a, b):
	return a + b*x

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams.update({'font.size': 20})
plt.rc('legend',fontsize=20)

ecolor = ['#cc2909', '#efc600', '#60af92', '#3e77af','#3e77af','#60af92','#efc600','#cc2909']
fcolor = ['#ea6868', '#eddea2', '#99d1b9', '#a3c1e5','#a3c1e5','#99d1b9','#eddea2','#ea6868']
marker = ['o','s','D','^','o','s','D','^']

def closest_k_point(L):
	k = []
	for i in range(L):
		for j in range(L):
			k.append(float(i) * b1 / L + float(j) * b2 / L)
	q = k[0]
	d_q = np.linalg.norm(q - K)
	for i in range(len(k)):
		d_k = np.linalg.norm(k[i] - K)
		if d_k < d_q and d_k > 1E-12:
			q = k[i]
			d_q = d_k
	return [q, d_q]

def E(k):
	return np.linalg.norm(1. + np.exp(1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a1 - a2)))

#filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/sp-s-*tprime*-theta*.txt")[0]
#filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/gapped_spectrum/delta_sp.txt")[0]
#filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/K_point/delta_sp.txt")[0]
filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/K_point/delta_sp_kappa^2=0.txt")[0]
#filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/science_paper/delta_sp.txt")[0]
#filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/science_paper/delta_sp_q.txt")[0]
with open(filename) as f:
	lines = (line for line in f if not line.startswith('L'))
	new_del = []
	for l in lines:
		new_del.append(','.join(l.split()))
	data = np.loadtxt(new_del, delimiter=',', skiprows=2)

cnt = 0
L_list = []
ax = plt.gca()

e_k = {"3" : "-1.73205", "4" : "-1", "5" : "-0.618034", "6" : "-1", "7" : "-0.554958", "8" : "-0.414214", "9" : "-0.68404", "10" : "-0.381966", "11" : "-0.309721", "12" : "0.517638", "13" : "0.29079", "14" : "0.24698"}
a1 = np.array([3./2., np.sqrt(3.)/2.])
a2 = np.array([3./2., -np.sqrt(3.)/2.])
delta = np.array([1./2., np.sqrt(3.)/2.])
b1 = np.array([2.*np.pi/3., 2.*np.pi/np.sqrt(3.)])
b2 = np.array([2.*np.pi/3., -2.*np.pi/np.sqrt(3.)])
K = np.array([2.*np.pi/3., 2.*np.pi/3./np.sqrt(3.)])

for i in range(len(data[:,0])):
	if len(L_list) == 0 or data[L_list[-1],0] != data[i,0]:
		L_list.append(i)
L_list.append(len(data[:,0]))

for i in range(len(L_list)-1):
	L = int(data[L_list[i], 0])
	
	#q, d_q = closest_k_point(L)
	#data[L_list[i]:L_list[i+1],2] *= 3./2. * d_q / np.abs(float(e_k[str(L)]))# / d_q / (2.*L*L)**0.5
	#data[L_list[i]:L_list[i+1],3] *= 3./2. * d_q / np.abs(float(e_k[str(L)]))# / d_q / (2.*L*L)**0.5
	
	#q = b1 / L + 0*b2 / L
	#d_q = np.linalg.norm(q)
	#data[L_list[i]:L_list[i+1],2] *= 3./2. * d_q / E(K + q)
	#data[L_list[i]:L_list[i+1],3] *= 3./2. * d_q / E(K + q)
	
	#data[L_list[i]:L_list[i+1],2] /= (2.*L*L)**0.5
	#data[L_list[i]:L_list[i+1],3] /= (2.*L*L)**0.5
	
	ax.errorbar(data[L_list[i]:L_list[i+1],1], data[L_list[i]:L_list[i+1],2], yerr=data[L_list[i]:L_list[i+1],3], markeredgewidth=1.6, capsize=8, fmt=marker[cnt]+'-', ms=10.0, color=ecolor[cnt], ecolor=ecolor[cnt], mfc=fcolor[cnt], label=f"$L={L}$")
	
	#x = data[L_list[i]:L_list[i+1],1]
	#y = data[L_list[i]:L_list[i+1],2]
	#xnew = np.linspace(x.min(), x.max(), 300)
	#f = scipy.interpolate.interp1d(x, y, kind=1)
	#ax.plot( xnew, f(xnew), color=color_cycle[cnt], markersize=0.0, linewidth=2.0)
	
	'''
	nmin = 0
	nmax = 7
	parameter, perr = fit_function( [6., 1.], x[nmin:nmax], y[nmin:nmax], LinearFunction, datayerrors=data[L_list[i]:L_list[i+1],3][nmin:nmax])
	px = np.linspace(x[0], x[-1], 1000)
	ax.plot(px, LinearFunction(px, *parameter), 'k-', linewidth=3.0)
	print "L = ", L, ", v = ", parameter[1], " +- ", perr[1], ", v0 = ", parameter[0], " +- ", perr[0]
	'''
	
	cnt += 1

pylab.xlabel(r"$V/t$", fontsize=18)
#pylab.ylabel(r"$\Delta_{\text{cdw}} \ \sqrt{N} \ v_F \ |\boldsymbol{q}| \ / \ |E(\boldsymbol{K}+\boldsymbol{q})|$", fontsize=18)
#pylab.ylabel(r"$\Delta_{\text{sp}} \ \sqrt{N}$", fontsize=18)
pylab.ylabel(r"$\Delta (\Psi_T)\ \times \sqrt{N_s}$", fontsize=18)
#pylab.ylabel(r"$\Delta (\vec{K} - \vec{b}_1 / L)\ \sqrt{N}$", fontsize=18)
pylab.ylim(ymax = 6)
ax.axvline(1.355, color='k', linestyle='--', lw=2.0)
ax.axhline(3.9, color='k', linestyle='--', lw=2.0)

plt.legend(borderpad=0.15, labelspacing=0.05, framealpha=1.0)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("pdf/gap_sp.pdf", bbox_inches='tight', pad_inches = 0.05)

plt.show()
