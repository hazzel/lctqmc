import re
import os
import glob
import sys
sys.path.append('/home/stephan/mc/ctqmc')
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/ctqmc")
import numpy as np
from cdecimal import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import pylab
from ParseDataOutput import *
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/qising-SSE")
sys.path.append("/home/stephan/mc/qising-SSE")
from Fit import *
from texify import *
import scipy.integrate
import scipy.interpolate

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams.update({'font.size': 20})
plt.rc('legend',fontsize=12)

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen', 'cyan', 'magenta']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']
cnt = 0

list_of_files_R = []

#list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw-L6-s-*tprime*-theta8*.txt")
#list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw-L9-s-*tprime*-theta8*.txt")
#list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw-L12-s-*tprime*-theta8*.txt")
#list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw-L15-s-*tprime*-theta8*.txt")
#list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw-L18-s-*tprime*-theta8*.txt")
#list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw-L21-s-*tprime*-theta8*.txt")
list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw_11-*-s-*tprime*-theta320*.txt")
#list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw-L18-s-*tprime*-theta8*.txt")
#list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/R_cdw-L21-s-*tprime*-theta8*.txt")

#list_of_files_R.sort()
data_list_R = [ ( pylab.loadtxt(filename), label ) for label, filename in enumerate(list_of_files_R) ]

nu = 1./1.101#1./1.14
Vc = 0.4
c = -3.
omega = 0.779

ax = plt.gca()

for data, label in data_list_R:
	L = list_of_files_R[label].split("-")[1].replace("L", "")
	theta = list_of_files_R[label].split("-")[4].replace("theta", "").split(".")[0]
	x = data[:,0]
	#x = float(L)**(1./nu) * (data[:,0] - Vc) / Vc
	y = data[:,1]
	#ax = plt.subplot(121)
	ax.plot( x, y, marker="o", color=color_cycle[cnt], markersize=10.0, linewidth=0.0, label=r"$L="+L+", \Theta=" + theta + "$")
	if len(data[0,:]) > 2:
		(_, caps, _) = ax.errorbar(x, y, yerr=data[:,2], marker='None', capsize=8, color=color_cycle[cnt], linewidth=0.0, )
		for cap in caps:
			cap.set_markeredgewidth(1.6)
	xnew = np.linspace(x.min(), x.max(), 300)
	f = scipy.interpolate.interp1d(x, y, kind=1)
	ax.plot( xnew, f(xnew), color=color_cycle[cnt], markersize=0.0, linewidth=2.0)
	ax.legend()
	cnt += 1
	
	'''
	#ax = plt.subplot(122)
	x = float(L)**(1./nu) * (1. + c * float(L)**(-omega)) * (data[:,0] - Vc) / Vc
	ax.plot( x, y, marker="o", color=color_cycle[cnt], markersize=10.0, linewidth=0.0, label=r"$L="+L+"$")
	if len(data[0,:]) > 2:
		(_, caps, _) = ax.errorbar(x, y, yerr=data[:,2], marker='None', capsize=8, color=color_cycle[cnt], linewidth=0.0, )
		for cap in caps:
			cap.set_markeredgewidth(1.6)
	xnew = np.linspace(x.min(), x.max(), 300)
	f = scipy.interpolate.interp1d(x, y, kind=3)
	ax.plot( xnew, f(xnew), color=color_cycle[cnt], markersize=0.0, linewidth=2.0)
	ax.legend()
	cnt += 1
	'''
'''
list_of_files_R = []
list_of_files_R += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/mqmc/plot/thomas/R_chern-L*theta20*.txt")
data_list_R = [ ( pylab.loadtxt(filename), label ) for label, filename in enumerate(list_of_files_R) ]
for data, label in data_list_R:
	L = list_of_files_R[label].split("-")[1].replace("L", "")
	ax.plot( data[:,0], data[:,1], marker="o", color=color_cycle[cnt], markersize=10.0, linewidth=0.0, label=r"$L="+L+"$")
	if len(data[0,:]) > 2:
		(_, caps, _) = ax.errorbar(data[:,0], data[:,1], yerr=data[:,2], marker='None', capsize=8, color=color_cycle[cnt], linewidth=2.0, )
		for cap in caps:
			cap.set_markeredgewidth(1.6)
	x = data[:,0]
	y = data[:,1]
	xnew = np.linspace(x.min(), x.max(), 300)
	f = scipy.interpolate.interp1d(x, y, kind=1)
	ax.plot( xnew, f(xnew), color=color_cycle[cnt], markersize=0.0, linewidth=2.0)
	cnt += 1
'''
ax.set_xlabel(r"$V$", fontsize=18)
ax.set_ylabel(r"$R_{cdw}$", fontsize=18)
#ax.set_xlim([0., 0.6])
#ax.set_ylim([0., 3.1])
ax.legend()

plt.tight_layout()
plt.show()
