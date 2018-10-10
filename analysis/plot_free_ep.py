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
plt.rc('legend',fontsize=18)

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']

filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/V=0/delta_ep_K.txt")[0]

with open(filename) as f:
	lines = (line for line in f if not line.startswith('L'))
	new_del = []
	for l in lines:
		new_del.append(','.join(l.split()))
	data = np.loadtxt(new_del, delimiter=',', skiprows=0)

cnt = 0
L_list = []
ax = plt.gca()

e_k = {"4" : "-1", "5" : "-0.618034", "7" : "-0.554958", "8" : "-0.414214", "10" : "-0.381966", "11" : "-0.309721", "13" : "0.29079", "14" : "0.24698"}
b1 = np.array([2.*np.pi/3., 2.*np.pi/np.sqrt(3.)])
b2 = np.array([2.*np.pi/3., -2.*np.pi/np.sqrt(3.)])
K = [2.*np.pi/3., 2.*np.pi/3./np.sqrt(3.)]

ax.plot( 1./data[:,0], data[:,2], marker="o", color=color_cycle[cnt], markersize=10.0, linewidth=0.0)
(_, caps, _) = ax.errorbar(1./data[:,0], data[:,2], yerr=data[:,3], marker='None', capsize=8, color=color_cycle[cnt], linewidth=0.0)
for cap in caps:
	cap.set_markeredgewidth(1.6)

x = 1./data[:,0]
y = data[:,2]
xnew = np.linspace(x.min(), x.max(), 300)
f = scipy.interpolate.interp1d(x, y, kind=1)
ax.plot( xnew, f(xnew), color=color_cycle[cnt], markersize=0.0, linewidth=2.0)

pylab.xlabel(r"$1 / L$", fontsize=18)
#pylab.ylabel(r"$\Delta_{\text{cdw}} \ \sqrt{N} \ v_F \ |\boldsymbol{q}| \ / \ |E(\boldsymbol{K}+\boldsymbol{q})|$", fontsize=18)
pylab.ylabel(r"$\Delta_{\text{ep}} \ \sqrt{N}$", fontsize=18)
pylab.xlim(left=0)

plt.tight_layout()
plt.savefig("free_ep.pdf", bbox_inches='tight', pad_inches = 0.1)

plt.show()
