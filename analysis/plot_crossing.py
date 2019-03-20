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

def FitFunction(x, a, b, c):
	return b*x + c

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams.update({'font.size': 20})
plt.rc('legend',fontsize=18)

prop_cycle = plt.rcParams['axes.prop_cycle']
ecolor = prop_cycle.by_key()['color']
marker = ['o','s','D','^','o','s','D','^']
fillstyle = ['none', 'full']

list_of_files = []
list_of_files += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/crossing.txt")
list_of_files.sort()
cnt = 0

tprime_list = [ 0.4, 0.45, 0.5, 0.55, 0.60 ]

for filename in list_of_files:
	with open(filename) as f:
		lines = (line for line in f if not line.startswith('L'))
		new_del = []
		for l in lines:
			new_del.append(','.join(l.split()))
		data = np.loadtxt(new_del, delimiter=',', skiprows=0)

	fig, ax1 = plt.subplots()
	plt.title(r"$\text{" + filename.split("/")[-1].replace("_", "-") + "}$")
	ax1.set_ylim(bottom = 0)
	ax1.set_xlabel(r"$1 / L$", fontsize=18)
	ax1.set_ylabel(r"$V_c$", fontsize=18)
	
	for tprime in tprime_list:
		data_filter = np.array(list(filter(lambda x : x[1] == tprime, data)))

		ax1.plot( 1./data_filter[:,0], data_filter[:,2], marker="o", color=ecolor[cnt], markersize=10.0, linewidth=2.0, label=f"$t^{{\prime \prime}}={tprime}$")
		(_, caps, _) = ax1.errorbar(1./data_filter[:,0], data_filter[:,2], yerr=data_filter[:,3], marker='None', capsize=8, color=ecolor[cnt], linewidth=0.0, )
		for cap in caps:
			cap.set_markeredgewidth(1.6)
	
		n = min(len(data_filter), 5)
		parameter, perr = scipy.optimize.curve_fit( FitFunction, 1./data_filter[-n:,0], data_filter[-n:,2], p0=[0.1, 0.1, 1.], method='trf')
		xnew = np.linspace(0, 1./data_filter[-n,0], 300)
		ax1.plot( xnew, FitFunction(xnew, *parameter), color=ecolor[cnt], markersize=0.0, linewidth=2.0, ls="--")
		print(f"R* = {parameter[2]} +- {perr[2,2]}")
	
		cnt += 1
	ax1.legend()
plt.tight_layout()
plt.show()
