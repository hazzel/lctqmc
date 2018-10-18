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
	return a*x*x + b*x + c

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams.update({'font.size': 20})
plt.rc('legend',fontsize=18)

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']

list_of_files = []
list_of_files += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/crossing.txt")
list_of_files.sort()
cnt = 0

for filename in list_of_files:
	with open(filename) as f:
		lines = (line for line in f if not line.startswith('L'))
		new_del = []
		for l in lines:
			new_del.append(','.join(l.split()))
		data = np.loadtxt(new_del, delimiter=',', skiprows=0)

	fig, ax1 = plt.subplots()
	plt.title(r"$\text{" + filename.split("/")[-1].replace("_", "-") + "}$")

	ax1.plot( 1./data[:,0], data[:,1], marker="o", color="b", markersize=10.0, linewidth=2.0, label=r"E")
	(_, caps, _) = ax1.errorbar(1./data[:,0], data[:,1], yerr=data[:,2], marker='None', capsize=8, color="b", linewidth=0.0, )
	for cap in caps:
		cap.set_markeredgewidth(1.6)
	
	n = 6
	parameter, perr = scipy.optimize.curve_fit( FitFunction, 1./data[-n:,0], data[-n:,1], p0=[0.1, 0.1, 1.], method='trf')
	xnew = np.linspace(0, 1./data[-n,0], 300)
	ax1.plot( xnew, FitFunction(xnew, *parameter), color="k", markersize=0.0, linewidth=2.0)
	print(f"R* = {parameter[2]} +- {perr[2,2]}")

	ax1.set_ylim(bottom = 0)
	ax1.set_xlabel(r"$1 / L$", fontsize=18)
	ax1.set_ylabel(r"$V_c$", fontsize=18)
	ax1.legend()
	cnt += 1
plt.tight_layout()
plt.show()
