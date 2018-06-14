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

def FitFunction(x, a, b, c):
	return b*np.exp(-c*x)+a

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams.update({'font.size': 20})
plt.rc('legend',fontsize=18)

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']

list_of_files = []
list_of_files += glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/benchmark/L3-s-tprime=0.5-theta.txt")
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
	print filename.split("/")[-1].replace("_", "-")
	plt.title(r"$\text{" + filename.split("/")[-1].replace("_", "-") + "}$")

	ax1.plot( data[:,0], data[:,1], marker="o", color="b", markersize=10.0, linewidth=2.0, label=r"E")
	(_, caps, _) = ax1.errorbar(data[:,0], data[:,1], yerr=data[:,2], marker='None', capsize=8, color="b", linewidth=0.0, )
	for cap in caps:
		cap.set_markeredgewidth(1.6)
		
	#parameter, perr = scipy.optimize.curve_fit( FitFunction, data[:,0], data[:,1], p0=[0.1, 0.1, 1.], method='trf')
	
	ax2 = ax1.twinx()
	ax2.plot( data[:,0], data[:,3], marker="o", color="g", markersize=10.0, linewidth=2.0, label=r"$R_{\text{cdw}}$")
	(_, caps, _) = ax2.errorbar(data[:,0], data[:,3], yerr=data[:,4], marker='None', capsize=8, color="g", linewidth=0.0, )
	for cap in caps:
		cap.set_markeredgewidth(1.6)

	#ax.set_xlabel(r"$\Delta \tau$", fontsize=18)
	ax1.set_xlabel(r"$\Theta$", fontsize=18)
	ax1.legend()
	ax2.legend()
	cnt += 1
plt.tight_layout()
plt.show()
