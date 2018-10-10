import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import scipy.optimize

plt.rc("mathtext", fontset="cm")
plt.rc("font", family="serif", size=16)
plt.rc('legend',fontsize=15)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('text', usetex=True)

prop_cycle = plt.rcParams['axes.prop_cycle']
color_cycle = prop_cycle.by_key()['color']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']
fillstyle = ['none', 'full']

# BZ:
#  __
# /  \
# \__/
a1 = np.array([1./2., -np.sqrt(3.)/2.])
a2 = np.array([1., 0.])
b1 = np.array([2.*np.pi, -2.*np.pi/3.*np.sqrt(3.)])
b2 = np.array([0., 4.*np.pi/np.sqrt(3.)])
K = { 6 : np.array([4.1887902, 0.]), 9 : np.array([4.1887902, 0.]), 12 : np.array([4.1887902, 0.]), 15 : np.array([4.1887902, 0.]) }

def FitFunction(x, a, b, c):
	return np.log(a*x)*b + c

def E(k):
	return np.linalg.norm(1. + np.exp(1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a1 - a2)))

def avg_group(vA0, vB0):
	vA, ind, counts = np.unique(vA0, return_index=True, return_counts=True) # get unique values in vA0
	vB = vB0[ind]
	for dup in vA[counts>1]: # store the average (one may change as wished) of original elements in vA0 reference by the unique elements in vB
		vB[np.where(vA==dup)] = np.average(vB0[np.where(vA0==dup)])
	return vA, vB

i_L = 0
i_gamma = 1
i_U = 2
i_UoverUc = 3
i_q = 4
i_kx = 5
i_ky = 6
i_delta = 7
i_sigma = 8

L_list = [ 6, 9, 12, 15 ]
gamma_list = [ 1. ]
U_list = [ 0.5 ]
v0 = 3.**0.5/2.

fig, ax = plt.subplots()
#plt.title(r"$\gamma = " + str(gamma_list[0]) + "$")
ax.set_xlabel(r"$a \Delta k$", fontsize=18)
ax.set_ylabel(r"$E / t$", fontsize=18)

filelist = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/analysis/data_gamma.txt")
filelist.sort()

x = np.linspace(0., 2., 1000)
ax.plot( x, v0 * x, color="k", ls="--", linewidth=2.0)
plt.axhline(y=1., ls="--", lw=2.0, color="grey")

cnt_U = 0
cnt_L = 0
for filename in filelist:
	with open(filename) as f:
		lines = (line for line in f if not line.startswith('L'))
		new_del = []
		for l in lines:
			new_del.append(','.join(l.split()))
		data = np.loadtxt(new_del, delimiter=',', skiprows=0)

	data_filter = list(filter(lambda x : x[i_L] in L_list and x[i_gamma] in gamma_list and x[i_U] in U_list, data))

	for U in U_list:
		data_filter_U = np.array(list(filter(lambda x : x[i_U] == U, data_filter)))
		for L in L_list:
			data_filter_L = np.array(list(filter(lambda x : x[i_L] == L, data_filter_U)))
			
			q_list = data_filter_L[:, i_q]
			kx_list = data_filter_L[:, i_kx]
			ky_list = data_filter_L[:, i_ky]
			delta_list = data_filter_L[:, i_delta]
			sigma_list = data_filter_L[:, i_sigma]
			
			#warping factor
			for i in range(1, len(delta_list)):
				q = q_list[i]
				k = np.array([kx_list[i], ky_list[i]])
				delta_list[i] *= v0 * q / E(k)
				sigma_list[i] *= v0 * q / E(k)
				
			x_grad, y_grad = avg_group(q_list, delta_list)
			v_q = np.diff(y_grad) / np.diff(x_grad) / v0
			avg_q = (x_grad[:-1] + x_grad[1:]) / 2.
			
			if U == 0.5:
				label = r"$U/t="+str(U)+",\ L="+str(L)+"$"
			else:
				label = r"$U/t="+str(U)+",\ L="+str(L)+"$"
			ax.plot( q_list, delta_list, marker=marker_cycle[cnt_L], color=color_cycle[cnt_L], markersize=12.0, fillstyle=fillstyle[cnt_U], linewidth=0.0, label=label)
			(_, caps, _) = ax.errorbar(q_list, delta_list, yerr=sigma_list, marker='None', capsize=8, color=color_cycle[cnt_L], linewidth=0.0)
			for cap in caps:
				cap.set_markeredgewidth(1.6)
			
			if L == 15:
				ax.plot( avg_q, v_q, marker=marker_cycle[cnt_L], color=color_cycle[cnt_L], markersize=14.0, fillstyle=fillstyle[cnt_U], linewidth=2.0, label=r"$v(\vec{q})$")
				parameter, perr = scipy.optimize.curve_fit( FitFunction, avg_q, v_q, p0=[1., 1., 1.], method='trf')
				px = np.linspace(avg_q[0], avg_q[-1], 1000)
				ax.plot(px, FitFunction(px, *parameter), 'k-', linewidth=2.0)
			
			cnt_L += 1
		cnt_U += 1

plt.legend(borderpad=0.05, labelspacing=0.075, facecolor='wheat', framealpha=0.5)
plt.tight_layout()
plt.savefig("gap_of_q.pdf", bbox_inches='tight', pad_inches = 0.1)

plt.show()
