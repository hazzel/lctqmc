import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math

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

def E(k):
	return np.linalg.norm(1. + np.exp(1.j * np.dot(k, a1)) + np.exp(1.j * np.dot(k, a1 - a2)))

i_L = 0
i_gamma = 1
i_U = 2
i_q = 3
i_kx = 4
i_ky = 5
i_delta = 6
i_sigma = 7

L_list = [ 6, 9, 12, 15 ]
gamma_list = [ 0. ]
U_list = [ 0.5, 3.75 ]
v0 = 3.**0.5/2.

fig, ax = plt.subplots()
#plt.title(r"$\gamma = " + str(gamma_list[0]) + "$")
ax.set_xlabel(r"$a \Delta k$", fontsize=18)
ax.set_ylabel(r"$E / t$", fontsize=18)

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.48, 0.95, r"$\gamma=0$", transform=ax.transAxes, fontsize=18, verticalalignment='top', bbox=props)
ax.text(0.60, 0.35, r"$\text{40\% reduced } v_0$", transform=ax.transAxes, fontsize=18, verticalalignment='top')
ax.text(0.60, 0.74, r"$v_0$", transform=ax.transAxes, fontsize=18, verticalalignment='top')

filelist = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/analysis/data_of_q.txt")
filelist.sort()

x = np.linspace(0., 2., 1000)
ax.plot( x, v0*0.6 * x, color="k", ls="-.", linewidth=2.0)
ax.plot( x, v0 * x, color="k", ls="--", linewidth=2.0)

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
			
			
			if U == 0.5:
				label = r"$U/t="+str(U)+",\ \, \ L="+str(L)+"$"
			else:
				label = r"$U/t="+str(U)+",\ L="+str(L)+"$"
			ax.plot( q_list, delta_list, marker=marker_cycle[cnt_L], color=color_cycle[cnt_L], markersize=12.0, fillstyle=fillstyle[cnt_U], linewidth=0.0, label=label)
			(_, caps, _) = ax.errorbar(q_list, delta_list, yerr=sigma_list, marker='None', capsize=8, color=color_cycle[cnt_L], linewidth=0.0)
			for cap in caps:
				cap.set_markeredgewidth(1.6)
			if L == 15:
				ax.plot( q_list[:2], delta_list[:2], color=color_cycle[cnt_L], linewidth=2.0, ls="--")
			cnt_L += 1
		cnt_U += 1

plt.legend(borderpad=0.05, labelspacing=0.075, facecolor='wheat', framealpha=0.5)
plt.tight_layout()
plt.savefig("gap_of_q.pdf", bbox_inches='tight', pad_inches = 0.1)

plt.show()
