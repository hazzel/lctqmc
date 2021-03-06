import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

plt.rc("mathtext", fontset="cm")
plt.rc("font", family="serif", size=16)
plt.rc('legend',fontsize=18)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('text', usetex=True)

prop_cycle = plt.rcParams['axes.prop_cycle']
color_cycle = prop_cycle.by_key()['color']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']

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
U_list = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.4, 3.6, 3.75, 4, 4.5 ]
#U_list = [ 0.1, 0.2, 0.3, 0.4, 0.5, 0.65, 0.8, 1., 1.5, 2, 2.5, 3, 3.2, 3.4, 3.55, 3.6, 3.65, 3.75, 4., 4.5, 5. ]
q_list = [ 0 ]
v0 = 3.**0.5/2.

fig, ax = plt.subplots()
ax.set_xlabel(r"$U / t$", fontsize=18)
ax.set_ylabel(r"$[ E_n / \Delta k - v_0 n ] / v_0 t$", fontsize=18)
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#ax.text(0.05, 0.95, r"$\gamma=0$", transform=ax.transAxes, fontsize=18, verticalalignment='top', bbox=props)
ax.text(0.53, 0.92, r"$U_c(0)/t = 3.85$", transform=ax.transAxes, fontsize=16, verticalalignment='top', color="grey")

plt.axvline(x=3.85, ls="--", lw=2.0, color="grey")

filelist = glob.glob("data_of_q.txt")
filelist.sort()

for filename in filelist:
	with open(filename) as f:
		lines = (line for line in f if not line.startswith('L'))
		new_del = []
		for l in lines:
			new_del.append(','.join(l.split()))
		data = np.loadtxt(new_del, delimiter=',', skiprows=0)

	data_filter = list(filter(lambda x : x[i_L] in L_list and x[i_gamma] in gamma_list and x[i_U] in U_list, data))

	cnt_q = 0
	for L in L_list:
		data_filter_L = np.array(list(filter(lambda x : x[i_L] == L and x[i_U] in U_list, data_filter)))
		
		for q_i in q_list:
			data_filter_q = np.array(list(filter(lambda x : (np.abs(K[L]-b1/L*q_i - np.array([x[i_kx], x[i_ky]])) < 1E-6).all().any(), data_filter_L)))

			if len(data_filter_q > 0):
				delta_k = np.linalg.norm(b1) / L
				offset = v0 * q_i
				delta_list = (data_filter_q[:, i_delta] / delta_k - offset) / v0
				sigma_list = data_filter_q[:, i_sigma] / delta_k / v0
				
				ax.plot( U_list, delta_list, marker=marker_cycle[cnt_q], color=color_cycle[cnt_q], markersize=10.0, linewidth=0.0, label=f"$n={q_i}, L={L}$")
				(_, caps, _) = ax.errorbar(U_list, delta_list, yerr=sigma_list, marker='None', capsize=9, color=color_cycle[cnt_q], linewidth=2.0)
				for cap in caps:
					cap.set_markeredgewidth(1.8)
				cnt_q += 1

leg = plt.legend(borderpad=0.05, labelspacing=0.075)
leg.get_frame().set_linewidth(2.)
plt.tight_layout()
plt.savefig(f"pdf/scaled_gap_n{q_list[0]}.pdf", bbox_inches='tight', pad_inches = 0.1)

plt.show()
