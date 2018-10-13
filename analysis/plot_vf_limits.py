import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import math
import scipy.optimize

def FitFunction(x, a, b, c):
	return a + b*x + c*x*x

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

gamma_list = [ 0. ]
U_list = [ 0.5, 1, 2, 2.5, 3, 3.2, 3.4, 3.6, 3.75 ]
#U_list = [ 2, 3, 3.75 ]
plot_tang, plot_tc = False, True

v0 = 3.**0.5/2.

fig, ax = plt.subplots()
ax.set_xlabel(r"$1 / L$", fontsize=18)
ax.set_ylabel(r"$v_F(U, L) / v_0$", fontsize=18)
ax.set_xlim([0., 0.28])
plt.axhline(y=1., lw=2., c="k")

props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(0.85, 0.96, r"$\gamma=0$", transform=ax.transAxes, fontsize=18, verticalalignment='top', bbox=props)

filelist = glob.glob("data_of_q.txt")
filelist.sort()

for filename in filelist:
	with open(filename) as f:
		lines = (line for line in f if not line.startswith('L'))
		new_del = []
		for l in lines:
			new_del.append(','.join(l.split()))
		data = np.loadtxt(new_del, delimiter=',', skiprows=0)

	data_filter = list(filter(lambda x : x[i_gamma] in gamma_list and x[i_U] in U_list, data))

	cnt_U = 0
	for U in U_list:
		data_filter_U = np.array(list(filter(lambda x : x[i_U] == U, data_filter)))
		data_split_L = np.array(np.split(data_filter_U, np.where(np.diff(data_filter_U[:, i_L]))[0] + 1))
		L_list = np.array([L_split[0][i_L] for L_split in data_split_L])
		
		if plot_tang:
			vf_list = np.array([ np.diff(L_split[:, i_delta]) / np.diff(L_split[:, i_q]) for L_split in data_split_L ]) / v0

			vf_min = np.array([ vf[0] for vf in vf_list ])
			vf_min_sigma = np.array([ (L_split[0, i_sigma]**2. + L_split[1, i_sigma]**2.)**0.5 / L_split[1, i_q] for L_split in data_split_L ]) / v0

			label = f"$U_{{\\text{{TANG}} }} = {U}$"
			ax.plot( 1./L_list, vf_min, marker=marker_cycle[cnt_U], color=color_cycle[cnt_U], markersize=12.0, fillstyle=fillstyle[0], linewidth=0.0, label=label)
			(_, caps, _) = ax.errorbar(1./L_list, vf_min, yerr=vf_min_sigma, marker='None', capsize=8, color=color_cycle[cnt_U], linewidth=0.0)
			for cap in caps:
				cap.set_markeredgewidth(2.0)

			parameter, perr = scipy.optimize.curve_fit( FitFunction, 1./L_list, vf_min, p0=[1., 1., 1.], method='trf')
			px = np.linspace(0, 1/L_list[0], 1000)
			ax.plot(px, FitFunction(px, *parameter), ls='--', color=color_cycle[cnt_U], linewidth=2.0)
			(_, caps, _) = ax.errorbar([0.], [parameter[0]], yerr=perr[0,0]**0.5, marker='None', capsize=8, color=color_cycle[cnt_U], linewidth=4.0)
			for cap in caps:
				cap.set_markeredgewidth(2.0)

		if plot_tc:
			vf_list = np.array([ L_split[1:, i_delta] / L_split[1:, i_q] for L_split in data_split_L ]) / v0

			vf_min = np.array([ vf[0] for vf in vf_list ])
			vf_min_sigma = np.array([ L_split[1, i_sigma] / L_split[1, i_q] for L_split in data_split_L ]) / v0

			label = f"$U_{{\\text{{TC}} }} = {U}$"
			ax.plot( 1./L_list, vf_min, marker=marker_cycle[cnt_U], color=color_cycle[cnt_U], markersize=12.0, fillstyle=fillstyle[1], linewidth=0.0, label=label)
			(_, caps, _) = ax.errorbar(1./L_list, vf_min, yerr=vf_min_sigma, marker='None', capsize=8, color=color_cycle[cnt_U], linewidth=0.0)
			for cap in caps:
				cap.set_markeredgewidth(2.0)

			parameter, perr = scipy.optimize.curve_fit( FitFunction, 1./L_list, vf_min, p0=[1., 1., 1.], method='trf')
			px = np.linspace(0, 1/L_list[0], 1000)
			ax.plot(px, FitFunction(px, *parameter), ls='-.', color=color_cycle[cnt_U], linewidth=2.0)
			(_, caps, _) = ax.errorbar([0.], [parameter[0]], yerr=perr[0,0]**0.5, marker='None', capsize=8, color=color_cycle[cnt_U], linewidth=4.0)
			for cap in caps:
				cap.set_markeredgewidth(2.0)

		cnt_U += 1

leg = plt.legend(borderpad=0.2, labelspacing=0.3, loc="lower right")
leg.get_frame().set_linewidth(2.)
plt.tight_layout()
if plot_tang and plot_tc:
	out_file = "pdf/vf_comp.pdf"
elif plot_tang:
	out_file = "pdf/vf_tang.pdf"
else:
	out_file = "pdf/vf_tc.pdf"
plt.savefig(out_file, bbox_inches='tight', pad_inches = 0.1)

plt.show()
