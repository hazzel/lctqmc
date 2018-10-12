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

filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/suscept/R_cdw.txt")[0]

with open(filename) as f:
	lines = (line for line in f if not line.startswith('L'))
	new_del = []
	for l in lines:
		new_del.append(','.join(l.split()))
	data = np.loadtxt(new_del, delimiter=',', skiprows=0)

fig, ax = plt.subplots()
ax.set_xlabel(r"$t^{''} / t$", fontsize=18)
ax.set_ylabel(r"$R_{CDW}$", fontsize=18)
#props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
#ax.text(0.05, 0.95, r"$L=6$", transform=ax.transAxes, fontsize=18, verticalalignment='top', bbox=props)

L_list = [ 6, 9, 12 ]
V_list = [ 0.1, 0.3, 0.5 ]

cnt_V, cnt_L = 0, 0
for V in V_list:
	cnt_L = 0
	for L in L_list:
		data_filter = np.array(list(filter(lambda x : x[1] == V and x[0] == L, data)))
		
		if data_filter.size > 0:
			tprime = data_filter[:,2]
			suscept = data_filter[:,3]
			suscept_sigma = data_filter[:,4]

			ax.plot( tprime, suscept, color=color_cycle[cnt_V], marker=marker_cycle[cnt_L], markersize=10.0, linewidth=0.0, label=f"$V={V}, L={L}$")
			(_, caps, _) = ax.errorbar(tprime, suscept, yerr=suscept_sigma, marker='None', capsize=10, color=color_cycle[cnt_V], linewidth=3.0)
			for cap in caps:
				cap.set_markeredgewidth(2.0)
			cnt_L += 1
		cnt_V += 1

plt.legend(borderpad=0.05, labelspacing=0.075)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.savefig("pdf/R_of_tprime.pdf", bbox_inches='tight', pad_inches = 0.1)

plt.show()
