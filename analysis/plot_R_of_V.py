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

ecolor = ['#cc2909', '#efc600', '#60af92', '#3e77af','#3e77af','#60af92','#efc600','#cc2909']
fcolor = ['#ea6868', '#eddea2', '#99d1b9', '#a3c1e5','#a3c1e5','#99d1b9','#eddea2','#ea6868']
marker = ['o','s','D','^','o','s','D','^']
fillstyle = ['none', 'full']

filename = glob.glob("/net/home/lxtsfs1/tpc/hesselmann/code/lctqmc/plot/tprime=0.5/scan/R_cdw.txt")[0]

with open(filename) as f:
	lines = (line for line in f if not line.startswith('L'))
	new_del = []
	for l in lines:
		new_del.append(','.join(l.split()))
	data = np.loadtxt(new_del, delimiter=',', skiprows=0)

cnt = 0
fig, ax = plt.subplots()
ax.set_xlabel(r"$V/t$", fontsize=18)
ax.set_ylabel(r"$\chi_{CDW}$", fontsize=18)

L_list = [ 6, 9 ]
tprime_list = [ 0.6 ]

for tprime in tprime_list:
	for L in L_list:
		data_filter = np.array(list(filter(lambda x : x[2] == tprime and x[0] == L, data)))
		
		if data_filter.size > 0:
			V = data_filter[:,1]
			suscept = data_filter[:,3]
			suscept_sigma = data_filter[:,4]

			ax.plot( V, suscept, color=ecolor[cnt], marker=marker[cnt], markersize=10.0, linewidth=0.0, label=f"$tprime={tprime}, L={L}$")
			(_, caps, _) = ax.errorbar(V, suscept, yerr=suscept_sigma, marker='None', capsize=10, color=ecolor[cnt], linewidth=3.0)
			for cap in caps:
				cap.set_markeredgewidth(2.0)
			cnt += 1

plt.legend(borderpad=0.05, labelspacing=0.075)
leg = plt.legend()
leg.get_frame().set_linewidth(0.0)
plt.tight_layout()
plt.savefig("pdf/suscept_of_tprime.pdf", bbox_inches='tight', pad_inches = 0.1)

plt.show()
