import numpy as np
import scipy as sp
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.ticker as ticker
import glob
import math

plt.rc("mathtext", fontset="cm")
plt.rc("font", family="serif")
plt.rc('legend')
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rc('text', usetex=True)

ecolor = ['#cc2909', '#efc600', '#60af92', '#3e77af','#3e77af','#60af92','#efc600','#cc2909']
fcolor = ['#ea6868', '#eddea2', '#99d1b9', '#a3c1e5','#a3c1e5','#99d1b9','#eddea2','#ea6868']
marker = ['o','s','D','^','o','s','D','^']
markerSize = np.array([5,4.5,4,5.5,5,4.5,4,5.5])*2.0

fsTicksLabel = 18
fsTicksLabelInset = 16
fsAxesLabel = 21
fsLegend = 16
fsLabel = 18

# BZ:
#  __
# /  \
# \__/
a1 = np.array([1./2., -np.sqrt(3.)/2.])
a2 = np.array([1., 0.])
b1 = np.array([2.*np.pi, -2.*np.pi/3.*np.sqrt(3.)])
b2 = np.array([0., 4.*np.pi/np.sqrt(3.)])
K = { 6 : np.array([4.1887902, 0.]), 9 : np.array([4.1887902, 0.]), 12 : np.array([4.1887902, 0.]), 15 : np.array([4.1887902, 0.]) }

i_L = 0
i_gamma = 1
i_U = 2
i_q = 3
i_kx = 4
i_ky = 5
i_delta = 6
i_sigma = 7

L_list = [ 15 ]
gamma_list = [ 0. ]
U_list = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.4, 3.6, 3.75, 4, 4.5 ]
q_list = [ 0, 1, 2]

fig = plt.figure(figsize=(11*0.7,8*0.7), dpi=96)
#------------------------------------------------------------------------------------

ax = plt.subplot(121)
ax.tick_params(which='both', axis='both', direction='in', bottom='on', top='on', left='on', right='on', labelsize=fsTicksLabel)
ax.set_xlabel(r"$U/t$", fontsize=fsAxesLabel)
ax.set_ylabel(r"$E/t$", fontsize=fsAxesLabel)
ax.set_yticks([0, 0.4, 0.8, 1.2])
ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(4))
ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax.set_xlim(0.3, 4.7)
ax.set_ylim(-0.05, 1.2)
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
ax.text(-0.31, 0.98, r'\textbf{A}', color='#000000', transform=ax.transAxes, fontSize=fsLabel, ha='left')
ax.text( 2.34, 0.98, r'\textbf{B}', color='#000000', transform=ax.transAxes, fontSize=fsLabel, ha='left')
ax.text(0.08, 0.91, r'$\gamma = 0, L=15$', color='#000000', transform=ax.transAxes, fontSize=fsLabel, ha='left', bbox=dict(boxstyle="round", ec=(219/255., 202/255., 175/255.), fc=(244/255., 232/255., 211/255.)))
ax.text(0.20, 0.75, r"$a \Delta k = 0.97$", transform=ax.transAxes, fontsize=16, verticalalignment='top')
ax.text(0.20, 0.44, r"$a \Delta k = 0.48$", transform=ax.transAxes, fontsize=16, verticalalignment='top')
ax.text(0.20, 0.12, r"$a \Delta k = 0$", transform=ax.transAxes, fontsize=16, verticalalignment='top')
ax.text(0.715, 0.985, r"$U_c/t = 3.85$", transform=ax.transAxes, fontsize=16, verticalalignment='top', color="#888888", rotation=90)

plt.axhline(0, color='#dddddd', zorder=-2, linewidth=0.5)
plt.axvline(3.85, linestyle='--', color='#888888', linewidth=1.0)

filelist = glob.glob('data_of_q.txt')
filelist.sort()

cnt_q = 3
for filename in filelist:
	with open(filename) as f:
		lines = (line for line in f if not line.startswith('L'))
		new_del = []
		for l in lines:
			new_del.append(','.join(l.split()))
		data = np.loadtxt(new_del, delimiter=',', skiprows=0)

	data_filter = list(filter(lambda x : x[i_L] in L_list and x[i_gamma] in gamma_list and x[i_U] in U_list, data))

	for L in L_list:
		data_filter_L = np.array(list(filter(lambda x : x[i_L] == L and x[i_U] in U_list, data_filter)))
		
		for q_i in q_list:
			data_filter_q = np.array(list(filter(lambda x : (np.abs(K[L]-b1/L*q_i - np.array([x[i_kx], x[i_ky]])) < 1E-6).all().any(), data_filter_L)))
			delta_list = data_filter_q[:, i_delta]
			sigma_list = data_filter_q[:, i_sigma]			
			ax.errorbar(U_list, delta_list, yerr=sigma_list, markeredgewidth=1.6, capsize=8, fmt=marker[cnt_q]+'-', markersize=markerSize[cnt_q], color=ecolor[cnt_q], ecolor=ecolor[cnt_q], mfc=fcolor[cnt_q])

#------------------------------------------------------------------------------------

L_list = [ 6, 9, 12, 15 ]
gamma_list = [ 0. ]
U_list = [ 0.5, 1., 1.5, 2., 2.5, 3., 3.4, 3.6, 3.75, 4 ]
q_list = [ 1 ]
v0 = 3.**0.5/2.

ax2 = plt.subplot(122)
ax2.tick_params(which='both', axis='both', direction='in', bottom='on', top='on', left='on', right='on', labelleft='off',  labelright='on', labelsize=fsTicksLabel)
ax2.set_xlabel(r"$U/t$", fontsize=fsAxesLabel)
ax2.set_ylabel(r"$[E / (a \Delta k) - v_0]/v_0$", fontsize=fsAxesLabel, labelpad=15)
ax2.yaxis.set_label_position('right')
ax2.set_yticks([-0.5, -0.4, -0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3])
ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax2.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))
ax2.set_xlim(0.3, 4.15)
ax2.set_ylim(-0.51, 0.11)

ax2.text(0.08, 0.91, r'$\gamma = 0$', color='#000000', transform=ax2.transAxes, fontSize=fsLabel, ha='left', bbox=dict(boxstyle="round", ec=(219/255., 202/255., 175/255.), fc=(244/255., 232/255., 211/255.)))
ax2.text(0.3, 0.12, r'40\% reduced $v_0$', color='#cc2909', transform=ax2.transAxes, fontSize=16, ha='left')
ax2.text(0.82, 0.585, r"$U_c/t = 3.85$", transform=ax2.transAxes, fontsize=16, verticalalignment='top', color="#888888", rotation=90)

plt.axhline(0, color='#dddddd', zorder=-2, linewidth=0.5)
plt.axvline(3.85, linestyle='--', color='#888888', linewidth=1.0)
plt.axhline(y=-0.4, color='#cc2909', zorder=-2, linewidth=2.)

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
			
			ax2.errorbar(U_list, delta_list, yerr=sigma_list, markeredgewidth=1.6, capsize=8, fmt=marker[cnt_q]+'-', markersize=markerSize[cnt_q], color=ecolor[cnt_q], ecolor=ecolor[cnt_q], mfc=fcolor[cnt_q], label=f"$L={L}$")
			cnt_q += 1

#plt.tight_layout()
fig.subplots_adjust(left=0.1, bottom=0.1,right=0.9,top=0.9,wspace=0.05)
leg = plt.legend(borderpad=0.15, labelspacing=0.4, loc=6, fontsize=16)
leg.get_frame().set_linewidth(0.)
plt.savefig("pdf/gap_of_U.pdf", bbox_inches='tight', pad_inches = 0.)
plt.show()
