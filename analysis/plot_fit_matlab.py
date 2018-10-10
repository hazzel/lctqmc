import numpy as np
import scipy.io
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import sys
sys.path.append("/net/home/lxtsfs1/tpc/hesselmann/mc/qising-SSE")
sys.path.append("/home/stephan/mc/qising-SSE")
from Fit import *

plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
plt.rcParams.update({'font.size': 20})
plt.rc('legend',fontsize=18)

def FitFunction(x, a, b):
	return a*np.exp(-b*x)
	
def fold_back_to_bz(K):
	dist = 1000
	xs = K[0]
	ys = K[1]
	xg = 0
	yg = 0
	G1x = b1[0]
	G1y = b1[1]
	G2x = b2[0]
	G2y = b2[1]

	for i in range(-1, 2):
		for j in range(-1, 2):
			x0 = i*G1x+j*G2x
			y0 = i*G1y+j*G2y

			dist2 = (xs-x0)**2. + (ys-y0)**2.
			if dist2 < dist:
				dist=dist2
				xg=xs-x0
				yg=ys-y0
	return np.array([xg, yg])

data = scipy.io.loadmat("../fahker_science_paper/datGF.mat")

i_gamma = 0
i_U = 1
i_alpha = 2
i_UoverUc = 3
i_L = 4
i_Gk = 5

j_kx = 0
j_ky = 1
j_Gc = 2
j_Gv = 3

k_tau = 0
k_g = 1
k_sigma = 2

# BZ:
#  __
# /  \
# \__/

b1 = np.array([2.*np.pi, -2.*np.pi/3.*np.sqrt(3.)])
b2 = np.array([0., 4.*np.pi/np.sqrt(3.)])
x0 = 2.0/(3.*3.**0.5) * np.pi * np.sqrt(3.)
y0 = 2.0/3.0 * np.pi * np.sqrt(3.)
corner = []
corner.append([-x0, -y0])		#down left
corner.append([x0, -y0])		#down right
corner.append([2.*x0, 0.])		#right
corner.append([x0, y0])			#up right
corner.append([-x0, y0])		#up left
corner.append([-2.*x0, 0.])	#left
corner.append([-x0, -y0])		#down left
corner = np.array(corner)

K = { 6 : np.array([4.1887902, 0.]), 9 : np.array([4.1887902, 0.]), 12 : np.array([4.1887902, 0.]), 15 : np.array([4.1887902, 0.]) }

gamma_list = [ 0. ]
L_list = [ 15 ]
#U_list = [ 0.5, 1., 1.5, 2., 3., 3.6, 3.75, 4, 4.5 ]
#L_list = [ 12 ]
U_list = [ 3.6 ]

color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']


### PLOT G_K(tau)

fig, ax = plt.subplots()
ax.set_xlabel(r"$\tau$", fontsize=16)
ax.set_ylabel(r"$G(K, \tau)$", fontsize=16)
ax.set_yscale("log")
#ax.set_xlim(0., 5.)
#ax.set_ylim(0.005, 1.05)

data_filter = list(filter(lambda x : x[i_gamma][0][0] in gamma_list and x[i_L][0][0] in L_list and x[i_U][0][0] in U_list, data["datGF"][1:]))
#data_filter = data["datGF"][1:]

cnt = 0
for data_G in data_filter:
	gamma = data_G[i_gamma][0][0]
	U = data_G[i_U][0][0]
	UoverUc = data_G[i_UoverUc][0][0]
	L = data_G[i_L][0][0]
	
	#if L == 6:
	#	print("gamma = ", gamma, ", U = ", U, ", U/Uc = ", UoverUc)
	
	output = []
	
	momenta = np.array([ K[L] - i * b1 / L for i in range(6) ])
	data_G_k_filter = list(filter(lambda x : ( np.abs(momenta - np.array([ x[j_kx][0][0], x[j_ky][0][0] ])) < 1E-6 ).all(1).any(), data_G[i_Gk][1:]))
	#data_G_k_filter = data_G[i_Gk][1:]

	for data_G_k in data_G_k_filter:
		kx = data_G_k[j_kx][0][0]
		ky = data_G_k[j_ky][0][0]
		q_K = np.linalg.norm( K[L] - np.array([ kx, ky ]) )
		q_Kp = min(np.linalg.norm( np.array([x0, y0]) - np.array([ kx, ky ]) ), np.linalg.norm( np.array([x0, -y0]) - np.array([ kx, ky ]) ))
		if q_K > 4. or q_K > q_Kp:
			continue
		tau = data_G_k[j_Gv][:,0]
		Gp = np.abs(data_G_k[j_Gv][:,1])# + np.abs(data_G_k[j_Gc][:,1])
		Gp_sigma = data_G_k[j_Gv][:,2]# + data_G_k[j_Gc][:,2]

		ax.plot(tau, Gp, marker="o", color=color_cycle[cnt%len(color_cycle)], markersize=10.0, linewidth=2.0)
		(_, caps, _) = ax.errorbar(tau, Gp, yerr=Gp_sigma, marker='None', capsize=8, color=color_cycle[cnt%len(color_cycle)])
		for cap in caps:
			cap.set_markeredgewidth(1.6)
		cnt += 1
			
		nmin = 20; nmax = 60
		parameter, perr = fit_function( [5., 0.5], tau[nmin:nmax], Gp[nmin:nmax], FitFunction, datayerrors=Gp_sigma[nmin:nmax])
		#parameter, perr = scipy.optimize.curve_fit( FitFunction, tau[nmin:nmax], Gp[nmin:nmax], p0=[5., 0.5], method='trf')

		px = np.linspace(tau[nmin], tau[nmax], 1000)
		ax.plot(px, FitFunction(px, *parameter), 'k-', linewidth=3.0)

		output.append([q_K, kx, ky, parameter[1], perr[1]])
	output = np.array(output)
	output = output[np.argsort(output[:,0])]
	for k in output:
		print(str(int(L)) + "\t" + str(gamma) + "\t" + str(U) + "\t" + str(UoverUc) + "\t" + str(round(k[0], 5)) + "\t" + str(k[1]) + "\t\t" + str(k[2]) + "\t\t" + str(round(k[3], 5)) + "\t\t\t\t\t" + str(round(k[4], 5)))
	print("")
	

'''
### PLOT RECIPROCAL LATTICE

L = 15
data_filter = list(filter(lambda x : x[i_gamma][0][0] == 0. and x[i_U][0][0] == 0.1 and x[i_L][0][0] == L, data["datGF"]))

x = []
y = []
for kx, ky in zip(data_filter[0][i_Gk][1:,j_kx], data_filter[0][i_Gk][1:,j_ky]):
	folded_back_k = fold_back_to_bz([kx[0][0], ky[0][0]])
	#x.append(folded_back_k[0])
	#y.append(folded_back_k[1])
	x.append(kx[0][0])
	y.append(ky[0][0])

fig, ax = plt.subplots()
plt.title(r"$L = " + str(L) + "$")
ax.set_xlabel(r"$k_x$", fontsize=16)
ax.set_ylabel(r"$k_y$", fontsize=16)
plt.xlim(-5, 8)
plt.ylim(-5, 8)
ax.scatter(x, y, s=125, c="k")
plt.plot(corner[:,0], corner[:,1], marker="o", lw=4.0, ms=15.0, fillstyle="top", c="b")
ax.plot([K[L][0]], [K[L][1]], c="r", marker="o", markersize=15, fillstyle="top", label="K")
ax.plot([K[L][0] - b1[0] / L], [K[L][1] - b1[1] / L], c="g", marker="o", markersize=15, fillstyle="top", label="Kq")
#ax.plot([M[0]], [M[1]], c="orange", marker="o", markersize=15, fillstyle="top", label="M")
#ax.plot([Gamma[0]], [Gamma[1]], c="b", marker="o", markersize=15, fillstyle="top", label="Gamma")
ax.arrow(0, 0, b1[0], b1[1], head_width=0.25, head_length=0.25, fc='k', ec='k')
ax.arrow(0, 0, b2[0], b2[1], head_width=0.25, head_length=0.25, fc='k', ec='k')
ax.arrow(b1[0], b1[1], b2[0], b2[1], head_width=0., head_length=0., fc='k', ec='k')
ax.arrow(b2[0], b2[1], b1[0], b1[1], head_width=0., head_length=0., fc='k', ec='k')
'''

plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig("matlab.pdf", bbox_inches='tight', pad_inches = 0.1)
plt.show()
