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

def FitFunctionL(x, a, b, c):
	return b*np.exp(-c*x)

def FitFunctionR(x, a, b, c):
	return a + b*np.exp(c*x)

def ExpSumFunction(x, a, b, c, d):
	return a + b*np.exp(-d*x) + c*np.exp(d*x)

def LinearFunction(x, a, b):
	return a - b*x

def parse_ed_file(filename):
	ed_data = []
	ed_lines = ed_file.read().splitlines()
	for line in ed_lines:
		ed_data.append([])
		values = filter(None, line.split("\t"))
		for x in values:
			try:
				num = float(x)
			except ValueError:
				num = 0
			ed_data[-1].append(num)
		ed_data[-1] = np.array(ed_data[-1])
	return ed_data

#latexify()
color_cycle = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'orange', 'darkgreen']
marker_cycle = ['o', 'D', '<', 'p', '>', 'v', '*', '^', 's']

filelist = []

#filelist.append(glob.glob("../job/*.out"))
filelist.append(glob.glob("/scratch/work/hesselmann/lctqmc/job/*.out"))
#filelist.append(glob.glob("/scratch/work/hesselmann/lctqmc/cluster/lctqmc_L7_theta40/*.out"))

filelist[0].sort()

filelist = [item for sublist in filelist for item in sublist]
figure, ax = plt.subplots(1, 1)
cnt = 0
for f in filelist:
	#figure, (ax1, ax, ax3) = plt.subplots(1, 3)
	plist = ParseParameters(f)
	elist = ParseEvalables(f)

	obs = "sp"
	if obs == "M2":
		ed_n = 1
		ax.set_ylabel(r"$\left \langle O_{cdw}(\tau) O_{cdw}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "epsilon_V":
		ed_n = 2
		ax.set_ylabel(r"$\left \langle O_{\epsilon_V}(\tau) O_{\epsilon_V}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "epsilon_as":
		ed_n = 3
		ax.set_ylabel(r"$\left \langle O_{\epsilon_as}(\tau) O_{\epsilon_as}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "kekule_s":
		ed_n = 4
		ax.set_ylabel(r"$\left \langle O_{kek_s}(\tau) O_{kek_s}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "kekule_as":
		ed_n = 5
		ax.set_ylabel(r"$\left \langle O_{kek_as}(\tau) O_{kek_as}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "kekule_K":
		ed_n = 6
		ax.set_ylabel(r"$\left \langle O_{kek_K}(\tau) O_{kek_K}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "chern":
		ed_n = 7
		ax.set_ylabel(r"$\left \langle O_{chern}(\tau) O_{chern}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "gamma_mod":
		ed_n = 8
		ax.set_ylabel(r"$\left \langle O_{phase}(\tau) O_{phase}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "gamma_mod_as":
		ed_n = 9
		ax.set_ylabel(r"$\left \langle O_{phase}(\tau) O_{phase}^{\dag} \right \rangle$", fontsize=16)
	elif obs == "epsilon":
		ed_n = 10
		#ax.set_ylabel(r"$\left \langle O_{\epsilon}(\tau) O_{\epsilon}^{\dag} \right \rangle - \left \langle O_{\epsilon}\right \rangle^2$", fontsize=16)
		ax.set_ylabel(r"$\left \langle O_{\epsilon}(\tau) O_{\epsilon}^{\dag} \right \rangle$", fontsize=16)
	elif "sp" in obs:
		ed_n = 11
		ax.set_ylabel(r"$\left \langle O_{sp}(\tau) O_{sp}^{\dag} \right \rangle$", fontsize=16)
	elif "tp" in obs:
		ed_n = 12
		ax.set_ylabel(r"$\left \langle O_{tp}(\tau) O_{tp}^{\dag} \right \rangle$", fontsize=16)
		
	for i in range(len(plist)):
		h = float(plist[i]["V"])
		theta = float(plist[i]["theta"])
		L = float(plist[i]["L"])
		dtau = float(plist[i]["dyn_delta_tau"])
		P = plist[i]["inv_symmetry"]
		
		#if h > 1.:
		#	continue
		
		ed_glob = glob.glob("../../ctint/data/ed_rhom_" + "*_Lx_" + str(int(L)) + "*V_" + format(h, '.6f') + "*GS*__c")
		if len(ed_glob) == 0:
			ed_glob = glob.glob("../../ctint/data/ed*_rhom_" + "*_Lx_" + str(int(L)) + "*V_" + format(h, '.6f') + "*GS*__gc")
		#ed_glob = glob.glob("../../ctint/data/ed*" + "L_" + str(int(L)) + "*V_" + format(h, '.6f') + "*T_" + format(0.01, '.6f') + "*")
		
		if len(ed_glob) > 0:
			ed_file = open(ed_glob[0])
			ed_data = parse_ed_file(ed_file)
			n_ed_tau = int(ed_data[0][9]) 
			n_ed_mat = int(ed_data[0][10])
			ed_tau = np.linspace(0., n_ed_tau * 0.2, n_ed_tau + 1)
		

		figure.suptitle(r"$L = " + str(L) + ", V = " + str(h) + ", 2 \Theta = " + str(theta) + "$", fontsize=16)# + str(1./T/2.) + "$")
		
		n_discrete_tau = float(plist[i]["dyn_tau_max"]) / dtau
		
		x_tau = np.linspace(0, theta, n_discrete_tau+1)
		y_tau = np.array(ArrangePlot(elist[i], "dyn_"+obs+"_tau")[0])
		err_tau = np.array(ArrangePlot(elist[i], "dyn_"+obs+"_tau")[1])

		if obs == "epsilon":
			y_tau = np.abs((y_tau[numpy.isfinite(y_tau)] - (ArrangePlot(elist[i], "epsilon")[0][0])**2.))
			err_tau = np.sqrt(err_tau[numpy.isfinite(err_tau)]**2. + (2.*ArrangePlot(elist[i], "epsilon")[0][0]*ArrangePlot(elist[i], "epsilon")[1][0]**2.)**2.)
			
			#y_tau = np.abs(y_tau[numpy.isfinite(y_tau)])
			#err_tau = err_tau[numpy.isfinite(err_tau)]
			#e = ArrangePlot(elist[i], "epsilon")[0][0]
			#s = ArrangePlot(elist[i], "epsilon")[1][0]
			#plt.axhline(e*e, color='r', linewidth=2.0, linestyle='--')
			#plt.axhline(e*e+2.*s*e, color='r', linewidth=1.0, linestyle='--')
			#plt.axhline(e*e-2.*s*e, color='r', linewidth=1.0, linestyle='--')
			
			
			#y_tau = np.abs(ArrangePlot(elist[i], "dyn_epjack_tau")[0])
			#err_tau = np.abs(ArrangePlot(elist[i], "dyn_epjack_tau")[1])
		elif obs == "epsilon_V":
			#y_tau = np.abs((y_tau[numpy.isfinite(y_tau)] - (ArrangePlot(elist[i], "epsilon")[0][0])**2.))
			#err_tau = np.sqrt(err_tau[numpy.isfinite(err_tau)]**2. + (2.*ArrangePlot(elist[i], "epsilon")[0][0]*ArrangePlot(elist[i], "epsilon")[1][0]**2.)**2.)
			
			#y_tau = np.abs((y_tau[numpy.isfinite(y_tau)] - 0.23882))
			#err_tau = err_tau[numpy.isfinite(err_tau)]
			
			
			y_tau = np.abs(y_tau[numpy.isfinite(y_tau)])
			err_tau = err_tau[numpy.isfinite(err_tau)]
			e = ArrangePlot(elist[i], "epsilon_V")[0][0]
			s = ArrangePlot(elist[i], "epsilon_V")[1][0]
			plt.axhline(e*e, color='r', linewidth=2.0, linestyle='--')
			plt.axhline(e*e+2.*s*e, color='r', linewidth=1.0, linestyle='--')
			plt.axhline(e*e-2.*s*e, color='r', linewidth=1.0, linestyle='--')
		elif obs == "sp_mat":
			y_tau = (np.array(ArrangePlot(elist[i], "dyn_sp_mat_0_tau")[0]) + np.array(ArrangePlot(elist[i], "dyn_sp_mat_3_tau")[0]))/2.
			err_tau = np.sqrt(np.array(ArrangePlot(elist[i], "dyn_sp_mat_0_tau")[1])**2. + np.array(ArrangePlot(elist[i], "dyn_sp_mat_3_tau")[1])**2.)/2.
		elif obs == "tp_mat":
			y_tau = (np.array(ArrangePlot(elist[i], "dyn_tp_mat_0_tau")[0]) + np.array(ArrangePlot(elist[i], "dyn_tp_mat_3_tau")[0]))/2.
			err_tau = np.sqrt(np.array(ArrangePlot(elist[i], "dyn_tp_mat_0_tau")[1])**2. + np.array(ArrangePlot(elist[i], "dyn_tp_mat_3_tau")[1])**2.)/2.
		elif obs == "Hv":
			y_tau = np.abs(y_tau[numpy.isfinite(y_tau)])
			err_tau = err_tau[numpy.isfinite(err_tau)]

			e = ArrangePlot(elist[i], "epsilon")[0][0]
			s = ArrangePlot(elist[i], "epsilon")[1][0]
			plt.axhline(e*e, color='r', linewidth=2.0, linestyle='--')
			plt.axhline(e*e+2.*s*e, color='r', linewidth=1.0, linestyle='--')
			plt.axhline(e*e-2.*s*e, color='r', linewidth=1.0, linestyle='--')
		else:
			y_tau = np.abs(y_tau[numpy.isfinite(y_tau)])
			err_tau = err_tau[numpy.isfinite(err_tau)]
		
		x_tau = np.array(range(0, len(y_tau))) * dtau
		y_tau = y_tau[:len(x_tau)]
		err_tau = err_tau[:len(x_tau)]
		
		'''
		y_log = y_tau
		err_log = err_tau
		for i in range(len(y_tau) - 1):
			y_log[i] = np.log(y_tau[i] / y_tau[i+1]) / dtau
			err_log[i] = ((err_tau[i]/y_tau[i+1])**2. + (err_tau[i+1]/y_tau[i])**2.)**0.5 / dtau
		x_tau = x_tau[:len(x_tau)-1]
		y_tau = y_log[:len(x_tau)]
		err_tau = err_log[:len(x_tau)]
		'''
		
		ax.set_xlabel(r"$\tau$", fontsize=16)
		ax.set_yscale("log")
		
		if cnt == 0:
			ax.plot(x_tau, y_tau, marker="o", color="green", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+', P = ' + P + '$')
			(_, caps, _) = ax.errorbar(x_tau, y_tau, yerr=err_tau, marker='None', capsize=8, color="green")
			for cap in caps:
				cap.set_markeredgewidth(1.6)
		else:
			ax.plot(x_tau, y_tau, marker="o", color="red", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+', P = ' + P + '$')
			(_, caps, _) = ax.errorbar(x_tau, y_tau, yerr=err_tau, marker='None', capsize=8, color="red")
			for cap in caps:
				cap.set_markeredgewidth(1.6)

		
		if len(ed_glob) > 0 and len(ed_data) > ed_n:
			ax.plot(ed_tau, ed_data[ed_n], marker='o', color="b", markersize=10.0, linewidth=0.0, label=r'$L='+str(int(L))+'$')
			#ax.plot(ed_tau, np.flipud(ed_data[ed_n]), marker='o', color="b", markersize=10.0, linewidth=2.0, label=r'$L='+str(int(L))+'$')
		
		'''
		j = 1
		#f_min = 150; f_max = 300
		f_min = 0
		#f_max = 15
		f_max = len(x_tau)-1
		step = 1
		fit_x = []
		fit_y = []
		fit_e = []
		fit_re = []
		while f_min + j*step <= f_max:
			nmin = f_min + (j-1) * step
			#nmax = f_min + j * step
			nmax = f_max
			
			try:
				#parameter, perr = fit_function( [0.1, 0.1, 1.], x_tau[nmin:nmax], y_tau[nmin:nmax], FitFunctionL, datayerrors=err_tau[nmin:nmax])
				#fit_x.append(nmin)
				#fit_y.append(parameter[2]*(2.*L*L)**0.5)
				#fit_e.append(perr[2]*(2.*L*L)**0.5)
				#fit_re.append(perr[2] / parameter[2]*(2.*L*L)**0.5)
				
				parameter, perr = scipy.optimize.curve_fit( FitFunctionL, x_tau[nmin:nmax], y_tau[nmin:nmax], p0=[0.1, 0.1, 1.], method='trf')
				fit_x.append(nmin)
				fit_y.append(parameter[2]*(2.*L*L)**0.5)
				fit_e.append(np.abs(np.sqrt(perr[2,2]))*(2.*L*L)**0.5)
				fit_re.append(np.abs(np.sqrt(perr[2,2]) / parameter[2])*(2.*L*L)**0.5)
			except RuntimeError:
				print "run time error during fit"
		
			px = np.linspace(x_tau[nmin], x_tau[nmax], 1000)
			#ax.plot(px, FitFunctionL(px, *parameter), 'k-', linewidth=3.0)
			#print "V = " + str(h)
			#print parameter
			#print perr
			#print parameter*(2.*L*L)**0.5
			#print perr*(2.*L*L)**0.5
			
			j+=1
		ax.axvline(f_min*dtau, color='k', linestyle='--')
		ax.axvline(f_max*dtau, color='k', linestyle='--')
		ax.legend()
		
		#f2, ax2 = plt.subplots(1, 1)
		#ax2.set_xlabel(r"$\tau_{min}$", fontsize=16)
		#ax2.set_ylabel(r"$\Delta_{sp} \sqrt{N}$", fontsize=16)
		#ax2.plot(fit_x, fit_y, marker='o', markersize=8.0, color="k")
		#(_, caps, _) = ax2.errorbar(fit_x, fit_y, yerr=fit_e, marker='None', capsize=8, color="k")
		#for cap in caps:
		#	cap.set_markeredgewidth(1.6)
		
		
		if len(fit_re) > 0:
			#ax2.axvline(fit_x[fit_re.index(min(np.abs(fit_re)))])
			#ax2.axhline(fit_y[fit_re.index(min(np.abs(fit_re)))])
			#f3, ax3 = plt.subplots(1, 1)
			#ax3.set_xlabel(r"$\tau_{min}$", fontsize=16)
			#ax3.set_ylabel(r"$\sigma / (\Delta / \sqrt{N})$", fontsize=16)
			#ax3.plot(fit_x, fit_re, marker='o', markersize=8.0, color="k")
			#ax3.axvline(fit_x[fit_re.index(min(np.abs(fit_re)))])
			#ax3.axhline(min(np.abs(fit_re)))

			
			#min = 0; nmax = 2*int(plist[i]["discrete_tau"])
			nmin = fit_x[fit_re.index(min(np.abs(fit_re)))]; nmax = f_max
			#nmin = 25; nmax = 50
			#parameter, perr = fit_function( [1., 6., 1.2], x_tau[nmin:nmax], y_tau[nmin:nmax], FitFunctionL, datayerrors=err_tau[nmin:nmax])
			parameter, perr = fit_function( [0.1, 0.1, 1.], x_tau[nmin:nmax], y_tau[nmin:nmax], FitFunctionL, datayerrors=err_tau[nmin:nmax])
			#parameter, perr = scipy.optimize.curve_fit( FitFunctionL, x_tau[nmin:nmax], y_tau[nmin:nmax], p0=[0.1, 0.1, 1.])
		
			px = np.linspace(x_tau[nmin], x_tau[nmax], 1000)
			ax.plot(px, FitFunctionL(px, *parameter), 'k-', linewidth=3.0)
			
			#print "V = " + str(h)
			#print "Nmin = ", fit_x[fit_re.index(min(np.abs(fit_re)))], " Nmax = ", f_max
			#print "Delta * sqrt(N) = ", fit_y[fit_re.index(min(np.abs(fit_re)))], " +- ", fit_e[fit_re.index(min(np.abs(fit_re)))]
			#print "Delta = ", fit_y[fit_re.index(min(np.abs(fit_re)))]/(2.*L*L)**0.5, " +- ", fit_e[fit_re.index(min(np.abs(fit_re)))]/(2.*L*L)**0.5
			#print "Delta * sqrt(N) = ", parameter[2] * (2.*L*L)**0.5, " +- ", perr[2]*(2.*L*L)**0.5
			#print "Delta = ", parameter[2], " +- ", perr[2]
			#print "------"
			
			#print str(int(L)) + "\t" + str(h) + "\t\t" + str(round(parameter[2] * (2.*L*L)**0.5, 5)) + "\t\t" + str(round(perr[2] * (2.*L*L)**0.5, 2)) + "\t\t\t0\t\t0"
			print str(int(L)) + "\t" + str(h) + "\t\t" + str(round(parameter[2], 5)) + "\t\t" + str(round(perr[2], 2)) + "\t\t\t0\t\t0"
		'''
		
		
		nmin = 4
		nmax = 50
		
		try:
			parameter, perr = scipy.optimize.curve_fit( FitFunctionL, x_tau[nmin:nmax], y_tau[nmin:nmax], p0=[1., 0.0004, 0.8], method='trf')
			px = np.linspace(x_tau[nmin], x_tau[nmax], 1000)
			ax.plot(px, FitFunctionL(px, *parameter), 'k-', linewidth=3.0)
			
			print "V = " + str(h)
			print parameter
			print parameter[2] * (2.*L*L)**0.5, " +- ", np.sqrt(perr[2,2]) * (2.*L*L)**0.5
			print parameter[2], " +- ", np.sqrt(perr[2,2])
		except RuntimeError:
			print "run time error during fit"
		
		
		if len(ed_glob) > 0 and len(ed_data) > ed_n:
			nmin = len(ed_tau) / 2
			nmax = len(ed_tau) - 1
			#nmin = 1
			#nmax = 6
			parameter_ed, perr_ed = scipy.optimize.curve_fit( FitFunctionL, ed_tau[nmin:nmax], ed_data[ed_n][nmin:nmax], p0=[0.1, 0.1, 1.])
			px = np.linspace(ed_tau[nmin], ed_tau[nmax], 1000)
			#px = np.linspace(ed_tau[len(ed_data[ed_n])/2], ed_tau[len(ed_data[ed_n])-1], 1000)
			ax.plot(px, FitFunctionL(px, *parameter_ed), 'k-', linewidth=3.0)
			print parameter_ed
		
		cnt += 1

plt.legend()
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.show()
