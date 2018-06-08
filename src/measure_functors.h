#pragma once
#include <ostream>
#include <vector>
#include <cmath>
#include "measurements.h"
#include "parser.h"
#include "parameters.h"
#include "green_function.h"

void eval_B_cdw(double& out,
	std::vector<std::valarray<double>*>& o)
{
	out = (*o[1])[0] / ((*o[0])[0] * (*o[0])[0]);
}

void eval_R_cdw(double& out,
	std::vector<std::valarray<double>*>& o)
{
	out = 1. - (*o[1])[0] / (*o[0])[0];
}

void eval_B_chern(double& out,
	std::vector<std::valarray<double>*>& o)
{
	out = (*o[1])[0] / ((*o[0])[0] * (*o[0])[0]);
}

void eval_R_chern(double& out,
	std::vector<std::valarray<double>*>& o)
{
	out = 1. - (*o[1])[0] / (*o[0])[0];
}

void eval_epsilon(std::valarray<double>& out,
	std::vector<std::valarray<double>*>& o)
{
	std::valarray<double>* ep_tau = o[0];
	double epsilon = (*o[1])[0];
	out.resize(ep_tau->size());
	for (int i = 0; i < ep_tau->size(); ++i)
		out[i] = (*ep_tau)[i] - epsilon * epsilon;
}

void eval_log_ratio(std::valarray<double>& out,
	std::vector<std::valarray<double>*>& o)
{
	int N = 1;
	std::valarray<double>* c_tau = o[0];
	out.resize(c_tau->size() - N);
	for (int i = 0; i < c_tau->size() - N; ++i)
		out[i] = std::log((*c_tau)[i] / (*c_tau)[i+N]);
}

void eval_n(double& out,
	std::vector<std::valarray<double>*>& o)
{
	double sign_re = (*o[0])[0];
	double sign_im = (*o[1])[0];
	double n_re = (*o[2])[0];
	double n_im = (*o[3])[0];
	out = (n_re * sign_re + n_im * sign_im)
		/ (sign_re * sign_re + sign_im * sign_im);
}

void eval_sign(double& out,
	std::vector<std::valarray<double>*>& o)
{
	double sign_re = (*o[0])[0];
	double sign_im = (*o[1])[0];
	out = std::sqrt(sign_re * sign_re + sign_im * sign_im);
}

void eval_fid_suscept(double& out, std::vector<std::valarray<double>*>& o, double* p)
{
	out = ((*o[0])[0] - (*o[1])[0] * (*o[2])[0]) / (p[0] * p[0]);
}

struct measure_M
{
	measurements& measure;
	parser& pars;
	parameters& param;
	green_function& gf;

	void init()
	{
		measure.add_observable("avg_norm_error", param.n_prebin);
		measure.add_observable("pert_order", param.n_prebin);
		measure.add_observable("k_L k_R", param.n_prebin);
		measure.add_observable("k_L", param.n_prebin);
		measure.add_observable("k_R", param.n_prebin);
		//measure.add_vectorobservable("dyn_Hv_tau", param.theta / param.block_size, param.n_prebin);
	}
	
	void perform()
	{
		if (std::abs(gf.tau() - param.theta/2.+param.block_size/2) < 1E-6)
		//if (std::abs(gf.tau() - param.theta/2.+param.block_size/2) < param.measure_window/2.)
		{
			unsigned k = gf.pert_order(), k_L = gf.pert_order(param.theta/2.), k_R = k - k_L;
			measure.add("pert_order", k);
			measure.add("k_L k_R", k_L * k_R);
			measure.add("k_L", k_L);
			measure.add("k_R", k_R);
			measure.add("avg_norm_error", gf.reset_norm_error());
			//std::vector<double> hv_tau = gf.measure_Hv_tau();
			//measure.add("dyn_Hv_tau", hv_tau);
		}
	}

	void collect(std::ostream& os)
	{
		double eval_param[] = {param.V};
		measure.add_evalable("fidelity susceptibility", "k_L k_R", "k_L", "k_R", eval_fid_suscept, eval_param);
		if (contains("M2") && contains("M4"))
			measure.add_evalable("B_cdw", "M2", "M4", eval_B_cdw);
		if (contains("M2") && contains("S_cdw_q"))
			measure.add_evalable("R_cdw", "M2", "S_cdw_q", eval_R_cdw);
		
		os << "PARAMETERS" << std::endl;
		pars.get_all(os);
		measure.get_statistics(os);
	}
	
	bool contains(const std::string& name)
	{
		return std::find(param.static_obs.begin(), param.static_obs.end(), name) != param.static_obs.end();
	}
};
