#pragma once
#include <ostream>
#include <vector>
#include <cmath>
#include "measurements.h"
#include "parser.h"

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

struct measure_M
{
	measurements& measure;
	parser& pars;

	void perform()
	{}

	void collect(std::ostream& os)
	{
		os << "PARAMETERS" << std::endl;
		pars.get_all(os);
		measure.get_statistics(os);
	}
	
	bool contains(const std::string& name)
	{
		//return std::find(config.param.static_obs.begin(), config.param.static_obs.end(), name) != config.param.static_obs.end();
		return false;
	}
};
