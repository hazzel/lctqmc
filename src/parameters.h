#pragma once
#include <complex>
#include <string>

struct parameters
{
	double theta, block_size, V, t, tprime, stag_mu, measure_window, dyn_tau_max, dyn_delta_tau;
	int L, n_updates_per_block, n_prebin, n_sweeps, dyn_tau_steps, direction, static_measure_interval, static_measure_cnt;
	std::vector<std::string> dyn_obs, static_obs;
	double sign_phase=1.;
};
