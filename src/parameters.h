#pragma once
#include <complex>
#include <string>

struct parameters
{
	double theta, block_size, V, t, tprime, stag_mu, measure_window, dyn_tau_max, dyn_delta_tau, inv_symmetry, epsilon;
	int L, n_updates_per_block, rebuild_interval, wrap_refresh_interval, n_prebin, n_sweeps, dyn_tau_steps, direction, static_measure_interval;
	int static_measure_cnt, rebuild_cnt, wrap_refresh_cnt;
	std::vector<std::string> dyn_obs, static_obs;
	double sign_phase=1.;
};
