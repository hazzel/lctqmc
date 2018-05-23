#pragma once
#include <complex>
#include <string>

struct parameters
{
	double theta, block_size, V, t, tprime, stag_mu;
	int L, n_updates_per_block, n_prebin;
	std::vector<std::string> obs, static_obs;
	double sign_phase=1.;
};
