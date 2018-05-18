#pragma once
#include <complex>
#include <string>

struct parameters
{
	double theta, block_size, V;
	int L;
	std::vector<std::string> obs, static_obs;
	double sign_phase=1.;
};
