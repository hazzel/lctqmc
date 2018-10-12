#pragma once
#include <vector>
#include <functional>
#include <utility>
#include <memory>
#include <ostream>
#include <iostream>
#include <boost/multi_array.hpp>
#include <Eigen/Dense>
#include "measurements.h"
#include "parameters.h"
#include "lattice.h"

typedef green_function::matrix_t matrix_t;
typedef green_function::numeric_t numeric_t;

struct wick_static_chi_cdw
{
	Random& rng;
	parameters& param;
	lattice& lat;
	std::vector<double> values;

	wick_static_chi_cdw(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(4);
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf.data();
		const int N = lat.n_sites(), Ns = N*N;
		std::fill(values.begin(), values.end(), 0.);
		for (int i = 0; i < N; i+=2)
			for (int j = 0; j < N; j+=2)
				for (int u = 0; u < 2; ++u)
					for (int v = 0; v < 2; ++v)
						values[u*2+v] += ca_et_gf_0[(j+v) * N + i+u] * ca_et_gf_0[(j+v) * N + i+u] / Ns;
		return values;
	}
};
