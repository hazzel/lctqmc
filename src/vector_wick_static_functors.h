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
	static constexpr int nq = 4;

	wick_static_chi_cdw(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(nq+1);
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf.data();
		const int N = lat.n_sites();
		std::fill(values.begin(), values.end(), 0.);
		for (int q = 0; q <= nq; ++q)
		{
			auto gamma_q = q * lat.b1 / lat.Lx;
			for (int i = 0; i < N; ++i)
			{
				auto& r_i = lat.real_space_coord((i/2)*2);
				for (int j = 0; j < N; ++j)
				{
					auto& r_j = lat.real_space_coord((j/2)*2);
					double qr = gamma_q.dot(r_i - r_j);
					values[q] += std::cos(qr) * std::real(ca_et_gf_0[j * N + i] * ca_et_gf_0[j * N + i]);
				}
			}
			values[q] /= N * N;
		}
		return values;
	}
};
