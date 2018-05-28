#pragma once
#include <vector>
#include <functional>
#include <utility>
#include <memory>
#include <ostream>
#include <iostream>
#include <boost/multi_array.hpp>
#include "measurements.h"
#include "parameters.h"
#include "lattice.h"

typedef green_function::matrix_t matrix_t;
typedef green_function::numeric_t numeric_t;

// M2(tau) = sum_ij <(n_i(tau) - 1/2)(n_j - 1/2)>
struct wick_M2
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_M2(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
		numeric_t M2 = 0.;
		const int N = lat.n_sites();
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
			{
				/*
				M2 += config.l.parity(i) * config.l.parity(j)
					* std::real((1. - et_gf_t(i, i)) * (1. - et_gf_0(j, j))
					+ config.l.parity(i) * config.l.parity(j) * td_gf(i, j) * td_gf(i, j)
					- (et_gf_t(i, i) + et_gf_0(j, j))/2. + 1./4.);
				*/
				//M2 += td_gf(i, j) * td_gf(i, j);
				M2 += ca_td_gf[j * N + i] * ca_td_gf[j * N + i];
			}
		return std::real(M2) / std::pow(N, 2.);
	}
};

// sp(tau) = sum_ij e^{-i K (r_i - r_j)} <c_i(tau) c_j^dag>
struct wick_sp
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_sp(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_td_gf = td_gf.data();
		numeric_t sp = 0.;
		auto& K = lat.symmetry_point("K");
		const int N = lat.n_sites();
		for (int o = 0; o < 2; ++o)
		for (int i = 0; i < N; i+=2)
		{
			auto& r_i = lat.real_space_coord(i);
			for (int j = 0; j < N; j+=2)
			{
				auto& r_j = lat.real_space_coord(j);
				double kdot = K.dot(r_i - r_j);
				sp += std::cos(kdot) * ca_td_gf[(j+o)*N+(i+o)] * lat.parity(i+o)*lat.parity(j+o);
				//sp += std::cos(kdot) * ca_td_gf[j*N+i] * (1. + lat.parity(i)*lat.parity(j));
			}
		}
		return std::real(sp) / N;
	}
};
