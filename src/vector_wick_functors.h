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

struct wick_sp_site
{
	Random& rng;
	parameters& param;
	lattice& lat;
	std::vector<double> values;

	wick_sp_site(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(lat.n_sites()*lat.n_sites()/4);
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_td_gf = td_gf.data();
		const int N = lat.n_sites();
		auto& K = lat.symmetry_point("K");
		std::fill(values.begin(), values.end(), 0.);
		for (int o = 0; o < 2; ++o)
		for (int i = 0; i < N; i+=2)
		{
			auto& r_i = lat.real_space_coord(i);
			for (int j = 0; j < N; j+=2)
			{
				auto& r_j = lat.real_space_coord(j);
				double kdot = K.dot(r_i - r_j);
				values[i/2*N/2 + j/2] += std::cos(kdot) * ca_td_gf[(j+o)*N+(i+o)] * lat.parity(i+o)*lat.parity(j+o) / 2.;
			}
		}
		return values;
	}
};


struct wick_tp_matrix
{
	Random& rng;
	parameters& param;
	lattice& lat;
	std::vector<double> values;

	wick_tp_matrix(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(4);
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_td_gf = td_gf.data();
		auto& K = lat.symmetry_point("K");
		auto& Kp = lat.symmetry_point("Kp");
		const int N = lat.n_sites();
		std::fill(values.begin(), values.end(), 0.);
		int i = 2 * static_cast<int>(rng() * N / 2);
		auto& r_i = lat.real_space_coord(i);
		for (int m = 0; m < N; m+=2)
			for (int u = 0; u < 2; ++u)
			{
				auto& r_m = lat.real_space_coord(m);
				int j = 2 * static_cast<int>(rng() * N / 2);
				auto& r_j = lat.real_space_coord(j);
				for (int n = 0; n < N; n+=2)
					for (int v = 0; v < 2; ++v)
					{
						auto& r_n = lat.real_space_coord(n);
						double kdot;
						if (param.L % 3 == 0)
							kdot = K.dot(r_i - r_j - r_m + r_n);
						else
							kdot = K.dot(r_i - r_m) + Kp.dot(r_j - r_n);
						values[u*2+v] += std::real(std::cos(kdot) * (ca_td_gf[(i+u)*N+(m+u)] * ca_td_gf[(j+v)*N+(n+v)] - ca_td_gf[(i+u)*N+(n+v)] * ca_td_gf[(j+v)*N+(m+u)]));
					}
			}
		return values;
	}
};
