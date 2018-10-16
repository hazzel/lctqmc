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

struct wick_chi_cdw
{
	Random& rng;
	parameters& param;
	lattice& lat;
	std::vector<double> values;

	wick_chi_cdw(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(4);
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
		const int N = lat.n_sites(), Ns = N*N;
		std::fill(values.begin(), values.end(), 0.);
		for (int i = 0; i < N; i+=2)
			for (int j = 0; j < N; j+=2)
				for (int u = 0; u < 2; ++u)
					for (int v = 0; v < 2; ++v)
						values[u*2+v] += ca_td_gf[(j+v) * N + i+u] * ca_td_gf[(j+v) * N + i+u] / Ns;
		return values;
	}
};

struct wick_sp_site
{
	Random& rng;
	parameters& param;
	lattice& lat;
	std::vector<double> values;

	wick_sp_site(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(lat.max_distance());
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_td_gf = td_gf.data();
		const int N = lat.n_sites();
		std::fill(values.begin(), values.end(), 0.);
		for (int d = 0; d < lat.max_distance(); ++d)
		{
			auto& d_bonds = lat.bonds("d" + std::to_string(d) + "_bonds");
			for (auto& b : d_bonds)
			{
				int i = b.first, j = b.second;
				auto& r_i = lat.real_space_coord(i);
				auto& r_j = lat.real_space_coord(j);
				values[d] += ca_td_gf[j*N+i] * lat.parity(i)*lat.parity(j);
			}
			values[d] /= d_bonds.size();
		}
		return values;
	}
};

struct wick_sp_q
{
	Random& rng;
	parameters& param;
	lattice& lat;
	std::vector<double> values;
	static constexpr int nq = 4;

	wick_sp_q(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(nq+1);
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_td_gf = td_gf.data();
		numeric_t sp = 0.;
		const int N = lat.n_sites();
		std::fill(values.begin(), values.end(), 0.);
		for (int q = 0; q <= nq; ++q)
		{
			auto K = lat.symmetry_point("K") + q * lat.b1 / lat.Lx;
			for (int i = 0; i < N; i+=2)
			{
				auto& r_i = lat.real_space_coord(i);
				for (int j = 0; j < N; j+=2)
				{
					auto& r_j = lat.real_space_coord(j);
					double kdot = K.dot(r_i - r_j);
					values[q] += std::cos(kdot) * (ca_td_gf[(j+0)*N+(i+0)] + ca_td_gf[(j+1)*N+(i+1)] + ca_td_gf[(j+0)*N+(i+1)] + ca_td_gf[(j+1)*N+(i+0)]);
				}
			}
			values[q] /= N;
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
