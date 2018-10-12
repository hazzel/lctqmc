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
#include "green_function.h"

typedef green_function::matrix_t matrix_t;
typedef green_function::numeric_t numeric_t;

struct wick_static_energy
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_static_energy(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf)
	{
		numeric_t energy = 0.;
		for (auto& a : lat.bonds("nearest neighbors"))
			energy += param.t * et_gf(a.second, a.first);
		for (auto& a : lat.bonds("nearest neighbors"))
			energy += param.V * ((1. - et_gf(a.first, a.first)) * (1. - et_gf(a.second, a.second))
				- et_gf(a.second, a.first) * et_gf(a.first, a.second) - (et_gf(a.first, a.first) + et_gf(a.second, a.second))/2. + 1./4.)/2.;
		for (auto& a : lat.bonds("t3_bonds"))
			energy += param.tprime * et_gf(a.second, a.first);
		for (int i = 0; i < lat.n_sites(); ++i)
			energy += -lat.parity(i) * param.stag_mu * et_gf(i, i);
		return std::real(energy);
	}
};

struct wick_static_epsilon
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_static_epsilon(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf)
	{
		numeric_t epsilon = 0.;
		for (auto& a : lat.bonds("nearest neighbors"))
			epsilon += std::real(et_gf(a.second, a.first));
		return std::real(epsilon) / lat.n_bonds();
	}
};

struct wick_static_Hv
{
	Random& rng;
	parameters& param;
	lattice& lat;
	green_function& gf;

	wick_static_Hv(Random& rng_, parameters& param_, lattice& lat_, green_function& gf_)
		: rng(rng_), param(param_), lat(lat_), gf(gf_)
	{}
	
	double get_obs(const matrix_t& et_gf)
	{
		numeric_t e = 0.;
		for (auto& a : lat.bonds("single_d1_bonds"))
			e -= param.V * std::real(et_gf(a.second, a.first) * et_gf(a.first, a.second));
		return std::real(e);
	}
};

// M2(tau) = sum_ij [ eta_i eta_j <(n_i - 1/2)(n_j - 1/2)> ]
struct wick_static_M2
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_static_M2(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf.data();
		numeric_t M2 = 0.;
		const int N = lat.n_sites();
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < N; ++j)
				M2 += ca_et_gf_0[j * N + i] * ca_et_gf_0[j * N + i];
		return std::real(M2) / std::pow(N, 2.);
	}
};

// S_cdw(q) = sum_ij [ <(n_i - 1/2)(n_j - 1/2)> e^(i q (r_i - r_j)) ]
struct wick_static_S_cdw_q
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_static_S_cdw_q(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf.data();
		numeric_t S = 0.;
		auto& q = lat.symmetry_point("q");
		const int N = lat.n_sites();
		for (int i = 0; i < N; ++i)
		{
			auto& r_i = lat.real_space_coord((i/2)*2);
			for (int j = 0; j < N; ++j)
			{
				auto& r_j = lat.real_space_coord((j/2)*2);
				double qr = q.dot(r_i - r_j);
				S += ca_et_gf_0[j * N + i] * ca_et_gf_0[j * N + i] * std::cos(qr);
			}
		}
		return std::real(S) / std::pow(lat.n_sites(), 2.);
	}
};

// S_cdw(q) = sum_ij [ <(n_i - 1/2)(n_j - 1/2)> e^(i q (r_i - r_j)) ]
struct wick_static_S_cdw_q20
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_static_S_cdw_q20(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf.data();
		numeric_t S = 0.;
		auto& q = lat.symmetry_point("q20");
		const int N = lat.n_sites();
		for (int i = 0; i < N; ++i)
		{
			auto& r_i = lat.real_space_coord((i/2)*2);
			for (int j = 0; j < N; ++j)
			{
				auto& r_j = lat.real_space_coord((j/2)*2);
				double qr = q.dot(r_i - r_j);
				S += ca_et_gf_0[j * N + i] * ca_et_gf_0[j * N + i] * std::cos(qr);
			}
		}
		return std::real(S) / std::pow(lat.n_sites(), 2.);
	}
};

// S_cdw(q) = sum_ij [ <(n_i - 1/2)(n_j - 1/2)> e^(i q (r_i - r_j)) ]
struct wick_static_S_cdw_q11
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_static_S_cdw_q11(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf.data();
		numeric_t S = 0.;
		auto& q = lat.symmetry_point("q11");
		const int N = lat.n_sites();
		for (int i = 0; i < N; ++i)
		{
			auto& r_i = lat.real_space_coord((i/2)*2);
			for (int j = 0; j < N; ++j)
			{
				auto& r_j = lat.real_space_coord((j/2)*2);
				double qr = q.dot(r_i - r_j);
				S += ca_et_gf_0[j * N + i] * ca_et_gf_0[j * N + i] * std::cos(qr);
			}
		}
		return std::real(S) / std::pow(lat.n_sites(), 2.);
	}
};

// M4(tau) = sum_ij sum_kl <(n_i(tau) - 1/2)(n_j - 1/2) (n_k(tau) - 1/2)(n_l - 1/2)>
struct wick_static_M4
{
	Random& rng;
	parameters& param;
	lattice& lat;
	const numeric_t* ca_et_gf_0;

	wick_static_M4(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	numeric_t evaluate(Eigen::Matrix<numeric_t, 4, 4>& mat44, int ns, int i, int j, int k, int l)
	{
		const double pi = lat.parity(i), pj = lat.parity(j), pk = lat.parity(k), pl = lat.parity(l);
		mat44(0, 1) = ca_et_gf_0[j*ns+i];
		mat44(1, 0) = -pi * pj * ca_et_gf_0[j*ns+i];
		mat44(0, 2) = ca_et_gf_0[k*ns+i];
		mat44(2, 0) = -pi * pk * ca_et_gf_0[k*ns+i];
		mat44(0, 3) = ca_et_gf_0[l*ns+i];
		mat44(3, 0) = -pi * pl * ca_et_gf_0[l*ns+i];
		mat44(1, 2) = ca_et_gf_0[k*ns+j];
		mat44(2, 1) = -pj * pk * ca_et_gf_0[k*ns+j];
		mat44(1, 3) = ca_et_gf_0[l*ns+j];
		mat44(3, 1) = -pj * pl * ca_et_gf_0[l*ns+j];
		mat44(2, 3) = ca_et_gf_0[l*ns+k];
		mat44(3, 2) = -pk * pl * ca_et_gf_0[l*ns+k];
		
		return pi * pj * pk * pl * mat44.determinant();
	}
	
	double get_obs(const matrix_t& et_gf)
	{
		ca_et_gf_0 = et_gf.data();
		numeric_t M4 = 0.;
		const int n = lat.n_sites();
		Eigen::Matrix<numeric_t, 4, 4> mat44 = Eigen::Matrix<numeric_t, 4, 4>::Zero();
		
		M4 += evaluate(mat44, n, 0, 0, 0, 0) * (3.*n*n - 2.*n);
		for (int i = 0; i < n; ++i)
			for (int j = i+1; j < n; ++j)
				M4 += evaluate(mat44, n, 0, 0, i, j) * (12.*n - 16.);
		for (int i = 0; i < n; ++i)
			for (int j = i+1; j < n; ++j)
				for (int k = j+1; k < n; ++k)
					for (int l = k+1; l < n; ++l)
						M4 += evaluate(mat44, n, i, j, k, l) * 24.;
		return std::real(M4) / std::pow(lat.n_sites(), 4.);
	}
};

struct wick_static_kek
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_static_kek(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf)
	{
		const numeric_t *ca_et_gf = et_gf.data();
		numeric_t kek = 0.;
		std::array<const std::vector<std::pair<int, int>>*, 3> single_kek_bonds =
			{&lat.bonds("single_kekule"), &lat.bonds("single_kekule_2"),
			&lat.bonds("single_kekule_3")};
		std::array<const std::vector<std::pair<int, int>>*, 3> kek_bonds =
			{&lat.bonds("kekule"), &lat.bonds("kekule_2"),
			&lat.bonds("kekule_3")};
		std::array<double, 3> factors = {-1., -1., 2.};
		const int N = kek_bonds.size(), M = single_kek_bonds[0]->size(), O = kek_bonds[0]->size(), ns = lat.n_sites();
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
			{
				auto& a = (*single_kek_bonds[i])[j];
				for (int m = 0; m < N; ++m)
					for (int n = 0; n < O; ++n)
					{
						auto& b = (*kek_bonds[m])[n];
						
						/*
						kek += factors[i] * factors[m]
								* (et_gf(a.second, a.first) * et_gf(b.first, b.second)
								+ lat.parity(a.first) * lat.parity(b.first) * et_gf(a.first, b.first) * et_gf(a.second, b.second));
						*/
								
						kek += factors[i] * factors[m]
							* (ca_et_gf[a.first*ns + a.second] * ca_et_gf[b.second*ns + b.first]
							+ lat.parity(a.first) * lat.parity(b.first) * ca_et_gf[b.first*ns + a.first] * ca_et_gf[b.second*ns + a.second]);
					}
			}
		return std::real(2.*kek) / std::pow(lat.n_bonds(), 2.);
	}
};
