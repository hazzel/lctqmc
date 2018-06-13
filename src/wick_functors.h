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
				M2 += lat.parity(i) * lat.parity(j)
					* std::real((1. - et_gf_t(i, i)) * (1. - et_gf_0(j, j))
					+ lat.parity(i) * lat.parity(j) * td_gf(i, j) * td_gf(i, j)
					- (et_gf_t(i, i) + et_gf_0(j, j))/2. + 1./4.);
				*/
				//M2 += td_gf(i, j) * td_gf(i, j);
				M2 += ca_td_gf[j * N + i] * ca_td_gf[j * N + i];
			}
		return std::real(M2) / std::pow(N, 2.);
	}
};

// ep(tau) = sum_{<ij>,<mn>} <c_i^dag(tau) c_j(tau) c_n^dag c_m>
struct wick_epsilon
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_epsilon(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
		numeric_t ep = 0.;
		auto& single_bonds = lat.bonds("single_d1_bonds");
		auto& bonds = lat.bonds("nearest neighbors");
		const int N = single_bonds.size(), M = bonds.size(), ns = lat.n_sites();
		for (int i = 0; i < N; ++i)
		{
			auto& a = single_bonds[i];
			for (int j = 0; j < M; ++j)
			{
				auto& b = bonds[j];
				
				//ep += et_gf_t(a.second, a.first) * et_gf_0(b.first, b.second)
				//	+ lat.parity(a.first) * lat.parity(b.first) * td_gf(a.first, b.first) * td_gf(a.second, b.second);
				
				ep += ca_et_gf_t[a.first*ns + a.second] * ca_et_gf_0[b.second*ns + b.first]
					+ lat.parity(a.first) * lat.parity(b.first) * ca_td_gf[b.first*ns + a.first] * ca_td_gf[b.second*ns + a.second];
				
			}
		}
		return std::real(2.*ep) / std::pow(N, 2.);
	}
};

struct wick_epsilon_as
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_epsilon_as(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
		numeric_t ep = 0.;
		std::vector<const std::vector<std::pair<int, int>>*> bonds =
			{&lat.bonds("nn_bond_1"), &lat.bonds("nn_bond_2"),
			&lat.bonds("nn_bond_3")};
		
		const int N = bonds.size(), M = bonds[0]->size(), ns = lat.n_sites();
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
			{
				auto& a = (*bonds[i])[j];
				for (int m = 0; m < N; ++m)
					for (int n = 0; n < M; ++n)
					{
						auto& b = (*bonds[m])[n];
						
						/*
						ep += 2.*(et_gf_t(a.second, a.first) * et_gf_0(b.first, b.second)
							+ lat.parity(a.first) * lat.parity(b.first) * td_gf(a.first, b.first) * td_gf(a.second, b.second));
							
						ep -= 2.*(et_gf_t(a.first, a.second) * et_gf_0(b.first, b.second)
							+ lat.parity(a.second) * lat.parity(b.first) * td_gf(a.second, b.first) * td_gf(a.first, b.second));
						*/
						ep += 2.*(ca_et_gf_t[a.first*ns+a.second] * ca_et_gf_0[b.second*ns+b.first]
							+ lat.parity(a.first) * lat.parity(b.first) * ca_td_gf[b.first*ns+a.first] * ca_td_gf[b.second*ns+a.second]);
							
						ep -= 2.*(ca_et_gf_t[a.second*ns+a.first] * ca_et_gf_0[b.second*ns+b.first]
							+ lat.parity(a.second) * lat.parity(b.first) * ca_td_gf[b.first*ns+a.second] * ca_td_gf[b.second*ns+a.first]);
						
						/*
						ep -= et_gf_t(a.second, a.first) * et_gf_0(b.second, b.first)
							+ lat.parity(a.first) * lat.parity(b.second) * td_gf(a.first, b.second) * td_gf(a.second, b.first);
							
						ep += et_gf_t(a.first, a.second) * et_gf_0(b.second, b.first)
							+ lat.parity(a.second) * lat.parity(b.second) * td_gf(a.second, b.second) * td_gf(a.first, b.first);
						*/
					}
			}
		return std::real(ep) / std::pow(N, 2.);
	}
};

// kekule(tau) = sum_{kekule} <c_i^dag(tau) c_j(tau) c_n^dag c_m>
struct wick_kekule_s
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_kekule_s(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
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
							* (et_gf_t(a.second, a.first) * et_gf_0(b.first, b.second)
							+ lat.parity(a.first) * lat.parity(b.first) * td_gf(a.first, b.first) * td_gf(a.second, b.second));
						*/
						
						kek += factors[i] * factors[m]
							* (ca_et_gf_t[a.first*ns + a.second] * ca_et_gf_0[b.second*ns + b.first]
							+ lat.parity(a.first) * lat.parity(b.first) * ca_td_gf[b.first*ns + a.first] * ca_td_gf[b.second*ns + a.second]);
					}
			}
		return std::real(2.*kek) / std::pow(lat.n_bonds(), 2.);
	}
};

// kekule(tau) = sum_{kekule} <c_i^dag(tau) c_j(tau) c_n^dag c_m>
struct wick_kekule_as
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_kekule_as(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
		numeric_t kek = 0.;
		std::array<const std::vector<std::pair<int, int>>*, 2> single_kek_bonds =
			{&lat.bonds("single_kekule"), &lat.bonds("single_kekule_2")};
		std::array<const std::vector<std::pair<int, int>>*, 2> kek_bonds =
			{&lat.bonds("kekule"), &lat.bonds("kekule_2")};
		std::array<double, 2> factors = {1., -1.};
		
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
							* (et_gf_t(a.second, a.first) * et_gf_0(b.first, b.second)
							+ lat.parity(a.first) * lat.parity(b.first) * td_gf(a.first, b.first) * td_gf(a.second, b.second));
						*/
						
						kek += factors[i] * factors[m]
							* (ca_et_gf_t[a.first*ns + a.second] * ca_et_gf_0[b.second*ns + b.first]
							+ lat.parity(a.first) * lat.parity(b.first) * ca_td_gf[b.first*ns + a.first] * ca_td_gf[b.second*ns + a.second]);
					}
			}
		return std::real(2.*kek) / std::pow(lat.n_bonds(), 2.);
	}
};

// kekule(tau) = sum_{kekule} <c_i^dag(tau) c_j(tau) c_n^dag c_m>
struct wick_kekule_K
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_kekule_K(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
		numeric_t kek = 0.;
		auto& K = lat.symmetry_point("K");
		//std::array<const std::vector<std::pair<int, int>>*, 3> nn_bonds =
		//	{&lat.bonds("nn_bond_1"), &lat.bonds("nn_bond_2"), &lat.bonds("nn_bond_3")};
		std::array<const std::vector<std::pair<int, int>>*, 1> nn_bonds = {&lat.bonds("nn_bond_2")};
		
		const int N = nn_bonds.size(), M = nn_bonds[0]->size(), ns = lat.n_sites();
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
			{
				auto& a = (*nn_bonds[i])[j];
				auto& r_a = lat.real_space_coord(a.second);
				//for (int m = 0; m < N; ++m)
				for (int n = 0; n < M; ++n)
				{
					auto& b = (*nn_bonds[i])[n];
					auto& r_b = lat.real_space_coord(b.second);
					
					double phase = std::cos(K.dot(r_a - r_b));
					
					kek += phase * (ca_et_gf_t[a.first*ns + a.second] * ca_et_gf_0[b.second*ns + b.first]
						+ lat.parity(a.first) * lat.parity(b.first) * ca_td_gf[b.first*ns + a.first] * ca_td_gf[b.second*ns + a.second]);
				}
			}
		return std::real(kek) / std::pow(lat.n_sites(), 2.);
	}
};

struct wick_gamma_mod
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_gamma_mod(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
		numeric_t gm = 0.;
		double pi = 4. * std::atan(1.);
		
		/*
		std::vector<const std::vector<std::pair<int, int>>*> bonds =
			{&lat.bonds("nn_bond_1"), &lat.bonds("nn_bond_2"),
			&lat.bonds("nn_bond_3")};
		std::vector<double> phases = {2.*std::sin(0. * pi), 2.*std::sin(2./3. * pi), 2.*std::sin(4./3. * pi)};
		*/
		std::vector<const std::vector<std::pair<int, int>>*> bonds =
			{&lat.bonds("nn_bond_2"), &lat.bonds("nn_bond_3")};
		std::vector<double> phases = {2.*std::sin(2./3. * pi), 2.*std::sin(4./3. * pi)};
		
		const int N = bonds.size(), M = bonds[0]->size(), ns = lat.n_sites();
		for (int i = 0; i < N; ++i)
			for (int j = 0; j < M; ++j)
			{
				auto& a = (*bonds[i])[j];
				for (int m = 0; m < N; ++m)
					for (int n = 0; n < M; ++n)
					{
						auto& b = (*bonds[m])[n];
						
						/*
						gm += 2.*phases[i] * phases[m]
							* (et_gf_t(a.second, a.first) * et_gf_0(b.first, b.second)
							+ lat.parity(a.first) * lat.parity(b.first) * td_gf(a.first, b.first) * td_gf(a.second, b.second));
							
						gm -= 2.*phases[i] * phases[m]
							* (et_gf_t(a.first, a.second) * et_gf_0(b.first, b.second)
							+ lat.parity(a.second) * lat.parity(b.first) * td_gf(a.second, b.first) * td_gf(a.first, b.second));
						*/
						gm += 2.*phases[i] * phases[m]
							* (ca_et_gf_t[a.first*ns+a.second] * ca_et_gf_0[b.second*ns+b.first]
							+ lat.parity(a.first) * lat.parity(b.first) * ca_td_gf[b.first*ns+a.first] * ca_td_gf[b.second*ns+a.second]);
							
						gm -= 2.*phases[i] * phases[m]
							* (ca_et_gf_t[a.second*ns+a.first] * ca_et_gf_0[b.second*ns+b.first]
							+ lat.parity(a.second) * lat.parity(b.first) * ca_td_gf[b.first*ns+a.second] * ca_td_gf[b.second*ns+a.first]);
						
						/*
						gm -= phases[i] * phases[m]
							* (et_gf_t(a.second, a.first) * et_gf_0(b.second, b.first)
							+ lat.parity(a.first) * lat.parity(b.second) * td_gf(a.first, b.second) * td_gf(a.second, b.first));
							
						gm += phases[i] * phases[m]
							* (et_gf_t(a.first, a.second) * et_gf_0(b.second, b.first)
							+ lat.parity(a.second) * lat.parity(b.second) * td_gf(a.second, b.second) * td_gf(a.first, b.first));
						*/
					}
			}
		return std::real(gm) / std::pow(lat.n_bonds(), 2.);
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

struct wick_tp
{
	Random& rng;
	parameters& param;
	lattice& lat;

	wick_tp(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{}
	
	double get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_td_gf = td_gf.data();
		auto& K = lat.symmetry_point("K");
		auto& Kp = lat.symmetry_point("Kp");
		const int N = lat.n_sites();
		numeric_t tp = 0.;
		int i = 2 * static_cast<int>(rng() * N / 2);
		auto& r_i = lat.real_space_coord(i);
		for (int m = 0; m < N; m+=2)
		{
			auto& r_m = lat.real_space_coord(m);
			int j = 2 * static_cast<int>(rng() * N / 2);
			auto& r_j = lat.real_space_coord(j);
			for (int n = 0; n < N; n+=2)
			{
				auto& r_n = lat.real_space_coord(n);
				//double kdot = K.dot(r_i - r_j - r_m + r_n);
				double kdot = K.dot(r_i - r_m) + Kp.dot(r_j - r_n);
				tp += std::real(std::cos(kdot) * (ca_td_gf[i*N+m] * ca_td_gf[(j+1)*N+(n+1)] - ca_td_gf[i*N+(n+1)] * ca_td_gf[(j+1)*N+m]));
			}
		}
		return tp / static_cast<double>(N*N/4);
	}
};
