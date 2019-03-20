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
	static constexpr int nq = 4;

	wick_chi_cdw(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(nq+1);
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_et_gf_0 = et_gf_0.data(), *ca_et_gf_t = et_gf_t.data(), *ca_td_gf = td_gf.data();
		const int N = lat.n_sites();
		std::fill(values.begin(), values.end(), 0.);
		for (int q = 0; q <= nq; ++q)
		{
			auto gamma_q = q * lat.b1 / lat.Lx;
			for (int i = 0; i < N; i+=2)
			{
				auto& r_i = lat.real_space_coord(i);
				for (int j = 0; j < N; j+=2)
				{
					auto& r_j = lat.real_space_coord(j);
					double qr = gamma_q.dot(r_i - r_j);
					numeric_t aa = ca_td_gf[j * N + i], ab = ca_td_gf[(j+1) * N + i], ba = ca_td_gf[j * N + i+1], bb = ca_td_gf[(j+1) * N + (i+1)];
					values[q] += std::cos(qr) * std::real(aa*aa + ab*ab + ba*ba + bb*bb);
				}
			}
			values[q] /= N * N;
		}
		return values;
	}
};

struct wick_sp_k
{
	Random& rng;
	parameters& param;
	lattice& lat;
	std::vector<double> values;
	std::vector<double> fourier_coeff;

	wick_sp_k(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize(lat.Lx * lat.Ly);
		
		double a1b1 = lat.a1.dot(lat.b1) / lat.Lx;
		double a1b2 = lat.a1.dot(lat.b2) / lat.Ly;
		double a2b1 = lat.a2.dot(lat.b1) / lat.Lx;
		double a2b2 = lat.a2.dot(lat.b2) / lat.Ly;
		for (int kx = 0; kx < lat.Lx; ++kx)
			for (int ky = 0; ky < lat.Ly; ++ky)
				for (int ix = -lat.Lx; ix <= lat.Lx; ++ix)
					for (int iy = -lat.Ly; iy <= lat.Ly; ++iy)
						fourier_coeff.emplace_back(std::cos(ix * (kx * a1b1 + ky * a1b2) + iy * (kx * a2b1 + ky * a2b2)));
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_td_gf = td_gf.data();
		const int N = lat.n_sites();
		std::fill(values.begin(), values.end(), 0.);
		
		for (int kx = 0; kx < lat.Lx; ++kx)
			for (int ky = 0; ky < lat.Ly; ++ky)
			{
				for (int ix = 0; ix < lat.Lx; ++ix)
					for (int iy = 0; iy < lat.Ly; ++iy)
					{
						int i = ix * lat.Ly + iy;
						for (int jx = 0; jx < lat.Lx; ++jx)
							for (int jy = 0; jy < lat.Ly; ++jy)
							{
								int j = jx * lat.Ly + jy;
								values[kx*lat.Ly+ky] += get_fourier_coeff(kx, ky, ix, jx, iy, jy) * std::real(ca_td_gf[(2*j+0)*N+(2*i+0)] + ca_td_gf[(2*j+1)*N+(2*i+1)]);
							}
					}
				values[kx*lat.Ly+ky] /= N;
			}
		return values;
	}
	
	double get_fourier_coeff(int kx, int ky, int ix, int jx, int iy, int jy)
	{
		int x = ix - jx + lat.Lx;
		int y = iy - jy + lat.Ly;
		return fourier_coeff[kx*lat.Ly*(2*lat.Lx+1)*(2*lat.Ly+1) + ky*(2*lat.Lx+1)*(2*lat.Ly+1) + x*(2*lat.Ly+1) + y];
	}
};

struct wick_sp_q
{
	Random& rng;
	parameters& param;
	lattice& lat;
	std::vector<double> values;
	static constexpr int nq = 2;

	wick_sp_q(Random& rng_, parameters& param_, lattice& lat_)
		: rng(rng_), param(param_), lat(lat_)
	{
		values.resize((nq+1)*(nq+1));
	}
	
	std::vector<double>& get_obs(const matrix_t& et_gf_0, const matrix_t& et_gf_t,
		const matrix_t& td_gf)
	{
		const numeric_t *ca_td_gf = td_gf.data();
		const int N = lat.n_sites();
		std::fill(values.begin(), values.end(), 0.);
		for (int p = 0; p <= nq; ++p)
		{
			for (int q = 0; q <= nq; ++q)
			{
				auto K = lat.symmetry_point("K") + q * lat.b1 / lat.Lx + p * lat.b2 / lat.Ly;
				for (int i = 0; i < N; i+=2)
				{
					auto& r_i = lat.real_space_coord(i);
					for (int j = 0; j < N; j+=2)
					{
						auto& r_j = lat.real_space_coord(j);
						double kdot = K.dot(r_i - r_j);
						//values[p*(nq+1)+q] += std::cos(kdot) * (ca_td_gf[(j+0)*N+(i+0)] + ca_td_gf[(j+1)*N+(i+1)] + ca_td_gf[(j+0)*N+(i+1)] + ca_td_gf[(j+1)*N+(i+0)]);
						values[p*(nq+1)+q] += std::cos(kdot) * std::real(ca_td_gf[(j+0)*N+(i+0)] + ca_td_gf[(j+1)*N+(i+1)]);
					}
				}
				values[p*(nq+1)+q] /= N;
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
