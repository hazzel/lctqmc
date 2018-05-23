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
