#pragma once
#include <map>
#include <vector>
#include <chrono>
#include "measurements.h"
#include "parameters.h"
#include "lattice.h"
#include "green_function.h"
#include "wick_static_base.h"
#include "vector_wick_static_base.h"
#include "wick_static_functors.h"

struct event_build
{
	Random& rng;
	parameters& param;
	lattice& lat;
	green_function& gf;

	void trigger()
	{
		green_function::matrix_t K(lat.n_sites(), lat.n_sites());
		for (auto& a : lat.bonds("nearest neighbors"))
			K(a.first, a.second) = -param.t;
		for (auto& a : lat.bonds("t3_bonds"))
			K(a.first, a.second) = -param.tprime;
		gf.set_K_matrix(K);
		unsigned Nv = 0.15 * (param.theta * lat.n_sites() * param.V);
		green_function::vlist_t vlist;
		for (int i = 0; i < Nv; ++i)
		{
			double tau = rng() * param.theta;
			auto& b = lat.bonds("nearest neighbors")[rng() * 2 * lat.n_bonds()];
			vlist.push_back({tau, b.first, b.second});
		}
		std::sort(vlist.begin(), vlist.end(), vertex::less());
		gf.initialize(vlist);
	}
	
	void init() {}
};

struct event_static_measurement
{
	typedef green_function::matrix_t matrix_t;

	Random& rng;
	measurements& measure;
	parameters& param;
	lattice& lat;
	green_function& gf;
	std::vector<wick_static_base<matrix_t>> obs;
	std::vector<vector_wick_static_base<matrix_t>> vec_obs;
	std::vector<std::string> names;
	std::vector<std::string> vec_names;

	event_static_measurement(Random& rng_, measurements& measure_, parameters& param_, lattice& lat_,
		green_function& gf_)
		: rng(rng_), measure(measure_), param(param_), lat(lat_), gf(gf_)
	{
		obs.reserve(param.static_obs.size());
		vec_obs.reserve(param.static_obs.size());
		for (int i = 0; i < param.static_obs.size(); ++i)
		{
			if (param.static_obs[i] == "energy")
				add_wick(wick_static_energy{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "M2")
				add_wick(wick_static_M2{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "S_cdw_q")
				add_wick(wick_static_S_cdw_q{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "M4")
				add_wick(wick_static_M4{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "epsilon")
				add_wick(wick_static_epsilon{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "kekule")
				add_wick(wick_static_kek{rng, param, lat}, param.static_obs[i]);
		}
	}

	template<typename T>
	void add_wick(T&& functor, const std::string& name)
	{
		obs.push_back(wick_static_base<matrix_t>(std::forward<T>(functor)));
		names.push_back(name);
	}
	
	template<typename T>
	void add_vector_wick(T&& functor, const std::string& name)
	{
		vec_obs.push_back(vector_wick_static_base<matrix_t>(std::forward<T>(functor)));
		vec_names.push_back(name);
	}

	void trigger()
	{
		if (std::abs(gf.tau() - param.theta/2.) < param.theta/8.)
			gf.measure_static_observable(measure, names, obs, vec_names, vec_obs);
	}
	
	void init()
	{
		for (int i = 0; i < obs.size(); ++i)
			measure.add_observable(names[i], param.n_prebin);
	}
};
