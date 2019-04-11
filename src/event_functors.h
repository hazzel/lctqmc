#pragma once
#include <vector>
#include <chrono>
#include "measurements.h"
#include "parameters.h"
#include "lattice.h"
#include "green_function.h"
#include "wick_static_base.h"
#include "vector_wick_static_base.h"
#include "wick_base.h"
#include "vector_wick_base.h"
#include "wick_static_functors.h"
#include "vector_wick_static_functors.h"
#include "wick_functors.h"
#include "vector_wick_functors.h"
#include "event_set_trial_wf.h"

struct event_build
{
	Random& rng;
	parameters& param;
	lattice& lat;
	green_function& gf;

	void trigger()
	{
		int Nv = 0.15 * std::abs(param.theta * lat.n_sites() * param.V);
		//int Nv = 2*param.theta;
		green_function::vlist_t vlist;
		for (int i = 0; i < Nv; ++i)
		{
			double tau = rng() * param.theta;
			auto& b = lat.bonds("nearest neighbors")[rng() * 2 * lat.n_bonds()];
			vlist.push_back({tau, b.first, b.second});
		}
		std::sort(vlist.begin(), vlist.end(), vertex::less());
		gf.initialize(0, vlist);
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
			else if (param.static_obs[i] == "S_cdw_q20")
				add_wick(wick_static_S_cdw_q20{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "S_cdw_q11")
				add_wick(wick_static_S_cdw_q11{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "chi_cdw")
				add_vector_wick(wick_static_chi_cdw{rng, param, lat}, param.static_obs[i], wick_static_chi_cdw::nq+1);
			else if (param.static_obs[i] == "M4")
				add_wick(wick_static_M4{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "epsilon")
				add_wick(wick_static_epsilon{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "Hv")
				add_wick(wick_static_Hv{rng, param, lat, gf}, param.static_obs[i]);
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
	void add_vector_wick(T&& functor, const std::string& name, int n_values)
	{
		vec_obs.push_back(vector_wick_static_base<matrix_t>(std::forward<T>(functor), n_values));
		vec_names.push_back(name);
	}

	void trigger()
	{
		if (obs.size() == 0 and vec_obs.size() == 0)
			return;
		++param.static_measure_cnt;
		if (param.static_measure_cnt >= param.static_measure_interval)
		if ((param.projective and std::abs(gf.tau() - param.theta/2.) < param.measure_window/2.)
			or (!param.projective))
		{
			//std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
			gf.measure_static_observables(measure, names, obs, vec_names, vec_obs);
			//std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
			//std::cout << "Time of static measurement: " << std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0).count() << std::endl;
			param.static_measure_cnt = 0;
		}
		if (param.ep_tau_steps > 0)
		{
			std::vector<double> hv_tau = gf.measure_Hv_tau();
			measure.add("dyn_Hv_tau", hv_tau);
		}
	}
	
	void init()
	{
		int static_prebin = param.n_prebin * (param.measure_window / param.block_size) / param.static_measure_interval;
		for (int i = 0; i < obs.size(); ++i)
			measure.add_observable(names[i], static_prebin);
		for (int i = 0; i < vec_obs.size(); ++i)
			measure.add_vectorobservable(vec_names[i], vec_obs[i].n_values, static_prebin);
		if (param.ep_tau_steps > 0)
			measure.add_vectorobservable("dyn_Hv_tau", param.ep_tau_steps, static_prebin);
	}
};

struct event_dynamic_measurement
{
	typedef green_function::matrix_t matrix_t;

	Random& rng;
	measurements& measure;
	parameters& param;
	lattice& lat;
	green_function& gf;
	std::vector<std::vector<double>> dyn_tau;
	std::vector<wick_base<matrix_t>> obs;
	std::vector<vector_wick_base<matrix_t>> vec_obs;
	std::vector<std::string> names;
	std::vector<std::string> vec_names;

	event_dynamic_measurement(Random& rng_, measurements& measure_,
		parameters& param_, lattice& lat_, green_function& gf_)
		: rng(rng_), measure(measure_), param(param_), lat(lat_), gf(gf_)
	{
		obs.reserve(param.dyn_obs.size());
		vec_obs.reserve(param.dyn_obs.size());
		for (int i = 0; i < param.dyn_obs.size(); ++i)
		{
			if (param.dyn_obs[i] == "M2")
				add_wick(wick_M2{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "chi_cdw")
				add_vector_wick(wick_chi_cdw{rng, param, lat}, param.dyn_obs[i], wick_chi_cdw::nq+1);
			else if (param.dyn_obs[i] == "epsilon")
				add_wick(wick_epsilon{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "epsilon_as")
				add_wick(wick_epsilon_as{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "kekule_s")
				add_wick(wick_kekule_s{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "kekule_as")
				add_wick(wick_kekule_as{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "kekule_K")
				add_wick(wick_kekule_K{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "gamma_mod")
				add_wick(wick_gamma_mod{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "2d_rep")
				add_wick(wick_2d_rep{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "chern")
				add_wick(wick_chern{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "sp")
				add_wick(wick_sp{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "sp_q")
				add_vector_wick(wick_sp_q{rng, param, lat}, param.dyn_obs[i], (wick_sp_q::nq+1)*(wick_sp_q::nq+1));
			else if (param.dyn_obs[i] == "sp_k")
				add_vector_wick(wick_sp_k{rng, param, lat}, param.dyn_obs[i], lat.Lx * lat.Ly);
			else if (param.dyn_obs[i] == "tp")
				add_wick(wick_tp{rng, param, lat}, param.dyn_obs[i]);
			else if (param.dyn_obs[i] == "tp_mat")
				add_vector_wick(wick_tp_matrix{rng, param, lat}, param.dyn_obs[i], 4);
		}
		for (int i = 0; i < obs.size(); ++i)
			dyn_tau.push_back(std::vector<double>(param.dyn_tau_steps + 1, 0.));
		for (int i = 0; i < vec_obs.size(); ++i)
			for (int j = 0; j < vec_obs[i].n_values; ++j)
				dyn_tau.push_back(std::vector<double>(param.dyn_tau_steps + 1, 0.));
	}

	template<typename T>
	void add_wick(T&& functor, const std::string& name)
	{
		obs.push_back(wick_base<matrix_t>(std::forward<T>(functor)));
		names.push_back(name);
	}
	
	template<typename T>
	void add_vector_wick(T&& functor, const std::string& name, int n_values)
	{
		vec_obs.push_back(vector_wick_base<matrix_t>(std::forward<T>(functor), n_values));
		vec_names.push_back(name);
	}

	void trigger()
	{
		if (obs.size() == 0 and vec_obs.size() == 0)
			return;
		
		if (param.projective and (std::abs(gf.tau() - (param.theta/2. + param.dyn_tau_max/2 + param.block_size/2.)) < 1E-6
			or std::abs(gf.tau() - (param.theta/2. - param.dyn_tau_max/2 + param.block_size/2.)) < 1E-6))
		//if (std::abs(gf.tau() - (param.theta/2. + param.dyn_tau_max/2 + param.block_size/2.)) <= 1E-6)
		//if (std::abs(gf.tau() - (param.theta/2. - param.dyn_tau_max/2 + param.block_size/2.)) <= 1E-6)
		{
			//std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
			for (int i = 0; i < dyn_tau.size(); ++i)
				std::fill(dyn_tau[i].begin(), dyn_tau[i].end(), 0.);
			gf.measure_dynamical_observables(dyn_tau, names, obs, vec_names, vec_obs);
			
			for (int i = 0; i < obs.size(); ++i)
				measure.add("dyn_"+names[i]+"_tau", dyn_tau[i]);
			int cnt = 0;
			for (int i = 0; i < vec_obs.size(); ++i)
				for (int j = 0; j < vec_obs[i].n_values; ++j)
				{
					measure.add("dyn_"+vec_names[i]+"_"+std::to_string(j)+"_tau", dyn_tau[obs.size()+cnt]);
					++cnt;
				}
			//std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
			//std::cout << "Time of dynamic measurement: " << std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0).count() << std::endl;
		}
	}
	
	void init()
	{
		for (int i = 0; i < obs.size(); ++i)
			measure.add_vectorobservable("dyn_"+names[i]+"_tau", param.dyn_tau_steps+1, param.n_prebin);
		for (int i = 0; i < vec_obs.size(); ++i)
			for (int j = 0; j < vec_obs[i].n_values; ++j)
				measure.add_vectorobservable("dyn_"+vec_names[i]+"_"+std::to_string(j)+"_tau", param.dyn_tau_steps+1, param.n_prebin);
	}
};
