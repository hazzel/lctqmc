#include <string>
#include <fstream>
#include <cmath>
#include <chrono>
#include <boost/algorithm/string.hpp>
#include "mc.h"
#include "honeycomb.h"
#include "move_functors.h"
#include "measure_functors.h"
#include "event_functors.h"
#ifdef PROFILER
	#include "gperftools/profiler.h"
#endif

mc::mc(const std::string& dir)
	: rng(Random()), qmc(rng), gf(rng, param, lat)
{
	//Read parameters
	pars.read_file(dir);
	sweep = 0;
	param.n_sweeps = pars.value_or_default<int>("SWEEPS", 0);
	n_static_cycles = pars.value_or_default<int>("static_cycles", 300);
	n_dyn_cycles = pars.value_or_default<int>("dyn_cycles", 300);
	n_warmup = pars.value_or_default<int>("warmup", 100000);
	param.n_prebin = pars.value_or_default<int>("prebin", 500);
	
	param.theta = pars.value_or_default<double>("theta", 40);
	param.block_size = pars.value_or_default<double>("block_size", 0.25);
	param.dyn_tau_max = pars.value_or_default<double>("dyn_tau_max", 10);
	param.dyn_delta_tau = pars.value_or_default<double>("dyn_delta_tau", 0.25);
	param.dyn_tau_steps = static_cast<int>((param.dyn_tau_max / param.dyn_delta_tau) + 0.5);
	
	param.L = pars.value_or_default<int>("L", 2);
	param.V = pars.value_or_default<double>("V", 1.0);
	param.t = pars.value_or_default<double>("t", 1.0);
	param.tprime = pars.value_or_default<double>("tprime", 0.0);
	param.stag_mu = pars.value_or_default<double>("stag_mu", 0.0);
	
	param.n_updates_per_block = pars.value_or_default<double>("updates_per_block", 1);
	param.static_measure_interval = pars.value_or_default<double>("static_measure_interval", 1);
	param.static_measure_cnt = 0;

	std::string static_obs_string = pars.value_or_default<std::string>("static_obs", "");
	boost::split(param.static_obs, static_obs_string, boost::is_any_of(","));
	std::string dyn_obs_string = pars.value_or_default<std::string>("obs", "");
	boost::split(param.dyn_obs, dyn_obs_string, boost::is_any_of(","));

	if (pars.defined("seed"))
		rng.NewRng(pars.value_of<int>("seed"));
	
	qmc.add_move(move_insert{rng, measure, param, lat, gf}, "insert", 1.0);
	qmc.add_move(move_remove{rng, measure, param, lat, gf}, "remove", 1.0);
	qmc.add_move(move_shift{rng, measure, param, lat, gf}, "shift", 1.0);
	qmc.add_measure(measure_M{measure, pars, param, gf}, "measurement");

	//Initialize lattice
	honeycomb hc(param.L, param.L);
	lat.generate_graph(hc);
	hc.generate_maps(lat);
	
	green_function::matrix_t K(lat.n_sites(), lat.n_sites());
	for (auto& a : lat.bonds("nearest neighbors"))
		K(a.first, a.second) = -param.t;
	for (auto& a : lat.bonds("t3_bonds"))
		K(a.first, a.second) = -param.tprime;
	gf.set_K_matrix(K);

	//Set up events
	qmc.add_event(event_build{rng, param, lat, gf}, "initial build");
	qmc.add_event(event_static_measurement{rng, measure, param, lat, gf}, "_static measure");
	qmc.add_event(event_dynamic_measurement{rng, measure, param, lat, gf}, "dynamic measure");

	#ifdef PROFILER
		ProfilerStart("/net/home/lxtsfs1/tpc/hesselmann/code/profiler/gperftools.prof");
	#endif
}

mc::~mc()
{
	#ifdef PROFILER
		ProfilerStop();
	#endif
}

void mc::random_write(odump& d)
{
	rng.RngHandle()->write(d);
}
void mc::seed_write(const std::string& fn)
{
	std::ofstream s;
	s.open(fn.c_str());
	s << rng.Seed() << std::endl;
	s.close();
}
void mc::random_read(idump& d)
{
	rng.NewRng();
	rng.RngHandle()->read(d);
}

void mc::init()
{
	qmc.init_moves();
	qmc.init_measurements();
	qmc.init_events();
	qmc.trigger_event("initial build");
}

void mc::write(const std::string& dir)
{
	odump d(dir+"dump");
	random_write(d);
	d.write(sweep);
	gf.serialize(d);
	d.close();
	seed_write(dir+"seed");
	std::ofstream f(dir+"bins");
	if (is_thermalized())
	{
		f << "Thermalization: Done." << std::endl
			<< "Sweeps: " << (sweep - n_warmup) << std::endl
			<< "Static bins: "
			<< std::endl
			<< "Dynamic bins: "
			<< std::endl;
	}
	else
	{
		f << "Thermalization: " << sweep << std::endl
			<< "Sweeps: 0" << std::endl
			<< "Static bins: 0" << std::endl
			<< "Dynamic bins: 0" << std::endl;
	}
	f.close();
}
bool mc::read(const std::string& dir)
{
	idump d(dir+"dump");
	if (!d)
	{
		std::cout << "read fail" << std::endl;
		return false;
	}
	else
	{
		random_read(d);
		d.read(sweep);
		gf.serialize(d);
		d.close();
		return true;
	}
}

void mc::write_output(const std::string& dir)
{
	std::ofstream f(dir);
	qmc.collect_results(f);
	f.close();
	/*
	const std::vector<std::pair<std::string, double>>& acc =
		qmc.acceptance_rates();
	for (auto a : acc)
		std::cout << a.first << " : " << a.second << std::endl;
	std::cout << "Average sign: " << qmc.average_sign() << std::endl;
	*/
}

bool mc::is_thermalized()
{
	return sweep >= n_warmup;
}

void mc::do_update()
{
	std::chrono::steady_clock::time_point t0, t1;
	for (int i = gf.tau() / param.block_size; i < param.theta / param.block_size; ++i)
	{
		t0 = std::chrono::steady_clock::now();
		gf.wrap((i + 0.5) * param.block_size);
		t1 = std::chrono::steady_clock::now();
		//std::cout << "Time of wrap: " << std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0).count() << std::endl;
		
		t0 = std::chrono::steady_clock::now();
		gf.rebuild();
		t1 = std::chrono::steady_clock::now();
		//std::cout << "Time of rebuild: " << std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0).count() << std::endl;
		if (is_thermalized())
		{
			qmc.do_measurement();
			qmc.trigger_event("_static measure");
			qmc.trigger_event("dynamic measure");
		}
		t0 = std::chrono::steady_clock::now();
		for (int n = 0; n < param.n_updates_per_block; ++n)
			qmc.do_update();
		t1 = std::chrono::steady_clock::now();
		//std::cout << "Time of updates: " << std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0).count() << std::endl;
		//std::cout << std::endl << "---" << std::endl;
	}
	for (int i = param.theta / param.block_size - 1; i >= 0; --i)
	{
		gf.wrap((i + 0.5) * param.block_size);
		gf.rebuild();
		if (is_thermalized())
		{
			qmc.do_measurement();
			qmc.trigger_event("_static measure");
			qmc.trigger_event("dynamic measure");
		}
		for (int n = 0; n < param.n_updates_per_block; ++n)
			qmc.do_update();
	}
	
	++sweep;
}

void mc::do_measurement()
{}

void mc::status()
{
	//if (sweep == n_warmup)
	//	std::cout << "Thermalization done." << std::endl;
	//if (is_thermalized() && sweep % (1000) == 0)
	//	std::cout << "sweep: " << sweep << std::endl;
}
