#include <string>
#include <fstream>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include "mc.h"
#include "honeycomb.h"
#include "move_functors.h"
#include "measure_functors.h"
#ifdef PROFILER
	#include "gperftools/profiler.h"
#endif

mc::mc(const std::string& dir)
	: rng(Random()), qmc(rng), gf(rng, param, lat)
{
	//Read parameters
	pars.read_file(dir);
	sweep = 0;
	int n_sweeps = pars.value_or_default<int>("SWEEPS", 0);
	n_static_cycles = pars.value_or_default<int>("static_cycles", 300);
	n_dyn_cycles = pars.value_or_default<int>("dyn_cycles", 300);
	n_warmup = pars.value_or_default<int>("warmup", 100000);
	n_prebin = pars.value_or_default<int>("prebin", 500);
	
	param.L = pars.value_or_default<int>("L", 2);
	param.theta = pars.value_or_default<double>("theta", 40);
	param.block_size = pars.value_or_default<double>("block_size", 1);
	param.V = pars.value_or_default<double>("V", 1.0);

	/*
	std::string static_obs_string = pars.value_or_default<std::string>("static_obs", "M2");
	boost::split(config.param.static_obs, static_obs_string, boost::is_any_of(","));
	std::string obs_string = pars.value_or_default<std::string>("obs", "M2");
	boost::split(config.param.obs, obs_string, boost::is_any_of(","));
	*/
	if (pars.defined("seed"))
		rng.NewRng(pars.value_of<int>("seed"));
	
	qmc.add_move(move_insert{rng, measure, param, lat, gf}, "insert", 1.0);
	qmc.add_move(move_remove{rng, measure, param, lat, gf}, "remove", 1.0);
	qmc.add_measure(measure_M{measure, pars}, "measurement");

	//Initialize lattice
	honeycomb hc(param.L, param.L);
	lat.generate_graph(hc);
	hc.generate_maps(lat);
	
	//Build H0
	green_function::matrix_t K(lat.n_sites(), lat.n_sites());
	for (auto& a : lat.bonds("nearest neighbors"))
		K(a.first, a.second) = -1.;
	//for (auto& a : lat.bonds("t3_bonds"))
	//	H0(a.first+as, a.second+as) = -param.tprime;
	gf.set_K_matrix(K);
	unsigned Nv = 40;
	green_function::vlist_t vlist;
	for (int i = 0; i < Nv; ++i)
	{
		double tau = rng() * param.theta;
		auto& b = lat.bonds("nearest neighbors")[rng() * 2 * lat.n_bonds()];
		vlist.push_back({tau, b.first, b.second});
	}
	std::sort(vlist.begin(), vlist.end(), vertex::less());
	gf.initialize(vlist);

	//Set up events

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
	qmc.init_events();
}

void mc::write(const std::string& dir)
{
	odump d(dir+"dump");
	random_write(d);
	d.write(sweep);
	d.write(static_bin_cnt);
	d.write(dyn_bin_cnt);
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
		d.read(static_bin_cnt);
		d.read(dyn_bin_cnt);
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
	for (int i = 0; i < param.theta / param.block_size - 1; ++i)
	{
		std::cout << "i = " << i << std::endl;
		gf.wrap(i * param.block_size);
		qmc.do_update();
		gf.rebuild();
		gf.measure();
		std::cout << "pert_order = " << gf.pert_order() << std::endl;
	}
	/*
	for (int i = param.block_size - 1; i > 0; --i)
	{
		gf.wrap(i * param.block_size);
		qmc.do_update();
		gf.rebuild();
		gf.measure();
		std::cout << "pert_order = " << gf.pert_order() << std::endl;
	}
	*/
	
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
