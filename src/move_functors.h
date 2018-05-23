#pragma once
#include <vector>
#include <algorithm>
#include "measurements.h"
#include "lattice.h"
#include "Random.h"
#include "parameters.h"
#include "green_function.h"

struct move_insert
{
	Random& rng;
	measurements& measure;
	parameters& param;
	lattice& lat;
	green_function& gf;
	
	vertex v;
	
	move_insert(Random& rng_, measurements& measure_, parameters& param_, lattice& lat_, green_function& gf_)
		: rng(rng_), measure(measure_), param(param_), lat(lat_), gf(gf_)
	{}

	double attempt()
	{
		v = gf.generate_random_vertex();
		return -0.25 * param.block_size * param.V * lat.n_bonds() * gf.add_vertex(v, true);
	}

	double accept()
	{
		gf.add_vertex(v, false);
		measure.add("insertion", 1.0);
		return 1.0;
	}

	void reject()
	{
		measure.add("insertion", 0.0);
	}
	
	void init()
	{
		measure.add_observable("insertion", 1000);
	}
};

struct move_remove
{
	Random& rng;
	measurements& measure;
	parameters& param;
	lattice& lat;
	green_function& gf;
	
	green_function::vlist_t::iterator vpos;
	
	move_remove(Random& rng_, measurements& measure_, parameters& param_, lattice& lat_, green_function& gf_)
		: rng(rng_), measure(measure_), param(param_), lat(lat_), gf(gf_)
	{}

	double attempt()
	{
		int num_vertices = gf.select_random_vertex(vpos);
		return -4 * num_vertices / (param.block_size * param.V * lat.n_bonds()) * gf.remove_vertex(vpos, true);
	}

	double accept()
	{
		gf.remove_vertex(vpos, false);
		measure.add("removal", 1.0);
		return 1.0;
	}

	void reject()
	{
		measure.add("removal", 0.0);
	}
	
	void init()
	{
		measure.add_observable("removal", 1000);
	}
};
