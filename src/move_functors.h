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
		measure.add_observable("insertion", param.theta / param.block_size * param.n_sweeps * param.n_updates_per_block / param.n_prebin);
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
		measure.add_observable("removal", param.theta / param.block_size * param.n_sweeps * param.n_updates_per_block / param.n_prebin);
	}
};

struct move_shift
{
	Random& rng;
	measurements& measure;
	parameters& param;
	lattice& lat;
	green_function& gf;
	
	green_function::vlist_t::iterator vpos;
	int si_prime;
	int sj_prime;
	
	move_shift(Random& rng_, measurements& measure_, parameters& param_, lattice& lat_, green_function& gf_)
		: rng(rng_), measure(measure_), param(param_), lat(lat_), gf(gf_)
	{}

	double attempt()
	{
		int num_vertices = gf.select_random_vertex(vpos);
		if (num_vertices > 0)
		{
			si_prime = vpos->si;
			sj_prime = vpos->sj;
			auto& neighbors = lat.neighbors(vpos->si, "nearest neighbors");
			while (sj_prime == vpos->sj)
				sj_prime = neighbors[rng() * neighbors.size()];
		}
		return gf.shift_vertex(vpos, si_prime, sj_prime, true);
	}

	double accept()
	{
		gf.shift_vertex(vpos, si_prime, sj_prime, false);
		measure.add("shift", 1.0);
		return 1.0;
	}

	void reject()
	{
		measure.add("shift", 0.0);
	}
	
	void init()
	{
		measure.add_observable("shift", param.theta / param.block_size * param.n_sweeps * param.n_updates_per_block / param.n_prebin);
	}
};
