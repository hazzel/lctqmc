#pragma once
#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include "lattice.h"


struct pi_flux_square
{
	//typedef lattice::graph_t graph_t;
	typedef boost::adjacency_list<boost::setS, boost::vecS,
		boost::undirectedS> graph_t;

	int Lx;
	int Ly;
	std::vector<Eigen::Vector2d> real_space_map;
	// Base vectors of Bravais lattice
	Eigen::Vector2d a1;
	Eigen::Vector2d a2;
	// Base vectors of reciprocal lattice
	Eigen::Vector2d b1;
	Eigen::Vector2d b2;
	// Vector to second sublattice point
	Eigen::Vector2d delta;
	// Center of inversion symmetry
	Eigen::Vector2d center;
	double pi = 4. * std::atan(1.);

	pi_flux_square(int Lx_ = 6, int Ly_ = 6)
		: Lx(Lx_), Ly(Ly_),
			a1(1., 1.), a2(1., -1.), delta(1., 0.)
	{
		b1 = Eigen::Vector2d(pi, pi);
		b2 = Eigen::Vector2d(pi, -pi);
		center = Eigen::Vector2d(0.5, 0.);
	}
	
	int neighbor_site(int site, int type)
	{
		int i = site % (2 * Lx) / 2, n_vertices = 2 * Lx * Ly;
		if (type == 0)
			return (site + 1 + n_vertices) % n_vertices;
		else if (type == 1)
		{
			if (i == 0)
				return (site - 1 + n_vertices) % n_vertices;
			else
				return (site - 2*Lx - 1 + n_vertices) % n_vertices;
		}
		else if (type == 2)
			return (site - 2*Lx + 1 + n_vertices) % n_vertices;
		else if (type == 3)
		{
			if (i == 0)
				return (site + 2*Lx - 1 + n_vertices) % n_vertices;
			else
				return (site - 1 + n_vertices) % n_vertices;
		}
	}

	graph_t* graph(lattice& l)
	{
		int n_sites = 2 * Lx * Ly;
		graph_t* g = new graph_t(n_sites);
		add_edges(g);
		
		l.a1 = a1;
		l.a2 = a2;
		l.b1 = b1;
		l.b2 = b2;
		l.center = center;
		l.delta = delta;
		l.Lx = Lx;
		l.Ly = Ly;
	
		//Symmetry points
		std::map<std::string, Eigen::Vector2d> points;

		points["K"] = closest_k_point({pi, 0.5*pi});
		points["Kp"] = closest_k_point({pi, -0.5*pi});
		points["Gamma"] = closest_k_point({0., 0.});
		points["q"] = closest_k_point(b1 / Lx);
		points["q20"] = closest_k_point(2. * b2 / Ly);
		points["q11"] = closest_k_point(b1 / Lx + b2 / Ly);
		points["Kq"] = closest_k_point(b1 / Lx + Eigen::Vector2d{pi, 0.5*pi});
		l.add_symmetry_points(points);
		
		return g;
	}

	void add_edges(graph_t* g)
	{
		typedef std::pair<int, int> edge_t;
		for (int j = 0; j < Ly; ++j)
			for (int i = 0; i < Lx; ++i)
			{
				int n = j * 2 * Lx + i * 2;
				
				boost::add_edge(n, neighbor_site(n, 0), *g);
				boost::add_edge(neighbor_site(n, 0), n, *g);
				
				boost::add_edge(n, neighbor_site(n, 1), *g);
				boost::add_edge(neighbor_site(n, 1), n, *g);
				
				boost::add_edge(n, neighbor_site(n, 2), *g);
				boost::add_edge(neighbor_site(n, 2), n, *g);
				
				boost::add_edge(n, neighbor_site(n, 3), *g);
				boost::add_edge(neighbor_site(n, 3), n, *g);
				
				real_space_map.push_back(Eigen::Vector2d{i * a1 + j * a2});
				real_space_map.push_back(Eigen::Vector2d{i * a1 + j * a2 + delta});
			}
	}

	Eigen::Vector2d closest_k_point(const Eigen::Vector2d& K)
	{
		Eigen::Vector2d x = {0., 0.};
		double dist = (x - K).norm();
		for (int i = 0; i < Lx; ++i)
			for (int j = 0; j < Ly; ++j)
			{
				Eigen::Vector2d y = static_cast<double>(i) / static_cast<double>(Lx)
					* b1 + static_cast<double>(j) / static_cast<double>(Ly) * b2;
				double d = (y - K).norm();
				if (d < dist)
				{
					x = y;
					dist = d;
				}
			}
		return x;
	}

	void generate_maps(lattice& l)
	{
		//Site maps
		l.generate_neighbor_map("nearest neighbors", [&]
			(lattice::vertex_t i, lattice::vertex_t j) {
			return l.distance(i, j) == 1; });
		l.generate_bond_map("nearest neighbors", [&]
			(lattice::vertex_t i, lattice::vertex_t j)
			{ return l.distance(i, j) == 1; });
		l.generate_bond_map("single_d1_bonds", [&]
			(lattice::vertex_t i, lattice::vertex_t j)
			{ return l.distance(i, j) == 1 && i < j; });
		for (int d = 0; d < l.max_distance(); ++d)
		{
			l.generate_bond_map("d" + std::to_string(d) + "_bonds", [&]
				(lattice::vertex_t i, lattice::vertex_t j)
				{ return l.distance(i, j) == d; });
		}
	}
};
