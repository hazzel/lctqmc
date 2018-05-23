#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <cassert>
#include <vector>
#include "Random.h"
#include "lattice.h"
#include "parameters.h"
#include "measurements.h"
#include "wick_static_base.h"
#include "vector_wick_static_base.h"

template<typename T, typename Pred>
typename std::vector<T>::iterator insert_sorted(std::vector<T>& vec, const T& item, Pred pred)
{
	return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, pred), item);
}

struct vertex
{
	double tau;
	int si;
	int sj;
	
	struct less
	{
		constexpr bool operator()(const vertex& lhs, const vertex& rhs) const
		{
			return lhs.tau < rhs.tau;
		}
		constexpr bool operator()(double lhs, const vertex& rhs) const
		{
			return lhs < rhs.tau;
		}
		constexpr bool operator()(const vertex& lhs, double rhs) const
		{
			return lhs.tau < rhs;
		}
	};

	struct less_equal
	{
		constexpr bool operator()(const vertex& lhs, const vertex& rhs) const
		{
			return lhs.tau <= rhs.tau;
		}
		constexpr bool operator()(double lhs, const vertex& rhs) const
		{
			return lhs <= rhs.tau;
		}
		constexpr bool operator()(const vertex& lhs, double rhs) const
		{
			return lhs.tau <= rhs;
		}
	};
};  

class green_function
{
	public:
		using numeric_t = double;
		using vector_t = Eigen::VectorXd;
		using matrix_t = Eigen::Matrix<numeric_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
		using vlist_t = std::vector<vertex>;

		green_function(Random& rng_, parameters& param_, lattice& lat_)
			: rng(rng_), param(param_), lat(lat_)
		{}
		
		void set_K_matrix(const matrix_t& K_)
		{
			K = K_;
			Eigen::SelfAdjointEigenSolver<matrix_t> solver(K);
			wK = solver.eigenvalues(); 
			uK = solver.eigenvectors();
			uKdag = solver.eigenvectors().adjoint();
			
			uKdagP = uKdag * solver.eigenvectors().leftCols(lat.n_sites()/2);
			Eigen::JacobiSVD<matrix_t> svd(uKdagP, Eigen::ComputeThinU);
			uKdagP = svd.matrixU();
		}
		
		unsigned int pert_order()
		{
			return vlist.size();
		}
		
		double tau()
		{
			return tpos;
		}
		
		// it can do  B(tau_m)... B(tau_n) * A  when sign = -1
		// or        A* B(tau_n)^{-1} ... B(tau_m)^{-1} when sign = 1 
		// Btau(tau_n) does not contain vertex contribution  
		void prop_from_left(const int sign, const double tau_m, const double tau_n, matrix_t& A) const
		{
			std::string side = (sign == -1) ? "L" : "R";
		
			//since we always point to the right of a vertex, we push downwards the lower boundary 
			//this is different with block convention (,]  
			vlist_t::const_iterator lower = std::lower_bound(vlist.begin(), vlist.end(), tau_n, vertex::less_equal());	//equal is exclude
			vlist_t::const_iterator upper = std::upper_bound(vlist.begin(), vlist.end(), tau_m, vertex::less());

			//upper > tau_m > lower > tau_n
			if (lower == upper)	// there is no vertex in between tau1 and tau2 
				K_prop(sign, tau_m - tau_n, side, A);
			else
			{
				K_prop(sign, lower->tau - tau_n, side, A);
				for (vlist_t::const_iterator it = lower; it != upper; ++it)
				{
					V_prop(it->si, it->sj, side, A);
					if (upper - it > 1)
						K_prop(sign, std::next(it)->tau - it->tau, side, A);
					else
						K_prop(sign, tau_m - it->tau, side, A);
				}
			}
		}

		// it can do A* B(tau_m)... B(tau_n)  when sign = -1
		// or        B(tau_n)^{-1} ... B(tau_m)^{-1}* A when sign = 1 
		void prop_from_right(const int sign, const double tau_m, const double tau_n, matrix_t& A) const
		{
			std::string side = (sign == 1) ? "L" : "R";
	
			vlist_t::const_iterator lower = std::lower_bound(vlist.begin(), vlist.end(), tau_n, vertex::less_equal());	//equal is exclude
			vlist_t::const_iterator upper = std::upper_bound(vlist.begin(), vlist.end(), tau_m, vertex::less());
			
			//upper > tau_m > lower > tau_n
			if (lower == upper)	// there is no vertex in between tau1 and tau2 
				K_prop(sign, tau_m - tau_n, side, A); 
			else
			{
				K_prop(sign, tau_m - (--upper)->tau, side, A);

				for (vlist_t::const_iterator it = upper; it != lower; --it)
				{
					V_prop(it->si, it->sj, side, A);
					K_prop(sign, it->tau - std::prev(it)->tau, side, A); 
				}
				//the last step by hand (it1 = lower, it2 point to tau2) 
				V_prop(lower->si, lower->sj, side, A);
				K_prop(sign, lower->tau - tau_n , side, A); 
			}
		}

		void K_prop(const int sign, double tau, const std::string& side, matrix_t& A) const
		{
			if (std::abs(tau) < 1E-15)
					return;

			if (side == "L")
			{	//  exp(tau * w) * A 
				for (int l=0; l<lat.n_sites(); ++l) 
					A.row(l) *= std::exp(sign * tau * wK[l]);
			}
			else if (side == "R")
			{	// A * exp(tau * w)
				for (int l=0; l<lat.n_sites(); ++l) 
					A.col(l) *= std::exp(sign * tau * wK[l]);
			}
		}

		void V_prop(int si, int sj, const std::string& side,  matrix_t& A) const // A*U^{dagger} V U or U^{dagger} V U * A 
		{
			if (side == "L")
				A.noalias() -= 2.* uKdag.col(si) * (uK.row(si)* A) + 2.* uKdag.col(sj) * (uK.row(sj)* A); 
			else if (side == "R")
				A.noalias() -= 2.* (A*uKdag.col(si))* uK.row(si) + 2.* (A*uKdag.col(sj)) * uK.row(sj);
		}

		void initialize(const std::vector<vertex>& vlist_)
		{
			storage.resize(param.theta / param.block_size + 1);

			vlist = vlist_;
			//tpos = vlist.begin()->tau;
			tpos = 0.;
			
			int block = tpos / param.block_size;
			storage[0] = uKdagP;
			
			for (int i=0; i < block; ++i)
			{
				matrix_t UR = storage[i];
				prop_from_left(-1, (i+1)*param.block_size, i*param.block_size, UR);

				Eigen::JacobiSVD<matrix_t> svd(UR, Eigen::ComputeThinU); 
				storage[i+1] = svd.matrixU();
			}

			storage.back() = uKdagP.adjoint();
			for (int i = storage.size()-2; i>block; --i)
			{
				matrix_t VL = storage[i+1];
				prop_from_right(-1, (i+1)*param.block_size, i*param.block_size, VL);

				Eigen::JacobiSVD<matrix_t> svd(VL, Eigen::ComputeThinV); 
				storage[i] = svd.matrixV().adjoint();
			}

			g_tau = g_stable();
		}
		
		void rebuild()
		{
			matrix_t g_stab = g_stable();

			double err = (g_tau - g_stab).norm();

			if (err > 1E-8)
				std::cout << "Error: " << err << std::endl;

			g_tau = g_stab;
		}
		
		matrix_t g_stable()
		{
			int block = tpos / param.block_size;

			matrix_t UR = storage[block];
			prop_from_left(-1, tpos, block*param.block_size, UR);
		
			matrix_t VL = storage[block+1];
			prop_from_right(-1, (block+1)*param.block_size, tpos, VL);
			
			matrix_t res = -UR * ((VL*UR).inverse() * VL);
			for (int l = 0; l < lat.n_sites(); ++l)
				res(l, l) += 1.0;

			return res; 
		}
		
		matrix_t g_exact()
		{
			matrix_t UR = uKdagP;
			prop_from_left(-1, tpos, 0, UR);
			matrix_t VL = uKdagP.adjoint();
			prop_from_right(-1, param.theta, tpos, VL);
			matrix_t id = matrix_t::Identity(lat.n_sites(), lat.n_sites());
			return id - UR * (VL * UR).inverse() * VL;
		}
		
		void wrap(double tau)
		{
			int old_block = tpos / param.block_size;
			int new_block = tau / param.block_size;
			
			//wrap Green's function 
			if (tau >= tpos)
			{
				// B G B^{-1}
				prop_from_left(-1, tau, tpos, g_tau);	// B(tau1) ... B(tau2) *U_  
				prop_from_left(1, tau, tpos, g_tau);	// V_ * B^{-1}(tau2) ... B^{-1}(tau1)
				tpos = tau;
			}
			else
			{
				// B^{-1} G B 
				prop_from_right(1, tpos, tau, g_tau);	//  B^{-1}(tau2) ... B^{-1}(tau1) * U_
				prop_from_right(-1, tpos, tau, g_tau);	//  V_ * B(tau1) ... B(tau2)
				tpos = tau;
			}
			
			//when we wrap to a new block we need to update storage 
			if (new_block > old_block)
			{// move to a larger block on the left
				matrix_t UR = storage[old_block];
				prop_from_left(-1, new_block*param.block_size, old_block*param.block_size, UR);
				Eigen::JacobiSVD<matrix_t> svd(UR, Eigen::ComputeThinU); 
				storage[new_block] = svd.matrixU();
			}
			else if (new_block < old_block)
			{// move to smaller block
				matrix_t VL = storage[old_block+1];
				prop_from_right(-1, (old_block+1)*param.block_size, old_block*param.block_size, VL);
				Eigen::JacobiSVD<matrix_t> svd(VL, Eigen::ComputeThinV);
				storage[new_block+1] = svd.matrixV().adjoint();
			}
		}
		
		double gij(const int si, const int sj) const
		{	// current g in the site basis 
			// (U gtau U^{dagger} )_ij 
			return  (uK.row(si) * g_tau) * uKdag.col(sj);  
		}
		
		//update changes g_tau 
		void update(const int si, const int sj, const double gij, const double gji)
		{
			//update g_tau
			Eigen::RowVectorXd ri = uK.row(si) * g_tau - uK.row(si); 
			Eigen::RowVectorXd rj = uK.row(sj) * g_tau - uK.row(sj); 

			g_tau.noalias() -= (g_tau*uKdag.col(sj)) * ri/gij + (g_tau*uKdag.col(si)) * rj/gji; 

			/*
			if (itau_%blocksize_==0)
			{ //special treatment when update the block starting point 
											//otherwise this vertex will be untreated in stablization 
				//std::cout << "update: special treatment because update on the block boundary" << std::endl; 
				Vprop(si, sj, "L",  Storage_[itau_/blocksize_]);// update U in Storage 
			}
			*/
		}
		
		vertex generate_random_vertex()
		{
			int block = tpos / param.block_size;
			double tau = (block + rng()) * param.block_size;
			auto& b = lat.bonds("nearest neighbors")[rng() * 2. * lat.n_bonds()];
			return vertex{tau, b.first, b.second};
		}
		
		unsigned int select_random_vertex(vlist_t::iterator& vpos)
		{
			int block = tpos / param.block_size;
			vlist_t::iterator lower = std::lower_bound(vlist.begin(), vlist.end(), block*param.block_size, vertex::less()); 
			vlist_t::iterator upper = std::upper_bound(vlist.begin(), vlist.end(), (block+1)*param.block_size, vertex::less_equal());  //equal is exclude
			unsigned num_vertices = std::distance(lower, upper); //number of vertices in this block
			if (num_vertices == 0)
				vpos = vlist.end();
			else
				vpos = std::next(lower, rng() * num_vertices);
			return num_vertices;
		}
		
		double add_vertex(const vertex& v, bool compute_only_weight)
		{
			int block = v.tau / param.block_size;
			vlist_t::iterator lower = std::lower_bound(vlist.begin(), vlist.end(), block*param.block_size, vertex::less()); 
			vlist_t::iterator upper = std::upper_bound(vlist.begin(), vlist.end(), (block+1)*param.block_size, vertex::less_equal());  //equal is exclude
			unsigned num_vertices = std::distance(lower, upper); //number of vertices in this block
			
			wrap(v.tau);
			double G = gij(v.si, v.sj); // rotate it to real space 
			double ratio = -4.* G * G / (num_vertices + 1.); // gji = gij when they belongs to different sublattice 

			if(!compute_only_weight)
			{
				update(v.si, v.sj, G, G);
				insert_sorted<vertex, vertex::less>(vlist, v, vertex::less());
			}
			return ratio;
		}

		double remove_vertex(vlist_t::iterator vpos, bool compute_only_weight)
		{
			if (vlist.empty()) //since we will get a empty list as use this opputunity to reset all memory 
			{
				rebuild();
				return 0.;
			}
			
			if (vpos == vlist.end())
				return 0.;

			wrap(vpos->tau);
			
			double G = gij(vpos->si, vpos->sj);  
			double ratio = -4.* G * G; // gji = gij when they belongs to different sublattice 

			if(!compute_only_weight)
			{
				update(vpos->si, vpos->sj, G, G); 
				vlist.erase(vpos);
			}
			return ratio; 
		}
		
		void measure_static_observable(measurements& measure, const std::vector<std::string>& names,
			const std::vector<wick_static_base<matrix_t>>& obs,
			const std::vector<std::string>& vec_names,
			const std::vector<vector_wick_static_base<matrix_t>>& vec_obs)
		{
			matrix_t g_site = uK * g_tau * uKdag;
			
			for (int i = 0; i < names.size(); ++i)
				measure.add(names[i], obs[i].get_obs(g_site));
			for (int i = 0; i < vec_names.size(); ++i)
				measure.add(vec_names[i], vec_obs[i].get_obs(g_site));
		}
		
		void print_matrix(const matrix_t& m)
		{
			for (int i = 0; i < m.rows(); ++i)
			{
				for (int j = 0; j < m.cols(); ++j)
					std::cout << m(i, j) << " ";
				std::cout << std::endl;
			}
		}
	private:
		Random& rng;
		parameters& param;
		lattice& lat;
		
		vlist_t vlist;
		double tpos;
		
		matrix_t g_tau;
		matrix_t K;
		vector_t wK;
		matrix_t uK;
		matrix_t uKdag;
		matrix_t uKdagP;
		
		std::vector<matrix_t> storage;
};
