#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <cassert>
#include <vector>
#include "lattice.h"

struct vertex
{
	double tau;
	unsigned int si;
	unsigned int sj;
	
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

		green_function(lattice& lat_)
			: lat(lat_)
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
			
			//upper > tau1 > lower > tau2 
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
				for (unsigned l=0; l<lat.n_sites(); ++l) 
					A.row(l) *= std::exp(sign * tau * wK(l));
			}
			else if (side == "R")
			{	// A * exp(tau * w)
				for (unsigned l=0; l<lat.n_sites(); ++l) 
					A.col(l) *= std::exp(sign * tau * wK(l));
			}
		}

		void V_prop(unsigned int si, unsigned int sj, const std::string& side,  matrix_t& A) const // A*U^{dagger} V U or U^{dagger} V U * A 
		{
			if (side == "L")
				A.noalias() -= 2.* uKdag.col(si) * (uK.row(si)* A) + 2.* uKdag.col(sj) * (uK.row(sj)* A); 
			else if (side == "R")
				A.noalias() -= 2.* (A*uKdag.col(si))* uK.row(si) + 2.* (A*uKdag.col(sj)) * uK.row(sj);
		}

		void initialize(const std::vector<vertex>& vlist_, double theta_, double block_size_)
		{
			block_size = block_size_;
			theta = theta_;
			storage.resize(theta / block_size + 1);

			vlist = vlist_;
			tpos = vlist.begin()->tau;
			
			int block = tpos / block_size;
			storage[0] = uKdagP;
			
			for (int i=0; i < block; ++i)
			{
				matrix_t UR = storage[i];
				prop_from_left(-1, (i+1)*block_size, i*block_size, UR);

				Eigen::JacobiSVD<matrix_t> svd(UR, Eigen::ComputeThinU); 
				storage[i+1] = svd.matrixU();
			}

			storage.back() = uKdagP.adjoint();
			for (int i = storage.size()-2; i>block; --i)
			{
				matrix_t VL = storage[i+1];
				prop_from_right(-1, (i+1)*block_size, i*block_size, VL);

				Eigen::JacobiSVD<matrix_t> svd(VL, Eigen::ComputeThinV); 
				storage[i] = svd.matrixV().adjoint();
			}

			g_tau = g_stable();
			measure();
		}
		
		matrix_t g_stable()
		{
			int block = tpos / block_size;

			matrix_t UR = storage[block];
			prop_from_left(-1, tpos, block*block_size, UR);
		
			matrix_t VL = storage[block+1];
			prop_from_right(-1, (block+1)*block_size, tpos, VL);
			
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
			prop_from_right(-1, theta, tpos, VL);
			matrix_t id = matrix_t::Identity(lat.n_sites(), lat.n_sites());
			return id - UR * (VL * UR).inverse() * VL;
		}
		
		void wrap(double tau)
		{
			int old_block = tpos / block_size;
			int new_block = tau / block_size;
			
			std::cout << "wrap from tau = " << tpos << " to " << tau << ", from block " << old_block << " to " << new_block << std::endl;
			
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
				prop_from_left(-1, new_block*block_size, old_block*block_size, storage[old_block]);
				Eigen::JacobiSVD<matrix_t> svd(storage[old_block], Eigen::ComputeThinU); 
				storage[new_block] = svd.matrixU();
			}
			else if (new_block < old_block)
			{// move to smaller block 
				prop_from_right(-1, (old_block+1)*block_size, old_block*block_size, storage[old_block+1]);
				Eigen::JacobiSVD<matrix_t> svd(storage[old_block+1], Eigen::ComputeThinV);
				storage[new_block+1] = svd.matrixV().adjoint();
			}
		}
		
		void measure()
		{
			matrix_t g_site = uK * g_tau * uKdag;
			
			std::cout << "g_site at tau = " << tpos << ":" << std::endl;
			print_matrix(g_site);
			std::cout << "g_exact at tau = " << tpos << ":" << std::endl;
			g_site = uK * g_exact() * uKdag;
			print_matrix(g_site);
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
		lattice& lat;
		
		double theta;
		double block_size;
		
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
