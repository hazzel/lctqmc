#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <vector>

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
	};

	struct less_equal
	{
		constexpr bool operator()(const vertex& lhs, const vertex& rhs) const
		{
			return lhs.tau <= rhs.tau;
		}
	};
};  

class green_function
{
	public:
		using numeric_t = double;
		using vector_t = Eigen::Vector<numeric_t, Eigen::Dynamic>;
		using matrix_t = Eigen::Matrix<numeric_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
		using vlist_t = std::vector<vertex>;

		green_function(const matrix_t& K_)
			: K(K_)
		{
			Eigen::SelfAdjointEigenSolver<matrix_t> solver(K);
			wK = solver.eigenvalues();  
			uK = solver.eigenvectors(); 
			uKdag = solver.eigenvectors().adjoint();
		}
		
		// it can do  B(tau_m)... B(tau_n) * A  when sign = -1
		// or        A* B(tau_n)^{-1} ... B(tau_m)^{-1} when sign = 1 
		// Btau(tau_n) does not contain vertex contribution  
		void prop_from_left(const int sign, const double tau_m, const double tau_n, matrix_t& A) const
		{
			std::assert(tau_m >= tau_n);

			std::string side = (sign == -1) ? "L" : "R";
		
			//since we always point to the right of a vertex, we push downwards the lower boundary 
			//this is different with block convention (,]  
			vlist_t::const_iterator lower = std::lower_bound(vlist.begin(), vlist.end(), tau_n, vertex::less_equal());	//equal is exclude
			vlist_t::const_iterator upper = std::upper_bound(vlist.begin(), vlist.end(), tau_m, vertex::less());

			//upper > tau_m > lower > tau_n
			if (lower == upper)	// there is no vertex in between tau1 and tau2 
				Kprop(sign, tau_m - tau_n, side, A);
			else
			{
				Kprop(sign, (*lower).tau - tau_n, side, A);
				for (vlist_t::const_iterator it = lower; it != upper; ++it)
				{
					Vprop(it->si, it->sj, side, A);
					if (upper - it > 1)
						Kprop(sign, std::next(it)->tau - it->tau, side, A);
					else
						Kprop(sign, tau_m - it->tau, side, A);
				}
			}
		}

		// it can do A* B(tau_m)... B(tau_n)  when sign = -1
		// or        B(tau_n)^{-1} ... B(tau_m)^{-1}* A when sign = 1 
		void prop_from_right(const int sign, const itime_type tau_m, const itime_type tau_n, matrix_t& A) const
		{ 
			std::assert(tau_m >= tau_n);

			std::string side = (sign == 1) ? "L" : "R";
	
			vlist_t::const_iterator lower = std::lower_bound(vlist.begin(), vlist.end(), tau_n, vertex::less_equal());	//equal is exclude
			vlist_t::const_iterator upper = std::upper_bound(vlist.begin(), vlist.end(), tau_m, vertex::less());
			
			//upper > tau1 > lower > tau2 
			if (lower == upper)	// there is no vertex in between tau1 and tau2 
				Kprop(sign, tau_m - tau_n, side, A); 
			else
			{
				Kprop(sign, tau_m - *(--upper), side, A);

				for (vlist_type::const_iterator it = upper; it != lower; --it)
				{
					Vprop(it->si, it->sj, side, A);
					Kprop(sign, it->tau - std::prev(it)->tau, side, A); 
				}
				//the last step by hand (it1 = lower, it2 point to tau2) 
				Vprop(lower->si, lower->sj, side, A);
				Kprop(sign, lower->tau - tau_n , side, A); 
			}
		}

		void K_prop(double tau, const std::string& side, matrix_t& A) const
		{
			if (std::abs(tau) < 1E-15)
					return;

			if (side == "L")
			{	//  exp(tau * w) * A 
				for (unsigned l=0; l<l.n_sites(); ++l) 
					A.row(l) *= std::exp(tau * wK(l));
			}
			else if (side == "R")
			{	// A * exp(tau * w)
				for (unsigned l=0; l<l.n_sites(); ++l) 
					A.col(l) *= std::exp(tau * wK(l));
			}
		}

		void V_prop(unsigned int si, unsigned int sj, const std::string& side,  matrix_t& A) const // A*U^{dagger} V U or U^{dagger} V U * A 
		{
			if (side == "L")
				A.noalias() -= 2.* uKdag.col(si) * (uK.row(si)* A) + 2.* uKdag.col(sj) * (uK.row(sj)* A); 
			else if (side == "R")
				A.noalias() -= 2.* (A*uKdag.col(si))* uK.row(si) + 2.* (A*uKdag.col(sj)) * uK.row(sj);
		}

		void initialize(const std::vector<vertex>& vlist_)
		{
			vlist = vlist_;
			
			int block = vpos->tau / n_block_size;
			storage[0] = uKdagP;
			for (int i=0; i < block; ++i)
			{
				matrix_t UR = storage[i];
				prop_from_left(-1, (i+1)*n_block_size, i*n_block_size, UR);

				Eigen::JacobiSVD<matrix_t> svd(UR, Eigen::ComputeThinU); 
				storage[i+1] = svd.matrixU(); 
			}

			storage.back() = uKdagP_.adjoint();
			for (int i = storage.size()-2; i>b; --i)
			{
				matrix_t VL = storage[i+1];
				prop_from_right(-1, (i+1)*n_block_size, i*n_block_size, VL);

				Eigen::JacobiSVD<Mat> svd(VL, Eigen::ComputeThinV); 
				storage[i] = svd.matrixV().adjoint();
			}

			g_tau = g_stable(); // from scratch
		}
		
		matrix_t g_stable()
		{
			int block = vpos->tau / n_block_size;

			matrix_t UR = Storage_[b];
			prop_from_left(-1, vpos->tau, block*n_block_size, UR);
		
			matrix_t VL = Storage_[b+1];
			prop_from_right(-1, (block+1)*n_block_size, vpos->tau, VL);
			
			matrix_t res = -UR * ((VL*UR).inverse() * VL);
			for (int l = 0; l < l.n_sites(); ++l)
				res(l, l) += 1.0;

			return res; 
		}
		
		void wrap() {}
		void measure() {}
	private:
		unsigned int n_matrix_size;
		unsigned int n_block_size;
		
		vlist_t vlist;
		vlist_t::iterator vpos;
		
		matrix_t g_tau;
		matrix_t K;
		vector_t wK;
		matrix_t uK;
		matrix_t uKdag;
		matrix_t uKdagP;
		
		std::vector<matrix_t> storage;
};
