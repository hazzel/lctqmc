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
#include "wick_base.h"
#include "vector_wick_base.h"

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
		//using numeric_t = std::complex<double>;
		using vector_t = Eigen::Matrix<numeric_t, Eigen::Dynamic, 1>;
		using row_vector_t = Eigen::Matrix<numeric_t, 1, Eigen::Dynamic>;
		using matrix_t = Eigen::Matrix<numeric_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
		using diag_matrix_t = Eigen::DiagonalMatrix<numeric_t, Eigen::Dynamic>;
		using vlist_t = std::vector<vertex>;

		green_function(Random& rng_, parameters& param_, lattice& lat_)
			: rng(rng_), param(param_), lat(lat_), norm_error_sum(0.), norm_error_cnt(0)
		{}
		
		void set_K_matrix(const matrix_t& K)
		{
			Eigen::SelfAdjointEigenSolver<matrix_t> solver(K);
			wK = solver.eigenvalues(); 
			uK = solver.eigenvectors();
			uKdag = solver.eigenvectors().adjoint();
			for (int i=0; i < uK.rows(); ++i)
				uKcr.push_back(uKdag.col(i) * uK.row(i));
		}
		
		void set_trial_wf(const matrix_t& P)
		{
			uKdagP = uKdag * P;
			
			//Eigen::JacobiSVD<matrix_t> svd(uKdagP, Eigen::ComputeThinU);
			//uKdagP = svd.matrixU();
	
			qr_solver.compute(uKdagP);
			matrix_t p_q = matrix_t::Identity(uKdagP.rows(), uKdagP.cols());
			uKdagP = qr_solver.matrixQ() * p_q;
		}
		
		void initialize(double tau, const std::vector<vertex>& vlist_)
		{
			if (param.projective)
			{
				storage.resize(param.theta / param.block_size + 1);
				tpos = tau;
				vlist = vlist_;
				
				int block = tpos / param.block_size;
				storage[0] = uKdagP;
				
				for (int i=0; i < block; ++i)
				{
					matrix_t UR = storage[i];
					prop_from_left(-1, (i+1)*param.block_size, i*param.block_size, UR);

					//Eigen::JacobiSVD<matrix_t> svd(UR, Eigen::ComputeThinU); 
					//storage[i+1] = svd.matrixU();

					qr_solver.compute(UR);
					matrix_t p_q = matrix_t::Identity(UR.rows(), UR.cols());
					storage[i+1] = qr_solver.matrixQ() * p_q;
				}

				storage.back() = uKdagP.adjoint();
				for (int i = storage.size()-2; i>block; --i)
				{
					matrix_t VL = storage[i+1];
					prop_from_right(-1, (i+1)*param.block_size, i*param.block_size, VL);
					
					//Eigen::JacobiSVD<matrix_t> svd(VL, Eigen::ComputeThinV); 
					//storage[i] = svd.matrixV().adjoint();

					qr_solver.compute(VL.adjoint());
					matrix_t p_q = matrix_t::Identity(VL.rows(), VL.cols());
					storage[i] = p_q * qr_solver.matrixQ().adjoint();

					/*
					qr_solver.compute(VL);
					matrix_t r = qr_solver.matrixQR().triangularView<Eigen::Upper>();
					storage[i] = r * qr_solver.colsPermutation().transpose();
					for (int i = 0; i < storage[i].rows(); ++i)
						storage[i].row(i) = 1./qr_solver.matrixQR()(i, i) * storage[i].row(i);
					*/
				}
			}
			else
			{
				storage_U.resize(param.theta / param.block_size + 1);
				storage_D.resize(param.theta / param.block_size + 1);
				storage_V.resize(param.theta / param.block_size + 1);
				
				tpos = tau;
				vlist = vlist_;
				
				int block = tpos / param.block_size;
				
				storage_U[0] = matrix_t::Identity(lat.n_sites(), lat.n_sites());
				storage_D[0] = diag_matrix_t(lat.n_sites());
				storage_D[0].setIdentity();
				storage_V[0] = matrix_t::Identity(lat.n_sites(), lat.n_sites());
				for (int i=0; i<block; ++i)
				{
					storage_U[i+1] = storage_U[i];
					prop_from_left(-1, (i+1)*param.block_size, i*param.block_size, storage_U[i+1]);

					/*
					Eigen::JacobiSVD<matrix_t> svd(storage_U[i+1] * storage_D[i], Eigen::ComputeThinU | Eigen::ComputeThinV);
					storage_U[i+1] = svd.matrixU();
					storage_D[i+1] = svd.singularValues().asDiagonal();
					storage_V[i+1] = svd.matrixV().adjoint() * storage_V[i];
					*/
					
					qr_solver.compute(storage_U[i+1] * storage_D[i]);
					matrix_t R = qr_solver.matrixQR().triangularView<Eigen::Upper>();
					storage_U[i+1] = qr_solver.matrixQ();
					storage_D[i+1] = qr_solver.matrixQR().diagonal().asDiagonal();
					storage_V[i+1] = R * (qr_solver.colsPermutation().transpose() * storage_V[i]);
					for (int p = 0; p < storage_V[i+1].rows(); ++p)
						storage_V[i+1].row(p) = 1./storage_D[i+1].diagonal()[p] * storage_V[i+1].row(p);
				}
			
				storage_U.back() = matrix_t::Identity(lat.n_sites(), lat.n_sites());
				storage_D.back() = diag_matrix_t(lat.n_sites());
				storage_D.back().setIdentity();
				storage_V.back() = matrix_t::Identity(lat.n_sites(), lat.n_sites());
				for (int i = storage_U.size()-2; i>block; --i)
				{
					storage_V[i] = storage_V[i+1];
					prop_from_right(-1, (i+1)*param.block_size, i*param.block_size, storage_V[i]);

					/*
					Eigen::JacobiSVD<matrix_t> svd(storage_D[i+1] * storage_V[i], Eigen::ComputeThinU | Eigen::ComputeThinV); 
					storage_U[i] = storage_U[i+1]*svd.matrixU();
					storage_D[i] = svd.singularValues().asDiagonal();
					storage_V[i] = svd.matrixV().adjoint();
					*/
					
					qr_solver.compute(storage_D[i+1] * storage_V[i]);
					matrix_t Q = qr_solver.matrixQ();
					matrix_t R = qr_solver.matrixQR().triangularView<Eigen::Upper>();
					matrix_t U_r = storage_U[i];
					matrix_t D_r = storage_D[i];
					matrix_t V_r = storage_V[i];
					storage_V[i] = storage_U[i] * Q;
					storage_D[i] = qr_solver.matrixQR().diagonal().asDiagonal();
					storage_U[i] = R * qr_solver.colsPermutation().transpose();
					for (int p = 0; p < storage_U[i].rows(); ++p)
						storage_U[i].row(p) = 1./storage_D[i].diagonal()[p] * storage_U[i].row(p);
				}
			}
			g_tau = g_stable();
			std::cout << "tau = " << tpos << std::endl;
			std::cout << "gf:" << std::endl;
			matrix_t g_site = uK * g_tau * uKdag;
			print_matrix(g_site);
			std::cout << std::endl;
			std::cout << "exact" << std::endl;
			print_matrix(uK * g_exact() * uKdag);
		}
		
		unsigned int pert_order()
		{
			return vlist.size();
		}
		
		unsigned int pert_order(double tau)
		{
			vlist_t::iterator lower = std::lower_bound(vlist.begin(), vlist.end(), tau, vertex::less_equal());	//equal is exclude
			return std::distance(vlist.begin(), lower);
		}
		
		double tau()
		{
			return tpos;
		}

		double reset_norm_error()
		{
			double avg_norm_error = norm_error_sum / norm_error_cnt;
			norm_error_sum = 0.;
			norm_error_cnt = 0;
			return avg_norm_error;
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
			/*
			if (side == "L")
				A.noalias() -= 2.* uKdag.col(si) * (uK.row(si)* A) + 2.* uKdag.col(sj) * (uK.row(sj)* A); 
			else if (side == "R")
				A.noalias() -= 2.* (A*uKdag.col(si))* uK.row(si) + 2.* (A*uKdag.col(sj)) * uK.row(sj);
			*/
			if (side == "L")
				A.noalias() -= 2.* (uKcr[si] + uKcr[sj]) * A;
			else if (side == "R")
				A.noalias() -= 2.* A * (uKcr[si] + uKcr[sj]);
		}
		
		void rebuild()
		{
			
			matrix_t g_stab = g_stable();
			double err = ((g_tau - g_stab).cwiseAbs()).maxCoeff();
			g_tau = g_stab;
			
			
			/*
			matrix_t g_prev;
			calculate_gf(g_prev);
			stabilize();
			matrix_t g_stab;
			calculate_gf(g_stab);
			double err = ((g_prev - g_stab).cwiseAbs()).maxCoeff();
			*/
			
			if (err > 1E-6)
				std::cout << "Error (tau = " << tpos << "): " << err << std::endl;
			norm_error_sum += err;
			++norm_error_cnt;
		}
		
		void calculate_gf(matrix_t& m)
		{
			m = matrix_t::Identity(lat.n_sites(), lat.n_sites());
			m.noalias() -= (R_tau * W_tau) * L_tau;
		}
		
		matrix_t g_stable()
		{
			if (param.projective)
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
			else
			{
				int block = tpos / param.block_size;
				
				matrix_t U1 = storage_U[block];
				vector_t D1 = storage_D[block];
				matrix_t V1 = storage_V[block]; 
				prop_from_left(-1, tpos, block*param.block_size, U1);
			
				matrix_t U2 = storage_U[block+1];
				vector_t D2 = storage_D[block+1];
				matrix_t V2 = storage_V[block+1]; 
				prop_from_right(-1, (block+1)*param.block_size, tpos, V2);
				
				Eigen::JacobiSVD<matrix_t> svd((D1.asDiagonal()*V1)*(U2*D2.asDiagonal()), Eigen::ComputeThinU | Eigen::ComputeThinV);
				matrix_t U, D, V;
				U = U1 * svd.matrixU();
				D = svd.singularValues().asDiagonal();
				V = svd.matrixV().adjoint()*V2;

				matrix_t res = D;
				res.noalias() += U.inverse()*V.inverse();

				Eigen::JacobiSVD<matrix_t> svd2(res, Eigen::ComputeThinU | Eigen::ComputeThinV);
				D = svd2.singularValues().asDiagonal();
				return (svd2.matrixV().adjoint()*V).inverse() *  D.inverse() *  (U*svd2.matrixU()).inverse();
			}
		}
		
		void stabilize()
		{
			int block = tpos / param.block_size;
			
			R_tau = storage[block];
			prop_from_left(-1, tpos, block*param.block_size, R_tau);
		
			L_tau = storage[block+1];
			prop_from_right(-1, (block+1)*param.block_size, tpos, L_tau);
			
			W_tau = (L_tau * R_tau).inverse();
		}
		
		matrix_t g_exact()
		{
			matrix_t UR = uKdagP;
			prop_from_left(-1, tpos, 0, UR);
			matrix_t VL = uKdagP.adjoint();
			prop_from_right(-1, param.theta, tpos, VL);
			matrix_t id = matrix_t::Identity(lat.n_sites(), lat.n_sites());
			std::cout << "UR" << std::endl;
			print_matrix(UR);
			std::cout << std::endl;
			std::cout << "VL" << std::endl;
			print_matrix(VL);
			std::cout << std::endl;
			matrix_t v = matrix_t::Identity(lat.n_sites(), lat.n_sites());
			V_prop(vlist[0].si, vlist[0].sj, "R", v);
			print_matrix(v);
			std::cout << std::endl;
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
				
				//prop_from_left(-1, tau, tpos, R_tau);	// B(tau1) ... B(tau2) *U_
				//prop_from_left(1, tau, tpos, L_tau);	// V_ * B^{-1}(tau2) ... B^{-1}(tau1)
			}
			else
			{
				// B^{-1} G B 
				prop_from_right(1, tpos, tau, g_tau);	//  B^{-1}(tau2) ... B^{-1}(tau1) * U_
				prop_from_right(-1, tpos, tau, g_tau);	//  V_ * B(tau1) ... B(tau2)
				
				//prop_from_right(1, tpos, tau, R_tau);	//  B^{-1}(tau2) ... B^{-1}(tau1) * U_
				//prop_from_right(-1, tpos, tau, L_tau);	//  V_ * B(tau1) ... B(tau2)
			}
			tpos = tau;
			++param.wrap_refresh_cnt;
			
			//when we wrap to a new block we need to update storage 
			if (new_block > old_block)
			{	// move to a larger block on the left
				if (param.projective)
				{
					matrix_t UR = storage[old_block];
					prop_from_left(-1, new_block*param.block_size, old_block*param.block_size, UR);
					
					if (param.wrap_refresh_cnt >= param.wrap_refresh_interval)
					{
						//Eigen::JacobiSVD<matrix_t> svd(UR, Eigen::ComputeThinU);
						//UR = svd.matrixU();

						qr_solver.compute(UR);
						matrix_t p_q = matrix_t::Identity(UR.rows(), UR.cols());
						UR = qr_solver.matrixQ() * p_q;
						
						param.wrap_refresh_cnt = 0;
					}
					storage[new_block] = UR;
				}
				else
				{
					storage_U[new_block] = storage_U[old_block];
					storage_D[new_block] = storage_D[old_block];
					storage_V[new_block] = storage_V[old_block];
					prop_from_left(-1, new_block*param.block_size, old_block*param.block_size, storage_U[new_block]);
					
					if (param.wrap_refresh_cnt >= param.wrap_refresh_interval)
					{
						Eigen::JacobiSVD<matrix_t> svd(storage_U[new_block] * storage_D[new_block], Eigen::ComputeThinU | Eigen::ComputeThinV);
						storage_U[new_block] = svd.matrixU();
						storage_D[new_block] = svd.singularValues().asDiagonal();
						storage_V[new_block] = svd.matrixV().adjoint() * storage_V[new_block];
						
						param.wrap_refresh_cnt = 0;
					}
				}
			}
			else if (new_block < old_block)
			{	// move to smaller block
				if (param.projective)
				{
					matrix_t VL = storage[old_block+1];
					prop_from_right(-1, (old_block+1)*param.block_size, old_block*param.block_size, VL);
					
					if (param.wrap_refresh_cnt >= param.wrap_refresh_interval)
					{
						//Eigen::JacobiSVD<matrix_t> svd(VL, Eigen::ComputeThinV);
						//VL = svd.matrixV().adjoint();

						qr_solver.compute(VL.adjoint());
						matrix_t p_q = matrix_t::Identity(VL.rows(), VL.cols());
						VL = p_q * qr_solver.matrixQ().adjoint();

						//qr_solver.compute(VL);
						//matrix_t r = qr_solver.matrixQR().triangularView<Eigen::Upper>();
						//storage[new_block+1] = r * qr_solver.colsPermutation().transpose();
						//for (int i = 0; i < storage[new_block+1].rows(); ++i)
						//	storage[new_block+1].row(i) = 1./qr_solver.matrixQR()(i, i) * storage[new_block+1].row(i);
						
						param.wrap_refresh_cnt = 0;
					}
					storage[new_block+1] = VL;
				}
				else
				{
					storage_U[new_block+1] = storage_U[old_block+1];
					storage_D[new_block+1] = storage_D[old_block+1];
					storage_V[new_block+1] = storage_V[old_block+1];
					prop_from_right(-1, (old_block+1)*param.block_size, old_block*param.block_size, storage_V[new_block+1]);
					
					if (param.wrap_refresh_cnt >= param.wrap_refresh_interval)
					{
						Eigen::JacobiSVD<matrix_t> svd(storage_D[new_block+1] * storage_V[new_block+1], Eigen::ComputeThinU | Eigen::ComputeThinV);
						storage_U[new_block+1] = storage_U[new_block+1]*svd.matrixU();
						storage_D[new_block+1] = svd.singularValues().asDiagonal();
						storage_V[new_block+1] = svd.matrixV().adjoint();
						
						param.wrap_refresh_cnt = 0;
					}
				}
			}
		}
		
		numeric_t gij(const int si, const int sj)
		{	// current g in the site basis 
			// (U gtau U^{dagger} )_ij 
			
			uK_si_gtau = uK.row(si) * g_tau;
			return  uK_si_gtau * uKdag.col(sj);
			
			//return  (uK.row(si) * g_tau) * uKdag.col(sj);
			
			/*
			uK_si_R_W = (uK.row(si) * R_tau) * W_tau;
			L_uKdag_sj = L_tau * uKdag.col(sj);
			return (si==sj) ? 1.0 : 0.0 - (uK_si_R_W * L_uKdag_sj);
			
			//return (si==sj) ? 1.0 : 0.0 - (uK.row(si)*R_tau) *  W_tau * (L_tau*uKdag.col(sj));  
			*/
		}
		
		//update changes g_tau 
		void update(const int si, const int sj, const numeric_t& gij, const numeric_t& gji)
		{
			
			//update g_tau
			row_vector_t ri = uK_si_gtau - uK.row(si);
			row_vector_t rj = -uK.row(sj);
			rj.noalias() += uK.row(sj) * g_tau;
			g_tau.noalias() -= (g_tau*uKdag.col(sj)) * ri / gij + (g_tau*uKdag.col(si)) * rj / gji;
			
			/*
			vector_t W_L_uKdag_sj = W_tau * L_uKdag_sj;
			vector_t L_uKdag_si = L_tau * uKdag.col(si);
			vector_t W_L_uKdag_si = W_tau * L_uKdag_si;
			row_vector_t uK_sj_R_W = (uK.row(sj) * R_tau) * W_tau;
			W_tau.noalias() += 1./gij * (W_L_uKdag_sj * uK_si_R_W);
			W_tau.noalias() += 1./gji * (W_L_uKdag_si * uK_sj_R_W);
			
			//W_tau.noalias() += ((W_tau * L_uKdag_sj) * uK_si_R_W) / gij + ((W_tau * (L_tau*uKdag.col(si))) * ((uK.row(sj)*R_tau)*W_tau)) / gji;
			
			//W_tau.noalias() += ((W_tau * (L_tau*uKdag.col(sj))) * ((uK.row(si)*R_tau)*W_tau)) / gij + ((W_tau * (L_tau*uKdag.col(si))) * ((uK.row(sj)*R_tau)*W_tau)) / gji;
			V_prop(si, sj, "L",  R_tau);
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
			numeric_t G = gij(v.si, v.sj); // rotate it to real space 
			double ratio = -4.* std::abs(G) * std::abs(G) / (num_vertices + 1.); // gji = gij when they belongs to different sublattice 

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
			
			numeric_t G = gij(vpos->si, vpos->sj);
			double ratio = -4.* std::abs(G) * std::abs(G); // gji = gij when they belongs to different sublattice 

			if(!compute_only_weight)
			{
				update(vpos->si, vpos->sj, G, G); 
				vlist.erase(vpos);
			}
			return ratio; 
		}

		double shift_vertex(vlist_t::iterator vpos, int si_prime, int sj_prime, bool compute_only_weight)
		{
			if (vlist.empty()) //since we will get a empty list as use this opputunity to reset all memory 
			{
				rebuild();
				return 0.;
			}
			
			if (vpos == vlist.end())
				return 0.;

			wrap(vpos->tau);
			
			numeric_t G = gij(vpos->sj, sj_prime);
			double ratio = 4.* std::abs(G) * std::abs(G); // gji = gij when they belongs to different sublattice 

			if(!compute_only_weight)
			{
				update(vpos->sj, sj_prime, G, -G);
				vpos->sj = sj_prime;
			}
			return ratio; 
		}
		
		void measure_static_observables(measurements& measure, const std::vector<std::string>& names,
			const std::vector<wick_static_base<matrix_t>>& obs,
			const std::vector<std::string>& vec_names,
			const std::vector<vector_wick_static_base<matrix_t>>& vec_obs)
		{
			matrix_t g_site = uK * g_tau * uKdag;
			
			/*
			matrix_t g_site;
			calculate_gf(g_site);
			g_site = uK * g_site * uKdag;
			*/
			
			for (int i = 0; i < names.size(); ++i)
				measure.add(names[i], obs[i].get_obs(g_site));
			for (int i = 0; i < vec_names.size(); ++i)
				measure.add(vec_names[i], vec_obs[i].get_obs(g_site));
		}
		
		void measure_dynamical_observables(std::vector<std::vector<double>>& dyn_tau, const std::vector<std::string>& names,
			const std::vector<wick_base<matrix_t>>& obs,
			const std::vector<std::string>& vec_names,
			const std::vector<vector_wick_base<matrix_t>>& vec_obs)
		{
			double tpos_buffer = tpos;
			matrix_t g_tau_buffer = g_tau;
			//matrix_t R_tau_buffer = R_tau;
			//matrix_t L_tau_buffer = L_tau;
			//matrix_t W_tau_buffer = W_tau;
			std::vector<matrix_t> storage_buffer = storage;
			std::vector<matrix_t> storage_buffer_U = storage_U;
			std::vector<diag_matrix_t> storage_buffer_D = storage_D;
			std::vector<matrix_t> storage_buffer_V = storage_V;
			
			std::vector<matrix_t> et_gf_L(param.dyn_tau_steps/2+1);
			std::vector<matrix_t> et_gf_R(param.dyn_tau_steps/2+1);
			matrix_t id = matrix_t::Identity(lat.n_sites(), lat.n_sites());
			
			matrix_t time_displaced_gf = id;
			
			if (tpos > param.theta/2)
				param.direction = -1;
			else
				param.direction = 1;

			et_gf_L[0] = g_tau;
			//calculate_gf(et_gf_L[0]);
			
			for (int n = 1; n <= param.dyn_tau_steps/2; ++n)
			{
				wrap(tpos + param.direction * param.dyn_delta_tau);
				//stabilize();
				if (n % param.wrap_refresh_interval == 0)
					rebuild();
				et_gf_L[n] = g_tau;
				//calculate_gf(et_gf_L[n]);
			}

			matrix_t& et_gf_0 = et_gf_L[param.dyn_tau_steps/2];
			et_gf_R[0] = et_gf_0;
			for (int n = 1; n <= param.dyn_tau_steps/2; ++n)
			{
				wrap(tpos + param.direction * param.dyn_delta_tau);
				//stabilize();
				if (n % param.wrap_refresh_interval == 0)
					rebuild();
				et_gf_R[n] = g_tau;
				//calculate_gf(et_gf_R[n]);
			}
			
			get_obs_values(dyn_tau, 0, et_gf_0, et_gf_0, et_gf_0, obs, vec_obs);
			for (int n = 1; n <= param.dyn_tau_steps/2; ++n)
			{
				if (param.direction == -1)
				{
					matrix_t g_l = et_gf_L[et_gf_L.size() - n];
					prop_from_left(-1, param.theta/2+param.block_size/2 + n*param.dyn_delta_tau, param.theta/2+param.block_size/2 + (n-1)*param.dyn_delta_tau, g_l);
					matrix_t g_r = et_gf_R[n];
					prop_from_left(-1, param.theta/2+param.block_size/2 - (n-1)*param.dyn_delta_tau, param.theta/2+param.block_size/2 - n*param.dyn_delta_tau, g_r);

					time_displaced_gf = g_l * time_displaced_gf;
					get_obs_values(dyn_tau, 2*n-1, et_gf_R[n-1], et_gf_L[et_gf_L.size() - n - 1], time_displaced_gf, obs, vec_obs);
					
					time_displaced_gf = time_displaced_gf * g_r;
					get_obs_values(dyn_tau, 2*n, et_gf_R[n], et_gf_L[et_gf_L.size() - n - 1], time_displaced_gf, obs, vec_obs);
				}
				else
				{
					matrix_t g_l = et_gf_R[n];
					prop_from_right(-1, param.theta/2+param.block_size/2 + n*param.dyn_delta_tau, param.theta/2+param.block_size/2 + (n-1)*param.dyn_delta_tau, g_l);
					matrix_t g_r = et_gf_L[et_gf_L.size() - n];
					prop_from_right(-1, param.theta/2+param.block_size/2 - (n-1)*param.dyn_delta_tau, param.theta/2+param.block_size/2 - n*param.dyn_delta_tau, g_r);
					
					time_displaced_gf = g_l * time_displaced_gf;
					get_obs_values(dyn_tau, 2*n-1, et_gf_L[et_gf_L.size() - n], et_gf_R[n], time_displaced_gf, obs, vec_obs);
					
					time_displaced_gf = time_displaced_gf * g_r;
					get_obs_values(dyn_tau, 2*n, et_gf_L[et_gf_L.size() - n - 1], et_gf_R[n], time_displaced_gf, obs, vec_obs);
				}
			}
			
			tpos = tpos_buffer;
			g_tau = g_tau_buffer;
			//R_tau = R_tau_buffer;
			//L_tau = L_tau_buffer;
			//W_tau = W_tau_buffer;
			storage = storage_buffer;
			storage_U = storage_buffer_U;
			storage_D = storage_buffer_D;
			storage_V = storage_buffer_V;
		}

		std::vector<double> measure_Hv_tau()
		{
			std::vector<double> hv_tau(param.ep_tau_steps, 0.);
			auto lower = std::lower_bound(vlist.begin(), vlist.end(), param.theta/2.-param.ep_window/2., vertex::less()); 
			auto upper = std::upper_bound(vlist.begin(), vlist.end(), param.theta/2.+param.ep_window/2., vertex::less_equal());  //equal is exclude
			
			int i_lower = std::distance(vlist.begin(), lower);
			int i_upper = std::distance(vlist.begin(), upper);
			
			int i = rng() * (i_upper - i_lower) + i_lower;
			for(int j = i_lower; j < i_upper; ++j)
			{
				if (i == j)
					continue;
				double delta_tau = vlist[j].tau - vlist[i].tau;
				int tau_block;
				if (delta_tau > 0.)
					tau_block = delta_tau / param.ep_delta_tau;
				else
					tau_block = (param.theta + delta_tau) / param.ep_delta_tau;
				int k = i_upper - i_lower;
				if (tau_block < param.ep_tau_steps)
					hv_tau[tau_block] += 1. / param.V / param.V / param.ep_delta_tau * param.ep_window / param.theta;
			}
				
			/*
			for(int i = i_lower; i < i_upper; ++i)
			//int i = rng() * (i_upper - i_lower) + i_lower;
				for(int j = i+1; j < i_upper; ++j)
				{
					double delta_tau = std::abs(vlist[i].tau - vlist[j].tau);
					int tau_block;
					if (param.trial_wave_function == "t_only")
						tau_block = (delta_tau < param.ep_window/2.) ? delta_tau / param.ep_delta_tau : (param.ep_window - delta_tau) / param.ep_delta_tau;
					else
						tau_block = delta_tau / param.ep_delta_tau;
					int k = i_upper - i_lower;
					if (tau_block < param.ep_tau_steps)
						hv_tau[tau_block] += 1. / param.V / param.V / param.ep_delta_tau * param.ep_window / param.theta / k;
				}
			*/

			return hv_tau;
		}
		
		void get_obs_values(std::vector<std::vector<double>>& dyn_tau, int tau,
			const matrix_t& et_gf_0, const matrix_t& et_gf_t, const matrix_t& td_gf, 
			const std::vector<wick_base<matrix_t>>& obs, const std::vector<vector_wick_base<matrix_t>>& vec_obs)
		{
			matrix_t et_gf_0_site = uK * et_gf_0 * uKdag;
			matrix_t et_gf_t_site = uK * et_gf_t * uKdag;
			matrix_t td_gf_site = uK * td_gf * uKdag;
			for (int i = 0; i < obs.size(); ++i)
				dyn_tau[i][tau] = obs[i].get_obs(et_gf_0_site, et_gf_t_site, td_gf_site);
			int cnt = 0;
			for (int i = 0; i < vec_obs.size(); ++i)
			{
				auto& values = vec_obs[i].get_obs(et_gf_0_site, et_gf_t_site, td_gf_site);
				for (int j = 0; j < vec_obs[i].n_values; ++j)
				{
					dyn_tau[obs.size()+cnt][tau] = values[j];
					++cnt;
				}
			}
		}
		
		void serialize(odump& out)
		{
			out.write(param.static_measure_cnt);
			out.write(tpos);
			int size = vlist.size();
			out.write(size);
			for (auto& v : vlist)
			{
				out.write(v.tau);
				out.write(v.si);
				out.write(v.sj);
			}
		}

		void serialize(idump& in)
		{
			in.read(param.static_measure_cnt);
			in.read(tpos);
			int size;
			in.read(size);
			for (int i = 0; i < size; ++i)
			{
				double tau;
				in.read(tau);
				int si;
				in.read(si);
				int sj;
				in.read(sj);
				vlist.push_back({tau, si, sj});
			}
			std::sort(vlist.begin(), vlist.end(), vertex::less());
			initialize(tpos, vlist);
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
		Eigen::ColPivHouseholderQR<matrix_t> qr_solver;
		
		vlist_t vlist;
		double tpos;
		double norm_error_sum;
		int norm_error_cnt;
		
		matrix_t g_tau;
		matrix_t L_tau;
		matrix_t W_tau;
		matrix_t R_tau;
		
		row_vector_t uK_si_gtau;
		row_vector_t uK_si_R_W;
		vector_t L_uKdag_sj;
		
		vector_t wK;
		matrix_t uK;
		matrix_t uKdag;
		std::vector<matrix_t> uKcr;
		matrix_t uKdagP;
		
		std::vector<matrix_t> storage;
		std::vector<matrix_t> storage_U;
		std::vector<diag_matrix_t> storage_D;
		std::vector<matrix_t> storage_V;
};
