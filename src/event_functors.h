#pragma once
#include <map>
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
#include "wick_functors.h"

struct event_set_trial_wf
{
	Random& rng;
	parameters& param;
	lattice& lat;
	green_function& gf;
	
	typedef green_function::matrix_t matrix_t;
	
	matrix_t symmetrize_EV(const matrix_t& S, const Eigen::VectorXd& en, const matrix_t& pm)
	{
		matrix_t S_s = S + pm * S;
		matrix_t S_a = S - pm * S;
		matrix_t S_so(lat.n_sites(), S_s.cols());
		matrix_t S_ao(lat.n_sites(), S_s.cols());
		matrix_t S_f = matrix_t::Zero(lat.n_sites(), 2*S_s.cols());

		for (int i = 0; i < S_s.cols(); ++i)
		{
			if (S_s.col(i).norm() > param.epsilon)
				S_s.col(i) /= S_s.col(i).norm();
			else
				S_s.col(i) *= 0.;
			if (S_a.col(i).norm() > param.epsilon)
				S_a.col(i) /= S_a.col(i).norm();
			else
				S_a.col(i) *= 0.;
		}

		int cnt = 0;
		for (int i = 0; i < S_s.cols(); ++i)
		{
			int j;
			for (j = i; j < S_s.cols() && std::abs(en(j)-en(i)) < param.epsilon ; ++j)
			{
				S_so.col(j) = S_s.col(j);
				S_ao.col(j) = S_a.col(j);
				for (int k = i; k < j; ++k)
				{
					S_so.col(j) -= S_so.col(k) * (S_so.col(k).adjoint() * S_s.col(j));
					S_ao.col(j) -= S_ao.col(k) * (S_ao.col(k).adjoint() * S_a.col(j));
				}
				//std::cout << "E=" << en(i) << ", orth: i=" << i << ", j=" << j << ": " << S_so.col(j).norm() << " " << S_ao.col(j).norm() << std::endl;
				if (S_so.col(j).norm() > param.epsilon)
				{
					S_so.col(j) /= S_so.col(j).norm();
					S_f.col(cnt) = S_so.col(j);
					++cnt;
				}
				if (S_ao.col(j).norm() > param.epsilon)
				{
					S_ao.col(j) /= S_ao.col(j).norm();
					S_f.col(cnt) = S_ao.col(j);
					++cnt;
				}
			}
			i = j - 1;
		}
		if (cnt != S.cols())
		{
			std::cout << "Error! Found " << cnt << " out of " << 2*S.cols() << std::endl;
			throw(std::runtime_error("Error in symmetrization. Wrong number of states."));
		}
		return S_f.leftCols(S.cols());
	}
	
	matrix_t ph_symmetrize_EV(const matrix_t& S, const matrix_t& pm, const matrix_t& inv_pm)
	{
		matrix_t S_s(lat.n_sites(), 2), S_a(lat.n_sites(), 2);
		S_s.col(0) = S.col(0) + S.col(1);
		S_s.col(1) = S.col(2) + S.col(3);
		S_a.col(0) = S.col(0) - S.col(1);
		S_a.col(1) = S.col(2) - S.col(3);
		matrix_t S_so(lat.n_sites(), S_s.cols());
		matrix_t S_ao(lat.n_sites(), S_s.cols());
		matrix_t S_sf = matrix_t::Zero(lat.n_sites(), S_s.cols());
		matrix_t S_af = matrix_t::Zero(lat.n_sites(), S_s.cols());

		for (int i = 0; i < S_s.cols(); ++i)
		{
			if (S_s.col(i).norm() > param.epsilon)
				S_s.col(i) /= S_s.col(i).norm();
			else
				S_s.col(i) *= 0.;
			if (S_a.col(i).norm() > param.epsilon)
				S_a.col(i) /= S_a.col(i).norm();
			else
				S_a.col(i) *= 0.;
		}

		int cnt_s = 0, cnt_a = 0;
		for (int i = 0; i < S_s.cols(); ++i)
		{
			int j;
			for (j = i; j < S_s.cols(); ++j)
			{
				S_so.col(j) = S_s.col(j);
				S_ao.col(j) = S_a.col(j);
				for (int k = i; k < j; ++k)
				{
					S_so.col(j) -= S_so.col(k) * (S_so.col(k).adjoint() * S_s.col(j));
					S_ao.col(j) -= S_ao.col(k) * (S_ao.col(k).adjoint() * S_a.col(j));
				}
				//std::cout << "orth: i=" << i << ", j=" << j << ": " << S_so.col(j).norm() << " " << S_ao.col(j).norm() << std::endl;
				if (S_so.col(j).norm() > param.epsilon)
				{
					S_so.col(j) /= S_so.col(j).norm();
					S_sf.col(cnt_s) = S_so.col(j);
					++cnt_s;
				}
				if (S_ao.col(j).norm() > param.epsilon)
				{
					S_ao.col(j) /= S_ao.col(j).norm();
					S_af.col(cnt_a) = S_ao.col(j);
					++cnt_a;
				}
			}
			i = j - 1;
		}
		matrix_t S_f(lat.n_sites(), cnt_s + cnt_a);
		S_f.leftCols(cnt_s) = S_sf.leftCols(cnt_s);
		S_f.rightCols(cnt_a) = S_af.leftCols(cnt_a);
		return S_f.leftCols(S.cols());
	}
	
	std::vector<std::vector<int>> get_energy_blocks(const Eigen::VectorXd& en)
	{
		std::vector<std::vector<int>> energy_blocks;
		energy_blocks.push_back({0, lat.n_sites()-1});
		for (int i = 1; i < lat.n_sites()/2; ++i)
		{
			if (std::abs(en(i) - en(energy_blocks.back()[0])) > param.epsilon)
				energy_blocks.push_back(std::vector<int>());
			energy_blocks.back().push_back(i);
			energy_blocks.back().push_back(lat.n_sites()-1-i);
		}
		for (int i = 0; i < energy_blocks.size(); ++i)
			std::sort(energy_blocks[i].begin(), energy_blocks[i].end());
		return energy_blocks;
	}
	
	matrix_t symmetrize_ph_blocks(const matrix_t& S, const std::vector<std::vector<int>>& energy_blocks, const Eigen::VectorXd& en, const matrix_t& pm)
	{
		matrix_t S_s = S + pm * S;
		matrix_t S_a = S - pm * S;
		matrix_t S_so(lat.n_sites(), S_s.cols());
		matrix_t S_ao(lat.n_sites(), S_s.cols());
		matrix_t S_f = matrix_t::Zero(lat.n_sites(), 2*S_s.cols());

		for (int i = 0; i < S_s.cols(); ++i)
		{
			if (S_s.col(i).norm() > param.epsilon)
				S_s.col(i) /= S_s.col(i).norm();
			else
				S_s.col(i) *= 0.;
			if (S_a.col(i).norm() > param.epsilon)
				S_a.col(i) /= S_a.col(i).norm();
			else
				S_a.col(i) *= 0.;
		}

		int cnt_s = 0, cnt_a = 0;
		for (int i = 0; i < energy_blocks.size(); ++i)
			for (int j = 0; j < energy_blocks[i].size(); ++j)
			{
				int b = energy_blocks[i][j];
				S_so.col(b) = S_s.col(b);
				S_ao.col(b) = S_a.col(b);
				for (int k = 0; k < j; ++k)
				{
					int a = energy_blocks[i][k];
					S_so.col(b) -= S_so.col(a) * (S_so.col(a).adjoint() * S_s.col(b));
					S_ao.col(b) -= S_ao.col(a) * (S_ao.col(a).adjoint() * S_a.col(b));
				}
				std::cout << "E=" << en(b) << ", orth: i=" << i << ", b=" << b << ": " << S_so.col(b).norm() << " " << S_ao.col(b).norm() << std::endl;
				if (S_so.col(b).norm() > param.epsilon)
				{
					S_so.col(b) /= S_so.col(b).norm();
					S_f.col(cnt_s) = S_so.col(b);
					++cnt_s;
				}
				if (S_ao.col(b).norm() > param.epsilon)
				{
					S_ao.col(b) /= S_ao.col(b).norm();
					S_f.col(lat.n_sites()/2+cnt_a) = S_ao.col(b);
					++cnt_a;
				}
			}
		if (cnt_s != S.cols()/2 || cnt_a != S.cols()/2)
		{
			std::cout << "Error! Found " << cnt_s << " symmetric states out of " << 2*S.cols() << std::endl;
			std::cout << "Error! Found " << cnt_a << " antisymmetric states out of " << 2*S.cols() << std::endl;
			throw(std::runtime_error("Error in symmetrization. Wrong number of states."));
		}
		return S_f.leftCols(S.cols());
	}
	
	std::vector<std::vector<int>> get_energy_levels(const Eigen::VectorXd& en)
	{
		std::vector<std::vector<int>> energy_levels;
		energy_levels.push_back({0});
		for (int i = 1; i < lat.n_sites(); ++i)
		{
			if (std::abs(en(i) - en(energy_levels.back()[0])) > param.epsilon)
				energy_levels.push_back(std::vector<int>());
			energy_levels.back().push_back(i);
		}
		for (int i = 0; i < energy_levels.size(); ++i)
			std::sort(energy_levels[i].begin(), energy_levels[i].end());
		return energy_levels;
	}
	
	matrix_t orthonormalize(const matrix_t& S)
	{
		matrix_t S_o = S, S_f = matrix_t::Zero(S.rows(), S.cols());

		std::cout << "start orthogonalize" << std::endl;
		
		for (int i = 0; i < S_o.cols(); ++i)
		{
			if (S_o.col(i).norm() > param.epsilon)
			{
				std::cout << "S_o norm = " << S_o.col(i).norm() << std::endl;
				S_o.col(i) /= S_o.col(i).norm();
			}
			else
				S_o.col(i) *= 0.;
		}

		for (int i = 0; i < S_o.cols(); ++i)
		{
			S_f.col(i) = S_o.col(i);
			for (int k = 0; k < i; ++k)
				S_f.col(i) -= S_f.col(k) * (S_f.col(k).adjoint() * S_o.col(i));
			
			std::cout << "S_f norm = " << S_f.col(i).norm() << std::endl;
			if (S_f.col(i).norm() > param.epsilon)
				S_f.col(i) /= S_f.col(i).norm();
		}
		
		std::cout << "end orthogonalize" << std::endl << std::endl;
		
		return S_f;
	}
	
	void split_quantum_numbers(std::vector<std::vector<int>>& energy_levels, const matrix_t& S, const matrix_t& pm)
	{
		for (int i = 0; i < energy_levels.size(); ++i)
		{
			std::vector<std::vector<int>> sub_levels;
			std::vector<numeric_t> quantum_numbers;
			sub_levels.push_back({energy_levels[i][0]});
			numeric_t q = S.col(energy_levels[i][0])
				.adjoint() * pm * S.col(energy_levels[i][0]);
			q = std::real(q);
			quantum_numbers.push_back(q);
			//std::cout << "level " << i << ", state " << 0
			//		<< ", q = " << q << std::endl;
			for (int j = 1; j < energy_levels[i].size(); ++j)
			{
				numeric_t q = S.col(energy_levels[i][j])
					.adjoint() * pm * S.col(energy_levels[i][j]);
				q = std::real(q);
				//std::cout << "level " << i << ", state " << j
				//	<< ", q = " << q << std::endl;
				int k;
				for (k = 0; k < quantum_numbers.size();)
					if (std::abs(quantum_numbers[k] - q) > param.epsilon)
						++k;
					else
						break;
				if (k == quantum_numbers.size())
				{
					quantum_numbers.push_back(q);
					sub_levels.push_back({energy_levels[i][j]});
				}
				else
					sub_levels[k].push_back(energy_levels[i][j]);
			}
			energy_levels.erase(energy_levels.begin() + i);
			for (int k = 0; k < quantum_numbers.size(); ++k)
				energy_levels.insert(energy_levels.begin()+i,
					sub_levels[sub_levels.size()-1-k]);
		}
		for (int i = 0; i < energy_levels.size(); ++i)
			std::sort(energy_levels[i].begin(), energy_levels[i].end());
	}
	
	/*
	matrix_t project_symmetry(const matrix_t& S, const std::vector<std::vector<int>>& energy_levels, const matrix_t& pm)
	{
		matrix_t S_f = matrix_t::Zero(lat.n_sites(), S.cols());
		
		for (int i = 0; i < energy_levels.size(); ++i)
		{
			int N = energy_levels[i].size();
			matrix_t projP(N, N);
			matrix_t S_proj = matrix_t::Zero(lat.n_sites(), N);
			for (int j = 0; j < N; ++j)
				for (int k = 0; k < N; ++k)
					projP(j, k) = S.col(energy_levels[i][j]).adjoint() * pm * S.col(energy_levels[i][k]);
				
			Eigen::ComplexEigenSolver<matrix_t> solver(projP);
			std::cout << "Projected eigenvalues: i = " << i << std::endl;
			for (int j = 0; j < N; ++j)
				std::cout << solver.eigenvalues()[j] << std::endl;
			for (int j = 0; j < N; ++j)
				for (int k = 0; k < N; ++k)
					S_proj.col(j) += solver.eigenvectors()(k, j) * S.col(energy_levels[i][k]);
			S_proj = orthonormalize(S_proj);
			for (int j = 0; j < N; ++j)
				S_f.col(energy_levels[i][j]) = S_proj.col(j);
		}
		return S_f;
	}
	
	matrix_t project_ph_symmetry(const matrix_t& S, const matrix_t& pm)
	{
		matrix_t S_f = matrix_t::Zero(lat.n_sites(), S.cols());
		
		int N = S.cols();
		matrix_t projP(N, N);
		matrix_t S_proj = matrix_t::Zero(lat.n_sites(), N);
		for (int j = 0; j < N; ++j)
			for (int k = 0; k < N; ++k)
				projP(j, k) = S.col(j).adjoint() * pm * S.col(k);
			
		Eigen::ComplexEigenSolver<matrix_t> solver(projP);
		std::cout << "PH Projected eigenvalues:" << std::endl;
		for (int j = 0; j < N; ++j)
			std::cout << solver.eigenvalues()[j] << std::endl;
		for (int j = 0; j < N; ++j)
			for (int k = 0; k < N; ++k)
			{
				S_proj.col(j) += solver.eigenvectors()(k, j) * S.col(k);
				if (std::imag(solver.eigenvectors()(k, j)) > param.epsilon)
					std::cout << "Imag value: " << solver.eigenvectors()(k, j) << std::endl;
			}
		S_proj = orthonormalize(S_proj);
		for (int j = 0; j < N; ++j)
			S_f.col(j) = S_proj.col(j);
		return S_f;
	}
	*/

	void print_representations(const matrix_t& S, const matrix_t& inv_pm, const matrix_t& sv_pm, const matrix_t& sh_pm,
		const matrix_t& rot60_pm, const matrix_t& rot120_pm, const matrix_t& ph_pm)
	{
		std::vector<matrix_t> rep(6, matrix_t::Zero(S.cols(), S.cols()));
		for (int i = 0; i < S.cols(); ++i)
			for (int j = 0; j < S.cols(); ++j)
			{
				rep[0](i, j) = S.col(i).adjoint() * inv_pm * S.col(j);
				rep[1](i, j) = S.col(i).adjoint() * sv_pm * S.col(j);
				rep[2](i, j) = S.col(i).adjoint() * sh_pm * S.col(j);
				rep[3](i, j) = S.col(i).adjoint() * rot60_pm * S.col(j);
				rep[4](i, j) = S.col(i).adjoint() * rot120_pm * S.col(j);
				rep[5](i, j) = S.col(i).adjoint() * ph_pm * S.col(j);
				
				for (int k = 0; k < 6; ++k)
					if (std::abs(rep[k](i, j)) < 1E-14)
						rep[k](i, j) = 0.;
			}
		std::cout << "rep inv_pm" << std::endl;
		for (int i = 0; i < rep[0].rows(); ++i)
		{
			for (int j = 0; j < rep[0].rows(); ++j)
				std::cout << rep[0](i, j) << " ";
			std::cout << std::endl;
		}
		std::cout << "rep sv_pm" << std::endl;
		for (int i = 0; i < rep[0].rows(); ++i)
		{
			for (int j = 0; j < rep[0].rows(); ++j)
				std::cout << rep[1](i, j) << " ";
			std::cout << std::endl;
		}
		std::cout << "rep sh_pm" << std::endl;
		for (int i = 0; i < rep[0].rows(); ++i)
		{
			for (int j = 0; j < rep[0].rows(); ++j)
				std::cout << rep[2](i, j) << " ";
			std::cout << std::endl;
		}
		std::cout << "rep rot60_pm" << std::endl;
		for (int i = 0; i < rep[0].rows(); ++i)
		{
			for (int j = 0; j < rep[0].rows(); ++j)
				std::cout << rep[3](i, j) << " ";
			std::cout << std::endl;
		}
		std::cout << "rep rot120_pm" << std::endl;
		for (int i = 0; i < rep[0].rows(); ++i)
		{
			for (int j = 0; j < rep[0].rows(); ++j)
				std::cout << rep[4](i, j) << " ";
			std::cout << std::endl;
		}
		std::cout << "rep ph_pm" << std::endl;
		for (int i = 0; i < rep[0].rows(); ++i)
		{
			for (int j = 0; j < rep[0].rows(); ++j)
				std::cout << rep[5](i, j) << " ";
			std::cout << std::endl;
		}
	}
	
	void print_energy_levels(const matrix_t& S, const Eigen::VectorXd& eigen_energies,
		const std::vector<std::vector<int>>& energy_levels,
		const matrix_t& inv_pm, const matrix_t& sv_pm, const matrix_t& sh_pm, const matrix_t& rot60_pm, const matrix_t& rot120_pm)
	{
		std::cout << "Single particle eigenvalues:" << std::endl;
		for (int k = 0; k < energy_levels.size(); ++k)
		{
			std::cout << "level " << k << ":" << std::endl;
			for (int j = 0; j < energy_levels[k].size(); ++j)
			{
				int i = energy_levels[k][j];
				std::cout << "E(" << i << ") = " << eigen_energies[i]
				<< ", P_inv = " << S.col(i).adjoint() * inv_pm * S.col(i)
				<< ", P_rot60 = " << S.col(i).adjoint() * rot60_pm * S.col(i)
				<< ", P_rot120 = " << S.conjugate().col(i).adjoint() * rot120_pm * S.conjugate().col(i)
				<< ", P_sv = " << S.conjugate().col(i).adjoint() * sv_pm * S.conjugate().col(i) 
				<< ", P_sh = " << S.conjugate().col(i).adjoint() * sh_pm * S.conjugate().col(i) << std::endl;
			}
			std::cout << "---" << std::endl;
		}
	}

	void trigger()
	{
		green_function::matrix_t K = green_function::matrix_t::Zero(lat.n_sites(), lat.n_sites());
		for (auto& a : lat.bonds("nearest neighbors"))
			K(a.first, a.second) = -param.t;
		for (auto& a : lat.bonds("t3_bonds"))
			K(a.first, a.second) = -param.tprime;
		gf.set_K_matrix(K);
		matrix_t P;
		
		Eigen::SelfAdjointEigenSolver<matrix_t> solver(K);
		matrix_t inv_pm = matrix_t::Zero(lat.n_sites(), lat.n_sites()),
			ph_pm = matrix_t::Zero(lat.n_sites(), lat.n_sites()),
			rot60_pm = matrix_t::Zero(lat.n_sites(), lat.n_sites()),
			rot120_pm = matrix_t::Zero(lat.n_sites(), lat.n_sites()),
			sv_pm = matrix_t::Zero(lat.n_sites(), lat.n_sites()),
			sh_pm = matrix_t::Zero(lat.n_sites(), lat.n_sites());
		for (int i = 0; i < lat.n_sites(); ++i)
		{
			inv_pm(i, lat.inverted_site(i)) = 1.;
			sv_pm(i, lat.reflected_v_site(i)) = 1.;
			sh_pm(i, lat.reflected_h_site(i)) = 1.;
			rot60_pm(i, lat.rotated_site(i, 60.)) = 1.;
			rot120_pm(i, lat.rotated_site(i, 120.)) = 1.;
			ph_pm(i, i) = lat.parity(i);
		}
		
		std::vector<numeric_t> total_quantum_numbers = {{1., 1., 1., 1., 1.}};
		std::vector<numeric_t> ph_2p_parity(4);
		std::vector<std::vector<int>> energy_levels = get_energy_levels(solver.eigenvalues());
		
		auto S_f = symmetrize_EV(solver.eigenvectors(), solver.eigenvalues(), inv_pm);
		
		if (lat.n_sites() % 3 != 0)
			P = S_f.leftCols(lat.n_sites()/2);
		else
		{
			for (int i = 0; i < lat.n_sites()/2-2; ++i)
			{
				total_quantum_numbers[0] *= (S_f.col(i).adjoint() * inv_pm * S_f.col(i)).trace();
				total_quantum_numbers[1] *= (S_f.col(i).adjoint() * sv_pm * S_f.col(i)).trace();
				total_quantum_numbers[2] *= (S_f.col(i).adjoint() * sh_pm * S_f.col(i)).trace();
				total_quantum_numbers[3] *= (S_f.col(i).adjoint() * rot60_pm * S_f.col(i)).trace();
				total_quantum_numbers[4] *= (S_f.col(i).adjoint() * rot120_pm * S_f.col(i)).trace();
			}
		
			matrix_t ph_1p_block = S_f.block(0, lat.n_sites()/2-2, lat.n_sites(), 4);
			Eigen::VectorXd ph_ev = Eigen::VectorXd::Zero(4);
			
			ph_1p_block = symmetrize_EV(ph_1p_block, ph_ev, ph_pm);
			//ph_1p_block = project_ph_symmetry(ph_1p_block, ph_pm);
			
			ph_1p_block = ph_symmetrize_EV(ph_1p_block, ph_pm, inv_pm);

			for (int i = 0; i < ph_1p_block.cols(); ++i)
				ph_2p_parity[i] = ph_1p_block.col(i).adjoint() * ph_pm * ph_1p_block.col(i);
			std::vector<matrix_t> ph_2p_block(4, matrix_t(lat.n_sites(), 2));
			
			//PH = -1
			ph_2p_block[0].col(0) = ph_1p_block.col(0);
			ph_2p_block[0].col(1) = ph_1p_block.col(3);
			ph_2p_block[1].col(0) = ph_1p_block.col(1);
			ph_2p_block[1].col(1) = ph_1p_block.col(2);
			//PH = 1
			ph_2p_block[2].col(0) = ph_1p_block.col(0);
			ph_2p_block[2].col(1) = ph_1p_block.col(1);
			ph_2p_block[3].col(0) = ph_1p_block.col(2);
			ph_2p_block[3].col(1) = ph_1p_block.col(3);
			
			ph_2p_parity = {-1., -1., 1., 1.};
			std::vector<std::vector<numeric_t>> e0_quantum_numbers = {4, std::vector<numeric_t>()};
			for (int i = 0; i < ph_2p_block.size(); ++i)
			{
				e0_quantum_numbers[i].push_back(ph_2p_block[i].col(0).adjoint() * inv_pm * ph_2p_block[i].col(0)
					* ph_2p_block[i].col(1).adjoint() * inv_pm * ph_2p_block[i].col(1));
				e0_quantum_numbers[i].push_back(ph_2p_block[i].col(0).adjoint() * sv_pm * ph_2p_block[i].col(0)
					* ph_2p_block[i].col(1).adjoint() * inv_pm * ph_2p_block[i].col(1));
				e0_quantum_numbers[i].push_back(ph_2p_block[i].col(0).adjoint() * sh_pm * ph_2p_block[i].col(0)
					* ph_2p_block[i].col(1).adjoint() * inv_pm * ph_2p_block[i].col(1));
				e0_quantum_numbers[i].push_back(ph_2p_block[i].col(0).adjoint() * rot60_pm * ph_2p_block[i].col(0)
					* ph_2p_block[i].col(1).adjoint() * inv_pm * ph_2p_block[i].col(1));
				e0_quantum_numbers[i].push_back(ph_2p_block[i].col(0).adjoint() * rot120_pm * ph_2p_block[i].col(0)
					* ph_2p_block[i].col(1).adjoint() * inv_pm * ph_2p_block[i].col(1));
			}
			
			/////////////////

			//auto energy_blocks = get_energy_blocks(solver.eigenvalues());
			//matrix_t ph_Np_block = S_f.block(0, 0, lat.n_sites(), lat.n_sites());
			//matrix_t ph_Np_block = solver.eigenvectors().block(0, 0, lat.n_sites(), lat.n_sites());
			//ph_Np_block = symmetrize_ph_blocks(ph_Np_block, energy_blocks, solver.eigenvalues(), ph_pm);
			//select_ph_states(ph_Np_block, S_f.block(0, 0, lat.n_sites(), lat.n_sites()), energy_blocks, H, ph_pm);

			/////////////////

			P.resize(lat.n_sites(), lat.n_sites() / 2);
			for (int i = 0; i < lat.n_sites()/2-2; ++i)
				P.col(i) = S_f.col(i);
			
			for (int i = 0; i < ph_2p_block.size(); ++i)
				if (std::abs(total_quantum_numbers[0] * e0_quantum_numbers[i][0] - param.inv_symmetry) < param.epsilon)
				{
					//std::cout << "Taken: i=" << i << std::endl;
					P.block(0, lat.n_sites()/2-2, lat.n_sites(), 2) = ph_2p_block[i];
					for (int j = 0; j < total_quantum_numbers.size(); ++j)
						total_quantum_numbers[j] *= e0_quantum_numbers[i][j];
					break;
				}
		}
		
		gf.set_trial_wf(P);
	}
	
	void init() {}
};

struct event_build
{
	Random& rng;
	parameters& param;
	lattice& lat;
	green_function& gf;

	void trigger()
	{
		unsigned Nv = 0.15 * (param.theta * lat.n_sites() * param.V);
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
			else if (param.static_obs[i] == "M4")
				add_wick(wick_static_M4{rng, param, lat}, param.static_obs[i]);
			else if (param.static_obs[i] == "epsilon")
				add_wick(wick_static_epsilon{rng, param, lat}, param.static_obs[i]);
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
	void add_vector_wick(T&& functor, const std::string& name)
	{
		vec_obs.push_back(vector_wick_static_base<matrix_t>(std::forward<T>(functor)));
		vec_names.push_back(name);
	}

	void trigger()
	{
		if (obs.size() == 0 and vec_obs.size() == 0)
			return;
		++param.static_measure_cnt;
		if (param.static_measure_cnt >= param.static_measure_interval)
		if (std::abs(gf.tau() - param.theta/2.) < param.measure_window/2.)
		{
			//std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
			gf.measure_static_observables(measure, names, obs, vec_names, vec_obs);
			//std::chrono::steady_clock::time_point t1 = std::chrono::steady_clock::now();
			//std::cout << "Time of static measurement: " << std::chrono::duration_cast<std::chrono::duration<float>>(t1 - t0).count() << std::endl;
			param.static_measure_cnt = 0;
		}
	}
	
	void init()
	{
		for (int i = 0; i < obs.size(); ++i)
			measure.add_observable(names[i], param.n_prebin);
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
			if (param.dyn_obs[i] == "epsilon")
				add_wick(wick_epsilon{rng, param, lat}, param.dyn_obs[i]);
			if (param.dyn_obs[i] == "kekule_s")
				add_wick(wick_kekule_s{rng, param, lat}, param.dyn_obs[i]);
			if (param.dyn_obs[i] == "kekule_as")
				add_wick(wick_kekule_as{rng, param, lat}, param.dyn_obs[i]);
			if (param.dyn_obs[i] == "kekule_K")
				add_wick(wick_kekule_K{rng, param, lat}, param.dyn_obs[i]);
			if (param.dyn_obs[i] == "gamma_mod")
				add_wick(wick_gamma_mod{rng, param, lat}, param.dyn_obs[i]);
			if (param.dyn_obs[i] == "sp")
				add_wick(wick_sp{rng, param, lat}, param.dyn_obs[i]);
			if (param.dyn_obs[i] == "tp")
				add_wick(wick_tp{rng, param, lat}, param.dyn_obs[i]);
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
		
		if (std::abs(gf.tau() - (param.theta/2. + param.dyn_tau_max/2 + param.block_size/2.)) < 1E-6
			or std::abs(gf.tau() - (param.theta/2. - param.dyn_tau_max/2 + param.block_size/2.)) < 1E-6)
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
