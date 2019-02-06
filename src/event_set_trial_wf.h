#pragma once
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
	
	matrix_t ph_symmetrize_EV(const matrix_t& S, const matrix_t& pm)
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
	
	matrix_t orthogonalize(const matrix_t& S)
	{
		Eigen::ColPivHouseholderQR<matrix_t> qr_solver(S);
		return qr_solver.matrixQ();
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
			//std::cout << "Projected eigenvalues: i = " << i << std::endl;
			//for (int j = 0; j < N; ++j)
			//	std::cout << solver.eigenvalues()[j] << std::endl;
			for (int j = 0; j < N; ++j)
				for (int k = 0; k < N; ++k)
					S_proj.col(j) += std::real(solver.eigenvectors()(k, j)) * S.col(energy_levels[i][k]);
			S_proj = orthogonalize(S_proj);
			for (int j = 0; j < N; ++j)
				S_f.col(energy_levels[i][j]) = S_proj.col(j);
		}
		return S_f;
	}

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
	
	matrix_t set_K_matrix()
	{
		matrix_t K = matrix_t::Zero(lat.n_sites(), lat.n_sites());
		if (param.geometry == "honeycomb")
		{
			for (auto& a : lat.bonds("nearest neighbors"))
				K(a.first, a.second) += -param.t;
			for (auto& a : lat.bonds("t3_bonds"))
				K(a.first, a.second) += -param.tprime;
		}
		else if (param.geometry == "square")
		{
			for (auto& a : lat.bonds("nearest neighbors"))
				K(a.first, a.second) += -param.t;
		}
		gf.set_K_matrix(K);
		return K;
	}
	
	matrix_t set_twf_matrix()
	{
		matrix_t tw = matrix_t::Zero(lat.n_sites(), lat.n_sites());
		if (param.geometry == "honeycomb")
		{
			if (param.trial_wave_function == "t_only")
			{
				for (auto& a : lat.bonds("nearest neighbors"))
					tw(a.first, a.second) += -param.t;

				/*
				std::complex<double> im(0., 1.);
				const double pi = std::atan(1.0)*4;
				double theta_x = pi, theta_y = pi;
				for (auto& a : lat.bonds("edge_x"))
				{
					tw(a.first, a.second) = -param.t * std::exp(im * theta_x);
					tw(a.second, a.first) = -param.t * std::exp(-im * theta_x);
				}
				for (auto& a : lat.bonds("edge_y"))
				{
					tw(a.first, a.second) = -param.t * std::exp(im * theta_y);
					tw(a.second, a.first) = -param.t * std::exp(-im * theta_y);
				}
				*/
				/*
				for (auto& a : lat.bonds("edge_x"))
				{
					tw(a.first, a.second) = +param.t;
					tw(a.second, a.first) = +param.t;
				}
				for (auto& a : lat.bonds("edge_y"))
				{
					tw(a.first, a.second) = +param.t;
					tw(a.second, a.first) = +param.t;
				}
				*/
			}
			else
			{
				for (auto& a : lat.bonds("nearest neighbors"))
					tw(a.first, a.second) += -param.t;
				for (auto& a : lat.bonds("t3_bonds"))
					tw(a.first, a.second) += -param.tprime;
			}
		}
		else if (param.geometry == "square")
		{
			for (auto& a : lat.bonds("nearest neighbors"))
				tw(a.first, a.second) += -param.t;
			
			int ns = lat.n_sites();
			for (int i = 0; i < lat.Lx; ++i)
				for (int j = 0; j < lat.Ly; ++j)
				{
					int n = j * lat.Lx + i;
					/*
					tw(n, (n+lat.Lx)%ns) += -param.t * (2*(i%2)-1);
					tw((n+lat.Lx)%ns, n) += -param.t * (2*(i%2)-1);
					
					if (i == lat.Lx - 1)
					{
						tw(n, (n - lat.Lx + 1)%ns) += -param.t;
						tw((n - lat.Lx + 1)%ns, n) += -param.t;
					}
					else
					{
						tw(n, (n + 1)%ns) += -param.t;
						tw((n + 1)%ns, n) += -param.t;
					}
					*/
					
					std::complex<double> im(0., 1.);
					const double pi = std::atan(1.0)*4;
					tw(n, (n+lat.Lx)%ns) += -param.t * std::exp(im * pi/4.);
					tw((n+lat.Lx)%ns, n) += -param.t * std::exp(-im * pi/4.);
					
					if (i == lat.Lx - 1)
					{
						tw(n, (n - lat.Lx + 1)%ns) += -param.t * std::exp(-im * pi/4.);
						tw((n - lat.Lx + 1)%ns, n) += -param.t * std::exp(im * pi/4.);
					}
					else
					{
						tw(n, (n + 1)%ns) += -param.t * std::exp(-im * pi/4.);
						tw((n + 1)%ns, n) += -param.t * std::exp(im * pi/4.);
					}
				}
		}
		return tw;
	}

	void trigger()
	{
		matrix_t K = set_K_matrix();
		if (!param.projective)
			return;
		matrix_t tw = set_twf_matrix();
		
		matrix_t P;
		Eigen::SelfAdjointEigenSolver<matrix_t> solver(tw);
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
			if (param.geometry == "honeycomb")
			{
				rot60_pm(i, lat.rotated_site(i, 60.)) = 1.;
				rot120_pm(i, lat.rotated_site(i, 120.)) = 1.;
			}
			ph_pm(i, i) = lat.parity(i);
		}
		
		if (param.geometry == "square")
			std::cout << solver.eigenvalues() << std::endl;
		
		if (param.geometry == "honeycomb" && lat.n_sites() % 3 != 0)
		{
			std::vector<std::vector<int>> energy_levels = get_energy_levels(solver.eigenvalues());
			auto S_f = solver.eigenvectors();
			S_f = project_symmetry(S_f, energy_levels, inv_pm);
			split_quantum_numbers(energy_levels, S_f, inv_pm);
			//print_energy_levels(S_f, solver.eigenvalues(), energy_levels, inv_pm, sv_pm, sh_pm, rot60_pm, rot120_pm);
			
			S_f = project_symmetry(S_f, energy_levels, sv_pm);
			split_quantum_numbers(energy_levels, S_f, sv_pm);
			//print_energy_levels(S_f, solver.eigenvalues(), energy_levels, inv_pm, sv_pm, sh_pm, rot60_pm, rot120_pm);
		
			S_f = project_symmetry(S_f, energy_levels, sh_pm);
			split_quantum_numbers(energy_levels, S_f, sh_pm);
			//print_energy_levels(S_f, solver.eigenvalues(), energy_levels, inv_pm, sv_pm, sh_pm, rot60_pm, rot120_pm);
			
			S_f = project_symmetry(S_f, energy_levels, rot60_pm);
			split_quantum_numbers(energy_levels, S_f, rot60_pm);
			//print_energy_levels(S_f, solver.eigenvalues(), energy_levels, inv_pm, sv_pm, sh_pm, rot60_pm, rot120_pm);
			
			P = S_f.leftCols(lat.n_sites()/2);
			
			//P = S_f.rightCols(lat.n_sites()/2);
			//P.col(P.cols()-1) = S_f.col(0);
			
			//double qn_P{1.};
			//for (int i = 0; i < P.cols(); ++i)
			//	qn_P *= (P.col(i).adjoint() * inv_pm * P.col(i)).trace();
			//std::cout << "qn_P = " << qn_P << std::endl;
			
			//P = solver.eigenvectors().leftCols(lat.n_sites()/2);
		}
		else
		{
			std::vector<numeric_t> total_quantum_numbers = {{1., 1., 1., 1., 1.}};
			std::vector<numeric_t> ph_2p_parity(4);
			std::vector<std::vector<int>> energy_levels = get_energy_levels(solver.eigenvalues());
			
			//auto S_f = symmetrize_EV(solver.eigenvectors(), solver.eigenvalues(), inv_pm);
			
			auto S_f = solver.eigenvectors();
			S_f = project_symmetry(S_f, energy_levels, inv_pm);
			split_quantum_numbers(energy_levels, S_f, inv_pm);
			
			S_f = project_symmetry(S_f, energy_levels, sv_pm);
			split_quantum_numbers(energy_levels, S_f, sv_pm);
		
			S_f = project_symmetry(S_f, energy_levels, sh_pm);
			split_quantum_numbers(energy_levels, S_f, sh_pm);
			if (param.geometry == "square")
				print_energy_levels(S_f, solver.eigenvalues(), energy_levels, inv_pm, sv_pm, sh_pm, rot60_pm, rot120_pm);
			
			if (param.geometry == "honeycomb")
			{
				S_f = project_symmetry(S_f, energy_levels, rot60_pm);
				split_quantum_numbers(energy_levels, S_f, rot60_pm);
			}
		
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
			ph_1p_block = ph_symmetrize_EV(ph_1p_block, ph_pm);

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
					* ph_2p_block[i].col(1).adjoint() * sv_pm * ph_2p_block[i].col(1));
				e0_quantum_numbers[i].push_back(ph_2p_block[i].col(0).adjoint() * sh_pm * ph_2p_block[i].col(0)
					* ph_2p_block[i].col(1).adjoint() * sh_pm * ph_2p_block[i].col(1));
				e0_quantum_numbers[i].push_back(ph_2p_block[i].col(0).adjoint() * rot60_pm * ph_2p_block[i].col(0)
					* ph_2p_block[i].col(1).adjoint() * rot60_pm * ph_2p_block[i].col(1));
				e0_quantum_numbers[i].push_back(ph_2p_block[i].col(0).adjoint() * rot120_pm * ph_2p_block[i].col(0)
					* ph_2p_block[i].col(1).adjoint() * rot120_pm * ph_2p_block[i].col(1));
			}

			P.resize(lat.n_sites(), lat.n_sites() / 2);
			for (int i = 0; i < lat.n_sites()/2-2; ++i)
				P.col(i) = S_f.col(i);
			
			//for (int i = 0; i < ph_2p_block.size(); ++i)
			//{
			//	std::cout << i << " " << total_quantum_numbers[0]*e0_quantum_numbers[i][0] << " " << total_quantum_numbers[1]*e0_quantum_numbers[i][1]
			//		<< " " << total_quantum_numbers[2]*e0_quantum_numbers[i][2] << std::endl;
			//}
			
			for (int i = 0; i < ph_2p_block.size(); ++i)
				if (std::abs(total_quantum_numbers[0] * e0_quantum_numbers[i][0] - param.inv_symmetry) < param.epsilon)
				{
					//std::cout << "Taken: i=" << i << std::endl;
					P.block(0, lat.n_sites()/2-2, lat.n_sites(), 2) = ph_2p_block[i];
					for (int j = 0; j < total_quantum_numbers.size(); ++j)
						total_quantum_numbers[j] *= e0_quantum_numbers[i][j];
					break;
				}
			//print_representations(P, inv_pm, sv_pm, sh_pm, rot60_pm, rot120_pm, ph_pm);
			
			if (std::abs(param.inv_symmetry - total_quantum_numbers[0]) > param.epsilon)
			{
				std::cout << "Error! Wrong parity of trial wave function." << std::endl;
				throw(std::runtime_error("Wrong parity in trial wave function."));
			}
		}
		gf.set_trial_wf(P);
	}
	
	void init() {}
};
