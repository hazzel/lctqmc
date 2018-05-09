#pragma once

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <vector>


struct vertex
{
	double tau;
	unsigned int x;
	unsigned int y;
};

class green_function
{
	public:
		using numeric_t = double;
		using matrix_t = Eigen::Matrix<numeric_t, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;

		green_function() {}

		matrix_t propagator(unsigned int tau_n, unsigned int tau_m)
		{
			dmatrix_t b = matrix_t::Identity(n_matrix_size, n_matrix_size);
			for (unsigned int n = tau_n; n > tau_m; --n)
				b *= K;

		}

		void wrap() {}
		void stabilize() {}
		void measure() {}
	private:
		unsigned int n_matrix_size;
		matrix_t g_tau;
		matrix_t K;
		std::vector<matrix_t> storage;
};
