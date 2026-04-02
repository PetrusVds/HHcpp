#include "utils/sparse-system-solver.hpp"

void SparseSystemSolver::thomas(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> c, std::vector<double> d, std::vector<double> &x, int x_offset) {
		int n = b.size();
		double t;

		// forward sweep
		c[0] /= b[0];
		d[0] /= b[0];

		for (int i = 1; i < n; ++i) {
			t = 1.0 / (b[i] - a[i] * c[i - 1]);

			if (i < n - 1) {
				c[i] = c[i] * t;
			}
			d[i] = (d[i] - a[i] * d[i - 1]) * t;
		}

		// back substitution
		x[x_offset + n - 1] = d[n - 1];
		for (int i = n - 2; i >= 0; --i) {
			x[x_offset + i] = d[i] - c[i] * x[x_offset + i + 1];
		}
	}