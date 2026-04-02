#ifndef SPARSE_SYSTEM_SOLVER_HPP
#define SPARSE_SYSTEM_SOLVER_HPP

#include "HHConfig.hpp"
#include <cassert>
#include <cmath>
#include <optional>
#include <functional>
#include <stdexcept>
#include <vector>

namespace SparseSystemSolver {
	/*
	* Solves a tridiagonal system of equations using the Thomas algorithm.
	*
	* Parameters:
	* - a: Sub-diagonal coefficients (length n, with a[0] unused)
	* - b: Main diagonal coefficients (length n)
	* - c: Super-diagonal coefficients (length n, with c[n-1] unused)
	* - d: Right-hand side vector (length n)
	* - x: Output vector containing the solution starting at x_offset
	* - x_offset: Index into x where the solution should be written
	*
	* Method follows the description found on wikipedia:
	* https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
	*/
	void thomas(const std::vector<double> &a, const std::vector<double> &b, std::vector<double> c, std::vector<double> d, std::vector<double> &x, int x_offset = 0);
}; // namespace SparseSystemSolver

#endif // SPARSE_SYSTEM_SOLVER_HPP