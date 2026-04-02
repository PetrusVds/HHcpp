#ifndef ROOT_FINDER_HPP
#define ROOT_FINDER_HPP

#include "HHConfig.hpp"
#include <cassert>
#include <cmath>
#include <optional>
#include <functional>
#include <stdexcept>
#include <vector>

namespace RootFinder {

	// Finds the root via Newton-Raphson iteration. See source: https://en.wikipedia.org/wiki/Newton%27s_method
	std::optional<double> newton_raphson(const std::function<double(double)> &f, double x0, double tol = 1e-10, int max_iter = 100, double dh = 1e-6);

} // namespace RootFinder

#endif // ROOT_FINDER_HPP