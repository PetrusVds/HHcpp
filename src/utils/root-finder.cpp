#include "utils/root-finder.hpp"

std::optional<double> RootFinder::newton_raphson(const std::function<double(double)> &f, double x0, double tol, int max_iter, double dh) {
    double x = x0;

    for (int iter = 0; iter < max_iter; ++iter) {
        double fx = f(x);
        if (std::abs(fx) < tol) {
            return x;
        }

        double dfdx = (f(x + dh) - f(x - dh)) / (2.0 * dh);
        if (std::abs(dfdx) < 1e-15) {
            return std::nullopt; // failed
        }

        double step = fx / dfdx;
        x -= step;

        if (std::abs(step) < tol) {
            return x;
        }
    }

    return std::nullopt; // failed
}