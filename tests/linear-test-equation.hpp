#ifndef LINEAR_TEST_EQUATION_HPP
#define LINEAR_TEST_EQUATION_HPP

#include "state.hpp"
#include "system.hpp"
#include "time-integrator.hpp"
#include <cmath>
#include <vector>

// dy/dt = lambda * y, exact solution: y(t) = y0 * exp(lambda * t)
class LinearTestEquation : public DiagonalSystem {
public:
	LinearTestEquation(double lambda, int size, int offset) : lambda(lambda), n(size), off(offset) {}

	void rhs(const State &state, std::vector<double> &dydt) override {
		for (int i = 0; i < n; ++i) {
			dydt[i] = lambda * state[off + i];
		}
	}

	void assemble_diagonal_system(const State &state, double dt, std::vector<double> &diag, std::vector<double> &rhs) override {
		for (int i = 0; i < n; ++i) {
			diag[i] = 1.0 - lambda * dt;
			rhs[i] = state[off + i];
		}
	}

	inline void accept(TimeIntegrator &integrator, State &state, double dt) override {
		integrator.visit(*this, state, dt);
	}

	int get_size() const override { return n; }
	int get_offset() const override { return off; }

private:
	double lambda;
	int n;
	int off;
};

#endif // LINEAR_TEST_EQUATION_HPP