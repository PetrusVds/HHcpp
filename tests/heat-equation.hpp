#ifndef HEAT_EQUATION_HPP
#define HEAT_EQUATION_HPP

#include "simulator.hpp"
#include "state.hpp"
#include "model.hpp"
#include "system.hpp"
#include "time-integrator.hpp"

// 1D heat equation on [0, L] with Dirichlet BCs u(0,t) = u(L,t) = 0
class HeatSystem : public TridiagonalSystem {
public:
	HeatSystem(int N, double L) : N(N), dx(L / (N + 1)) {}

	void rhs(const State &state, std::vector<double> &dydt) override {
		for (int i = 0; i < N; ++i) {
			double u_left = (i > 0) ? state[off + i - 1] : 0.0;
			double u_right = (i < N - 1) ? state[off + i + 1] : 0.0;
			double u_center = state[off + i];
			dydt[i] = factor * (u_left - 2.0 * u_center + u_right);
		}
	}

	void assemble_tridiagonal_system(const State &state, double dt, std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &d) override {
		double mu = -dt * factor; // mu = dt / dx^2
		for (int i = 0; i < N; ++i) {
			if (i > 0) a[i] = mu;
			b[i] = 1.0 - 2.0 * mu;
			if (i < N - 1) c[i] = mu;
			d[i] = state[off + i];
		}
	}

	inline void accept(TimeIntegrator &integrator, State &state, double dt) override {
		integrator.visit(*this, state, dt);
	}

	int get_size() const override { return N; }
	int get_offset() const override { return off; }

private:
	int N;
	int off = 0;
	double dx;
	double factor = 1 / (dx * dx);
};

// Wrapper for the Heat Equation model
class HeatModel : public Model {
public:
	HeatModel(int N, double L, std::shared_ptr<TimeIntegrator> integrator = std::make_shared<ImplicitEuler>()) : Model(N, make_init_state(N, L)) {		
		sim.add_subsystem(std::make_shared<HeatSystem>(N, L), integrator);
	}

	std::vector<std::vector<double>> temperatures() const { return slice(0); }

	static State make_init_state(int N, double L) {
		State state;
		state.y.resize(N);
		double dx = L / (N + 1);
		for (int i = 0; i < N; ++i) {
			double x = (i + 1) * dx;
			state.y[i] = std::sin(M_PI * x / L);
		}
		return state;
	}

    // exact solution for the heat equation with initial condition u(x,0) = sin(pi*x/L) and Dirichlet BCs u(0,t) = u(L,t) = 0
    static double exact_solution(double x, double t, double L) {
        return std::sin(M_PI * x / L) * std::exp(-std::pow(M_PI / L, 2) * t);
    }
};

#endif // HEAT_EQUATION_HPP