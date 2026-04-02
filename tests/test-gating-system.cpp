#define BOOST_TEST_MODULE GatingSystemTests
#include <boost/test/included/unit_test.hpp>

#include "HHConfig.hpp"
#include "gating-system.hpp"
#include "state.hpp"
#include "model.hpp"
#include "time-integrator.hpp"
#include <cmath>
#include <vector>

const int N = 3; // number of compartments for testing
const HHConfig config;
const double V_rest = Compartment::compute_V_rest(config).value();

// Helper 1: Build a State forced away from equilibrium
State make_unsteady_state() {
	State state(4 * N);

	for (int i = 0; i < N; ++i) {
		state[i] = V_rest;     // V
		state[N + i] = -1;     // m
		state[2 * N + i] = -1; // h
		state[3 * N + i] = -1; // n
	}
	return state;
}

// At steady state, rhs() must return ~0 for all gating variables
BOOST_AUTO_TEST_CASE(gating_rhs_zero_at_steady_state) {
	State state = HHModel::make_steady_state(N, config);
	GatingSystem sys(N);

	std::vector<double> dydt(3 * N, 0.0);
	sys.rhs(state, dydt);

	for (int i = 0; i < 3 * N; ++i) {
		BOOST_CHECK_SMALL(dydt[i], 1e-10);
	}
}

// for small time step, forward euler should cause convergence to steady state
BOOST_AUTO_TEST_CASE(forward_euler_converges_to_steady_state) {
	State state = make_unsteady_state();
	GatingSystem sys(N);
	ForwardEuler fe;

	std::vector<double> ss(3);
	Compartment::get_steady_state_gating_variables(V_rest, ss);
	double m_inf = ss[0], h_inf = ss[1], n_inf = ss[2];

	double dt = 0.01; // Small enough to remain stable
	double t = 0.0;
	for (int step = 0; step < 20000; ++step) {
		fe.step(sys, state, dt);
		t += dt;
	}

	for (int i = 0; i < N; ++i) {
		BOOST_CHECK_CLOSE(state[N + i], m_inf, 0.1);
		BOOST_CHECK_CLOSE(state[2 * N + i], h_inf, 0.1);
		BOOST_CHECK_CLOSE(state[3 * N + i], n_inf, 0.1);
	}
}

// stiff problem, large time step should cause instability for forward euler
BOOST_AUTO_TEST_CASE(forward_euler_instability) {
	State state = make_unsteady_state();
	GatingSystem sys(N);
	ForwardEuler fe;

	std::vector<double> ss(3);
	Compartment::get_steady_state_gating_variables(V_rest, ss);
	double m_inf = ss[0], h_inf = ss[1], n_inf = ss[2];

	double dt = 1.0; // much bigger than the stability limit for Forward Euler on this system, should diverge
	double t = 0.0;
	for (int step = 0; step < 100000; ++step) {
		fe.step(sys, state, dt);
		t += dt;
	}

	for (int i = 0; i < N; ++i) {
		double m = state[N + i];
		double h = state[2 * N + i];
		double n = state[3 * N + i];

		BOOST_CHECK(!std::isfinite(m) || std::abs(m - m_inf) > 0.1);

		// Note: h and n do not diverge since the time constant here is slower than m, but m alone is sufficient to demonstrate instability of the method, so we do not require h and n to diverge in this test.
	}
}

// Implicit Euler should remain stable even with large time step, but should still converge to steady state
BOOST_AUTO_TEST_CASE(implicit_euler_converges_to_steady_state) {
	State state = make_unsteady_state();
	GatingSystem sys(N);
	ImplicitEuler ie;

	std::vector<double> ss(3);
	Compartment::get_steady_state_gating_variables(V_rest, ss);
	double m_inf = ss[0], h_inf = ss[1], n_inf = ss[2];

	double dt = 1.0;
	double t = 0.0;
	for (int step = 0; step < 500; ++step) {
		ie.step(sys, state, dt);
		t += dt;
	}

	for (int i = 0; i < N; ++i) {
		BOOST_CHECK_CLOSE(state[N + i], m_inf, 0.1);
		BOOST_CHECK_CLOSE(state[2 * N + i], h_inf, 0.1);
		BOOST_CHECK_CLOSE(state[3 * N + i], n_inf, 0.1);
	}
}

// Implicit Euler should keep gates in [0, 1] even with absurdly large dt
BOOST_AUTO_TEST_CASE(implicit_euler_gates_bounded) {
	State state = make_unsteady_state();
	GatingSystem sys(N);
	ImplicitEuler ie;

	double dt = 100.0;
	double t = 0.0;
	for (int step = 0; step < 10; ++step) {
		ie.step(sys, state, dt);
		t += dt;
	}

	for (int i = 0; i < N; ++i) {
		BOOST_CHECK(state[N + i] >= 0.0 && state[N + i] <= 1.0);
		BOOST_CHECK(state[2 * N + i] >= 0.0 && state[2 * N + i] <= 1.0);
		BOOST_CHECK(state[3 * N + i] >= 0.0 && state[3 * N + i] <= 1.0);
	}
}

// NOTE: this is not required by the assignment, but it is commonly used in the literature.
BOOST_AUTO_TEST_CASE(exponential_euler_converges) {
	State state = make_unsteady_state();
	GatingSystem sys(N);
	ExponentialEuler ee;

	std::vector<double> ss(3);
	Compartment::get_steady_state_gating_variables(V_rest, ss);
	double m_inf = ss[0], h_inf = ss[1], n_inf = ss[2];

	double dt = 1.0; // this is the timestep at which Forward Euler becomes unstable
	double t = 0.0;
	for (int step = 0; step < 2000; ++step) {
		ee.step(sys, state, dt);
		t += dt;
	}

	for (int i = 0; i < N; ++i) {
		BOOST_CHECK_CLOSE(state[N + i], m_inf, 0.1);
		BOOST_CHECK_CLOSE(state[2 * N + i], h_inf, 0.1);
		BOOST_CHECK_CLOSE(state[3 * N + i], n_inf, 0.1);
	}
}