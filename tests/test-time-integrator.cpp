#define BOOST_TEST_MODULE TimeIntegratorTests
#include <boost/test/included/unit_test.hpp>

#include "state.hpp"
#include "linear-test-equation.hpp"
#include "time-integrator.hpp"
#include <cmath>
#include <vector>

class MockDiagonalSystem : public DiagonalSystem {
public:
    void rhs(const State&, std::vector<double>&) override {}
    int get_size() const override { return 1; }
    int get_offset() const override { return 0; }

    void assemble_diagonal_system(const State&, double, std::vector<double>&, std::vector<double>&) override {}

    void accept(TimeIntegrator &integrator, State &state, double dt) override {
        integrator.visit(*this, state, dt);
    }
};

class MockTridiagonalSystem : public TridiagonalSystem {
public:
    void rhs(const State&, std::vector<double>&) override {}
    int get_size() const override { return 1; }
    int get_offset() const override { return 0; }

    void assemble_tridiagonal_system(const State&, double, std::vector<double>&,  std::vector<double>&, std::vector<double>&, std::vector<double>&) override {}

    void accept(TimeIntegrator &integrator, State &state, double dt) override {
        integrator.visit(*this, state, dt);
    }
};

BOOST_AUTO_TEST_CASE(visit_methods_must_be_implemented) {
    MockDiagonalSystem diag_sys;
    MockTridiagonalSystem tri_sys;
    State state;

    ForwardEuler fe;
    ExponentialEuler ee;

    BOOST_CHECK_THROW(fe.visit(diag_sys, state, 0.1), std::runtime_error);
    BOOST_CHECK_THROW(fe.visit(tri_sys, state, 0.1), std::runtime_error);
    BOOST_CHECK_THROW(ee.visit(tri_sys, state, 0.1), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(forward_euler_single_variable) {
	State state;
	state.y = {1.0};

	LinearTestEquation sys(-1.0, 1, 0);
	ForwardEuler fe;

	double dt = 0.001;
	double t = 0.0;
	int steps = 1000;

	for (int i = 0; i < steps; ++i) {
		fe.step(sys, state, dt);
		t += dt;
	}

	double exact = std::exp(-1.0);
	BOOST_CHECK_CLOSE(state[0], exact, 0.1);
}

BOOST_AUTO_TEST_CASE(implicit_euler_single_variable) {
	State state;
	state.y = {1.0};

	LinearTestEquation sys(-1.0, 1, 0);
	ImplicitEuler ie;

	double dt = 0.001;
	double t = 0.0;
	int steps = 1000;

	for (int i = 0; i < steps; ++i) {
		ie.step(sys, state, dt);
		t += dt;
	}

	double exact = std::exp(-1.0);
	BOOST_CHECK_CLOSE(state[0], exact, 0.1);
}

BOOST_AUTO_TEST_CASE(forward_euler_with_offset) {
	State state;
	state.y = {0.0, 0.0, 2.0};

	LinearTestEquation sys(-1.0, 1, 2);
	ForwardEuler fe;

	double dt = 0.001;
	double t = 0.0;
	int steps = 1000;

	for (int i = 0; i < steps; ++i) {
		fe.step(sys, state, dt);
		t += dt;
	}

	double exact = 2.0 * std::exp(-1.0);
	BOOST_CHECK_CLOSE(state[2], exact, 0.1);
	BOOST_CHECK_SMALL(state[0], 1e-12);
	BOOST_CHECK_SMALL(state[1], 1e-12);
}

BOOST_AUTO_TEST_CASE(forward_euler_multiple_variables) {
	State state;
	state.y = {1.0, 2.0, 3.0};

	LinearTestEquation sys(-1.0, 3, 0);
	ForwardEuler fe;

	double dt = 0.001;
	double t = 0.0;
	int steps = 1000;

	for (int i = 0; i < steps; ++i) {
		fe.step(sys, state, dt);
		t += dt;
	}

	double factor = std::exp(-1.0);
	BOOST_CHECK_CLOSE(state[0], 1.0 * factor, 0.1);
	BOOST_CHECK_CLOSE(state[1], 2.0 * factor, 0.1);
	BOOST_CHECK_CLOSE(state[2], 3.0 * factor, 0.1);
}

BOOST_AUTO_TEST_CASE(implicit_euler_multiple_variables) {
	State state;
	state.y = {1.0, 2.0, 3.0};

	LinearTestEquation sys(-1.0, 3, 0);
	ImplicitEuler ie;

	double dt = 0.001;
	double t = 0.0;
	int steps = 1000;

	for (int i = 0; i < steps; ++i) {
		ie.step(sys, state, dt);
		t += dt;
	}

	double factor = std::exp(-1.0);
	BOOST_CHECK_CLOSE(state[0], 1.0 * factor, 0.1);
	BOOST_CHECK_CLOSE(state[1], 2.0 * factor, 0.1);
	BOOST_CHECK_CLOSE(state[2], 3.0 * factor, 0.1);
}

BOOST_AUTO_TEST_CASE(forward_euler_unstable) {
	State state;
	state.y = {1.0};

	LinearTestEquation sys(-1.0, 1, 0);
	ForwardEuler fe;

	double dt = 2.5;
	double t = 0.0;
	int steps = 50;

	for (int i = 0; i < steps; ++i) {
		fe.step(sys, state, dt);
		t += dt;
	}

	BOOST_CHECK(std::abs(state[0]) > 1e6); // verify that solution has blown up
}

BOOST_AUTO_TEST_CASE(implicit_euler_stays_stable) {
	State state;
	state.y = {1.0};

	LinearTestEquation sys(-1.0, 1, 0);
	ImplicitEuler ie;

	double dt = 2.5;
	double t = 0.0;
	int steps = 50;

	for (int i = 0; i < steps; ++i) {
		ie.step(sys, state, dt);
		t += dt;
	}

	BOOST_CHECK(std::abs(state[0]) < 1.0); // verify that solution has decayed, not blown up
}
