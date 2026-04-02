#define BOOST_TEST_MODULE HeatEquationTests
#include <boost/test/included/unit_test.hpp>

#include "heat-equation.hpp"
#include <cmath>
#include <memory>
#include <vector>

// Constants for the tests
static constexpr int N_test = 50;
static constexpr double L_test = 1.0;

// Forward Euler should converge when mu <= 1/2 (CFL condition)
BOOST_AUTO_TEST_CASE(forward_euler_heat_converges) {
	HeatModel model(N_test, L_test, std::make_shared<ForwardEuler>());

	double dx = L_test / (N_test + 1);
	double dt = 0.4 * dx * dx; // stay inside the CFL limit
	double t_end = 1.0;

	model.run(t_end, dt);

	auto temps = model.temperatures();
	const auto &last_temps = temps.back();
	double t_final = model.times().back();

	for (int i = 0; i < N_test; ++i) {
		double x = (i + 1) * dx;
		double exact = HeatModel::exact_solution(x, t_final, L_test);
		BOOST_CHECK_CLOSE(last_temps[i], exact, 1.0);
	}
}

// Forward Euler should blow up when mu > 1/2
BOOST_AUTO_TEST_CASE(forward_euler_heat_unstable_large_dt) {
	HeatModel model(N_test, L_test, std::make_shared<ForwardEuler>());

	double dx = L_test / (N_test + 1);
	double dt = 5.0 * dx * dx; // well above CFL limit
	double t_end = 50 * dt;

	model.run(t_end, dt);

	auto temps = model.temperatures();
	const auto &last_temps = temps.back();
	double max_val = 0.0;
	for (int i = 0; i < N_test; ++i) {
		max_val = std::max(max_val, std::abs(last_temps[i]));
	}
	BOOST_CHECK(max_val > 1e6); // verify that solution has blown up
}

// Implicit Euler should converge for any timestep (tests successful initialization of default argument)
BOOST_AUTO_TEST_CASE(implicit_euler_heat_converges) {
	HeatModel model(N_test, L_test); // Defaults to Implicit Euler

	double dx = L_test / (N_test + 1);
	double dt = 0.01;
	double t_end = 100.0;

	model.run(t_end, dt);

	auto temps = model.temperatures();
	const auto &last_temps = temps.back();
	double t_final = model.times().back();
	for (int i = 0; i < N_test; ++i) {
		double x = (i + 1) * dx;
		double exact = HeatModel::exact_solution(x, t_final, L_test);
		BOOST_CHECK_CLOSE(last_temps[i], exact, 1.0);
	}
}

// Implicit Euler is unconditionally stable, so it should converge even with large timesteps
BOOST_AUTO_TEST_CASE(implicit_euler_heat_stable_large_dt) {
	HeatModel model(N_test, L_test);

	double dt = 10.0;
	double t_end = 20 * dt;

	model.run(t_end, dt);

	auto temps = model.temperatures();
	const auto &last_temps = temps.back();

	for (int i = 0; i < N_test; ++i) {
		BOOST_CHECK(std::abs(last_temps[i]) < 2.0); // initial peak was 1.0
	}
}

// Full end-to-end test of the framework (Simulator trajectory storage
// capability verification)
BOOST_AUTO_TEST_CASE(heat_equation_end_to_end_simulator) {
	HeatModel model(N_test, L_test);

	double dt = 0.01;
	double t_end = 100.0;
	model.run(t_end, dt);

	// Verify that time vector was populated accurately by Model
	const auto &times = model.times();
	BOOST_CHECK(times.size() > 1);
	BOOST_CHECK_SMALL(times.front(), 1e-12);

	// Verify that slice() successfully pulled trajectory length
	auto temps = model.temperatures();
	BOOST_CHECK_EQUAL(temps.size(), times.size());
}

// Full end-to-end test through the Simulator, reset_before_run = false
BOOST_AUTO_TEST_CASE(heat_equation_end_to_end_simulator_no_reset) {
	HeatModel model(N_test, L_test);

	double dt = 0.01;
	model.run(50.0, dt, false); // current_time = 50.0

	// Run again without reset to continue from previous state
	model.run(50.0, dt, false); // current_time should now be 100.0

	const auto &times = model.times();
	BOOST_CHECK(times.size() > 1);

	auto temps = model.temperatures();
	BOOST_CHECK_EQUAL(temps.size(), times.size());

	// Last snapshot should match the exact solution at t_final = 100
	const auto &last_temps = temps.back();
	double t_final = times.back();
	double dx = L_test / (N_test + 1);
	for (int i = 0; i < N_test; ++i) {
		double x = (i + 1) * dx;
		double exact = HeatModel::exact_solution(x, t_final, L_test);
		BOOST_CHECK_CLOSE(last_temps[i], exact, 2.0);
	}
}
