#define BOOST_TEST_MODULE HHModelTests
#include <boost/test/included/unit_test.hpp>

#include "model.hpp"
#include "HHConfig.hpp"
#include "utils/compartment.hpp"

constexpr HHConfig default_config;

constexpr int N_test = 10;
constexpr double L_test = 100.0;
constexpr double dt_test = 0.01;
constexpr double LOWER_VOLTAGE_BOUND = default_config.E_K;
constexpr double UPPER_VOLTAGE_BOUND = default_config.E_Na;
constexpr double LOWER_GATING_BOUND = 0.0;
constexpr double UPPER_GATING_BOUND = 1.0;

// Trivial test to verify that the initial state of the HHModel is set up correctly, with voltages at the expected resting potential and gating variables at their steady-state values for that voltage.
BOOST_AUTO_TEST_CASE(initial_state) {
	HHConfig config;
	State hs = HHModel::make_steady_state(N_test, config);

	// Check that the state has the expected size
	BOOST_CHECK(hs.total_size() == 4 * N_test);

	// Check that the initial voltage is close to the expected resting potential
	auto V_rest = Compartment::compute_V_rest(config);
	BOOST_REQUIRE(V_rest.has_value());
	for (int i = 0; i < N_test; ++i) {
		BOOST_CHECK_CLOSE(hs.y[i], *V_rest, 1e-5);
	}

	// Check that the initial gating variables are close to their steady-state values at V_rest
	std::vector<double> ss(3);
	Compartment::get_steady_state_gating_variables(*V_rest, ss);
	for (int i = 0; i < N_test; ++i) {
		BOOST_CHECK_CLOSE(hs.y[N_test + i], ss[0], 1e-5);
		BOOST_CHECK_CLOSE(hs.y[2*N_test + i], ss[1], 1e-5);
		BOOST_CHECK_CLOSE(hs.y[3*N_test + i], ss[2], 1e-5);
	}
}

BOOST_AUTO_TEST_CASE(run_simulation) {
	HHModel model(N_test, L_test);
	model.run(0.1, dt_test);

	// Check that the time vector was populated
	const auto &times = model.times();
	BOOST_CHECK(times.size() > 1);
	BOOST_CHECK_SMALL(times.front(), 1e-12);

	// Check that the voltage trajectory was populated and has the expected shape
	const auto &voltages = model.voltages();
	BOOST_CHECK(voltages.size() == times.size());
	for (const auto &snap : voltages) {
		BOOST_CHECK(snap.size() == N_test);
	}

	// Check that the gating variable trajectories were populated and have the
	// expected shape
	const auto &ms = model.ms();
	const auto &hs = model.hs();
	const auto &ns = model.ns();
	BOOST_CHECK(ms.size() == times.size());
	BOOST_CHECK(hs.size() == times.size());
	BOOST_CHECK(ns.size() == times.size());
	for (const auto &snap : ms) {
		BOOST_CHECK(snap.size() == N_test);
	}
	for (const auto &snap : hs) {
		BOOST_CHECK(snap.size() == N_test);
	}
	for (const auto &snap : ns) {
		BOOST_CHECK(snap.size() == N_test);
	}

	// All gating variables should be between 0 and 1, and voltages should be
	// within a reasonable range (e.g. -100 mV to 100 mV) for the duration of the
	// simulation (note: this is of course only true when the integrators are
	// chosen to be stable for the given timestep, but since we're using Implicit
	// Euler by default, it should be fine for this test).
	for (const auto &snap : voltages) {
		for (double V : snap) {
			BOOST_CHECK_LE(V, UPPER_VOLTAGE_BOUND);
			BOOST_CHECK_GE(V, LOWER_VOLTAGE_BOUND);
		}
	}
	for (const auto &snap : ms) {
		for (double m : snap) {
			BOOST_CHECK_LE(m, UPPER_GATING_BOUND);
			BOOST_CHECK_GE(m, LOWER_GATING_BOUND);
		}
	}
	for (const auto &snap : hs) {
		for (double h : snap) {
			BOOST_CHECK_LE(h, UPPER_GATING_BOUND);
			BOOST_CHECK_GE(h, LOWER_GATING_BOUND);
		}
	}
	for (const auto &snap : ns) {
		for (double n : snap) {
			BOOST_CHECK_LE(n, UPPER_GATING_BOUND);
			BOOST_CHECK_GE(n, LOWER_GATING_BOUND);
		}
	}
}

BOOST_AUTO_TEST_CASE(config_i_ext_changes_simulation) {
	const int N = 3;
	const double L = 100.0;

	auto config1 = std::make_shared<HHConfig>();
	HHModel model1(N, L, std::make_shared<ImplicitEuler>(), std::make_shared<ImplicitEuler>(), config1);
	model1.run(1.0, 0.01);
	auto V1 = model1.voltages();

	auto config2 = std::make_shared<HHConfig>();
	config2->i_ext = 5.0; // much larger current
	HHModel model2(N, L, std::make_shared<ImplicitEuler>(), std::make_shared<ImplicitEuler>(), config2);
	model2.run(1.0, 0.01);
	auto V2 = model2.voltages();

	// Final voltages should differ due to different i_ext
	BOOST_CHECK_NE(V1.back()[0], V2.back()[0]);
}

BOOST_AUTO_TEST_CASE(config_r_a_changes_simulation) {
	const int N = 3;
	const double L = 100.0;

	auto config1 = std::make_shared<HHConfig>();
	HHModel model1(N, L, std::make_shared<ImplicitEuler>(), std::make_shared<ImplicitEuler>(), config1);
	model1.run(1.0, 0.01);
	auto V1 = model1.voltages();

	auto config2 = std::make_shared<HHConfig>();
	config2->r_a = 500.0; // much higher axial resistance
	HHModel model2(N, L, std::make_shared<ImplicitEuler>(), std::make_shared<ImplicitEuler>(), config2);
	model2.run(1.0, 0.01);
	auto V2 = model2.voltages();

	// Distal compartment should differ (less current reaches it with higher r_a)
	BOOST_CHECK_NE(V1.back()[N - 1], V2.back()[N - 1]);
}

BOOST_AUTO_TEST_CASE(config_a_changes_simulation) {
	const int N = 3;
	const double L = 100.0;

	auto config1 = std::make_shared<HHConfig>();
	HHModel model1(N, L, std::make_shared<ImplicitEuler>(), std::make_shared<ImplicitEuler>(), config1);
	model1.run(1.0, 0.01);
	auto V1 = model1.voltages();

	auto config2 = std::make_shared<HHConfig>();
	config2->a = 5.0; // much larger radius
	HHModel model2(N, L, std::make_shared<ImplicitEuler>(), std::make_shared<ImplicitEuler>(), config2);
	model2.run(1.0, 0.01);
	auto V2 = model2.voltages();

	BOOST_CHECK_NE(V1.back()[0], V2.back()[0]);
}

BOOST_AUTO_TEST_CASE(config_change_between_runs) {
	const int N = 3;
	const double L = 100.0;

	auto config = std::make_shared<HHConfig>();
	HHModel model(N, L, std::make_shared<ImplicitEuler>(), std::make_shared<ImplicitEuler>(), config);
	model.run(1.0, 0.01);
	auto V_before = model.voltages();

	// Change config between runs — VoltageSystem should auto-detect
	config->i_ext = 5.0;
	model.run(1.0, 0.01);
	auto V_after = model.voltages();

	BOOST_CHECK_NE(V_before.back()[0], V_after.back()[0]);
}

BOOST_AUTO_TEST_CASE(throws_error_when_compute_V_rest_fails) {
	auto config = std::make_shared<HHConfig>();
	config->E_K = config->E_Na; // E_K and E_Na must be different for compute_V_rest to work
	BOOST_CHECK_THROW(HHModel::make_steady_state(N_test, *config), std::runtime_error);
}