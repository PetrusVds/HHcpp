#define BOOST_TEST_MODULE VoltageSystemTests
#include <boost/test/included/unit_test.hpp>

#include "HHConfig.hpp"
#include "model.hpp"
#include "state.hpp"
#include "time-integrator.hpp"
#include "voltage-system.hpp"
#include <cmath>
#include <vector>

static constexpr HHConfig config; // default parameters for tests
static const double V_rest = Compartment::compute_V_rest(config).value();
static constexpr double L_test = 100.0;

// Helper: compute rhs for compartment 0 with a given VoltageSystem
static double rhs_compartment0(VoltageSystem &sys, const State &state, int N) {
	std::vector<double> dydt(N, 0.0);
	sys.rhs(state, dydt);
	return dydt[0];
}

BOOST_AUTO_TEST_CASE(voltage_rhs_formula_single_compartment) {
	const int N = 1;
	State state = HHModel::make_steady_state(N, config);

	double m = state[N];
	double h = state[2 * N];
	double n = state[3 * N];
	double membrane_current_density = Compartment::get_membrane_current_density(V_rest, m, h, n, config);
	double injected_current_density = Compartment::get_injected_current_density(config.i_ext, config.a, L_test);
	double expected_dvdt = (-membrane_current_density + injected_current_density) / config.C;

	VoltageSystem sys(N, L_test, config);
	std::vector<double> dydt(N, 0.0);
	sys.rhs(state, dydt);

	BOOST_CHECK_CLOSE(dydt[0], expected_dvdt, 1e-9);
}

BOOST_AUTO_TEST_CASE(voltage_rhs_axial_coupling_direction) {
	const int N = 2;
	State state = HHModel::make_steady_state(N, config);

	state[0] = V_rest + 20.0;

	VoltageSystem sys(N, L_test, config);
	std::vector<double> dydt(N, 0.0);
	sys.rhs(state, dydt);

	double injected_current_density = Compartment::get_injected_current_density(config.i_ext, config.a, L_test);

	double m1 = state[N + 1], h1 = state[2 * N + 1], n1 = state[3 * N + 1];
	double isolated_dvdt = (-Compartment::get_membrane_current_density(V_rest, m1, h1, n1, config) + 0.0) / config.C;
	BOOST_CHECK_GT(dydt[1], isolated_dvdt);

	double m0 = state[N], h0 = state[2 * N], n0 = state[3 * N];
	double isolated_dvdt0 = (-Compartment::get_membrane_current_density(state[0], m0, h0, n0, config) + injected_current_density) / config.C;
	BOOST_CHECK_LT(dydt[0], isolated_dvdt0);
}

BOOST_AUTO_TEST_CASE(voltage_rhs_uniform_state_no_net_axial) {
	const int N = 5;
	State state = HHModel::make_steady_state(N, config);

	VoltageSystem sys(N, L_test, config);
	std::vector<double> dydt(N, 0.0);
	sys.rhs(state, dydt);

	for (int i = 2; i < N; ++i) {
		BOOST_CHECK_CLOSE(dydt[i], dydt[1], 1e-9);
	}
	BOOST_CHECK_GT(dydt[0], dydt[1]);
}

BOOST_AUTO_TEST_CASE(forward_euler_voltage_converges) {
	const int N = 1;
	State state = HHModel::make_steady_state(N, config);

	VoltageSystem sys(N, L_test, config);
	std::vector<double> dydt_initial(N, 0.0);
	sys.rhs(state, dydt_initial);
	BOOST_CHECK_GT(dydt_initial[0], 0.0);

	ForwardEuler fe;
	double dt = 0.001;
	double t = 0.0;
	for (int step = 0; step < 50000; ++step) {
		fe.step(sys, state, dt);
		t += dt;
	}

	BOOST_CHECK_GT(state[0], V_rest);
}

BOOST_AUTO_TEST_CASE(implicit_euler_voltage_stable_large_dt) {
	const int N = 3;
	State state = HHModel::make_steady_state(N, config);

	VoltageSystem sys(N, L_test, config);
	ImplicitEuler ie;

	double dt = 100.0;
	double t = 0.0;
	for (int step = 0; step < 20; ++step) {
		ie.step(sys, state, dt);
		t += dt;
	}

	for (int i = 0; i < N; ++i) {
		BOOST_CHECK(state[i] > -200.0 && state[i] < 200.0);
	}
}

BOOST_AUTO_TEST_CASE(implicit_solve_single_compartment_formula) {
	const int N = 1;
	const double dt = 0.5;
	State state = HHModel::make_steady_state(N, config);

	double V_old = state[0];
	double m = state[N];
	double h = state[2 * N];
	double n = state[3 * N];

	double membrane_current_density = Compartment::get_membrane_current_density(V_old, m, h, n, config);
	double membrane_conductance_density = Compartment::get_membrane_conductance_density(m, h, n, config);
	double injected_current_density = Compartment::get_injected_current_density(config.i_ext, config.a, L_test);
	double b = config.C / dt + membrane_conductance_density;
	double d = (config.C / dt) * V_old + injected_current_density - membrane_current_density + membrane_conductance_density * V_old;
	double V_expected = d / b;

	VoltageSystem sys(N, L_test, config);
	ImplicitEuler ie;
	ie.step(sys, state, dt);

	BOOST_CHECK_CLOSE(state[0], V_expected, 1e-9);
}

BOOST_AUTO_TEST_CASE(implicit_solve_spatial_decay) {
	const int N = 5;
	const double dt = 1.0;
	State state = HHModel::make_steady_state(N, config);

	VoltageSystem sys(N, L_test, config);
	ImplicitEuler ie;
	ie.step(sys, state, dt);

	for (int k = 0; k < N - 1; ++k) {
		BOOST_CHECK_GT(state[k], state[k + 1]);
	}
}

BOOST_AUTO_TEST_CASE(cached_i_ext_density_auto_updates) {
	const int N = 1;
	HHConfig mutable_config;
	VoltageSystem sys(N, L_test, mutable_config);
	State state = HHModel::make_steady_state(N, config);

	double rhs_before = rhs_compartment0(sys, state, N);

	// Modify i_ext directly, cache should auto-refresh on next rhs() call
	mutable_config.i_ext *= 10.0;

	double rhs_after = rhs_compartment0(sys, state, N);

	// rhs should differ because the lazy cache detected the change
	BOOST_CHECK_NE(rhs_before, rhs_after);
}

BOOST_AUTO_TEST_CASE(cached_g_axial_auto_updates) {
	const int N = 2;
	HHConfig mutable_config;
	VoltageSystem sys(N, L_test, mutable_config);

	// Non-uniform voltage so axial current is nonzero
	State state = HHModel::make_steady_state(N, config);
	state[0] = V_rest + 20.0;

	std::vector<double> dydt_before(N, 0.0);
	sys.rhs(state, dydt_before);

	// Modify r_a directly, cache should auto-refresh
	mutable_config.r_a *= 5.0;

	std::vector<double> dydt_after(N, 0.0);
	sys.rhs(state, dydt_after);

	// rhs should now differ for the compartment receiving axial current
	BOOST_CHECK_NE(dydt_before[1], dydt_after[1]);
}

BOOST_AUTO_TEST_CASE(cached_values_auto_update_when_a_changes) {
	const int N = 2;
	HHConfig mutable_config;
	VoltageSystem sys(N, L_test, mutable_config);

	// Non-uniform voltage so axial current is nonzero
	State state = HHModel::make_steady_state(N, config);
	state[0] = V_rest + 20.0;

	std::vector<double> dydt_before(N, 0.0);
	sys.rhs(state, dydt_before);

	// Change radius directly, both g_axial and i_ext_density should auto-refresh
	mutable_config.a *= 2.0;

	std::vector<double> dydt_after(N, 0.0);
	sys.rhs(state, dydt_after);

	// Both compartments should be affected (comp 0 via i_ext + axial, comp 1 via axial)
	BOOST_CHECK_NE(dydt_before[0], dydt_after[0]);
	BOOST_CHECK_NE(dydt_before[1], dydt_after[1]);
}

BOOST_AUTO_TEST_CASE(implicit_solve_auto_uses_updated_cache) {
	const int N = 1;
	const double dt = 0.5;
	HHConfig mutable_config;
	VoltageSystem sys(N, L_test, mutable_config);

	State state1 = HHModel::make_steady_state(N, config);
	ImplicitEuler ie;
	ie.step(sys, state1, dt);
	double V_implicit_before = state1[0];

	// Change i_ext directly, implicit_solve should auto-refresh cache
	mutable_config.i_ext *= 10.0;

	State state2 = HHModel::make_steady_state(N, config);
	ie.step(sys, state2, dt);
	double V_implicit_after = state2[0];

	BOOST_CHECK_NE(V_implicit_before, V_implicit_after);
}

BOOST_AUTO_TEST_CASE(cache_not_refreshed_when_config_unchanged) {
	const int N = 1;
	HHConfig mutable_config;
	VoltageSystem sys(N, L_test, mutable_config);
	State state = HHModel::make_steady_state(N, config);

	double rhs1 = rhs_compartment0(sys, state, N);
	double rhs2 = rhs_compartment0(sys, state, N);

	BOOST_CHECK_CLOSE(rhs1, rhs2, 1e-12);
}