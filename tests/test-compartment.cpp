#define BOOST_TEST_MODULE CompartmentTests
#include <boost/test/included/unit_test.hpp>

#include "HHConfig.hpp"
#include "utils/compartment.hpp"
#include <cmath>
#include <vector>

// Test for calculating membrane current
BOOST_AUTO_TEST_CASE(test_membrane_current) {
	HHConfig config;
	auto V = Compartment::compute_V_rest(config);
	BOOST_REQUIRE(V.has_value());

	double am = (0.1 * (*V + 40)) / (1 - exp(-(*V + 40) / 10));
	double ah = 0.07 * exp(-(*V + 65) / 20);
	double an = (0.01 * (*V + 55)) / (1 - exp(-(*V + 55) / 10));
	double bm = 4.0 * exp(-(*V + 65) / 18);
	double bh = 1 / (1 + exp(-(*V + 35) / 10));
	double bn = 0.125 * exp(-(*V + 65) / 80);

	double m = am / (am + bm);
	double h = ah / (ah + bh);
	double n = an / (an + bn);

	double I_Na = config.g_Na * pow(m, 3) * h * (*V - config.E_Na);
	double I_K = config.g_K * pow(n, 4) * (*V - config.E_K);
	double I_leak = config.g_leak * (*V - config.E_leak);
	double I_total = I_Na + I_K + I_leak;

	BOOST_CHECK_SMALL(Compartment::get_membrane_current_density(*V, m, h, n, config) - I_total, 1e-12);
}

BOOST_AUTO_TEST_CASE(test_steady_state_gating_variables) {
	HHConfig config;
	auto V = Compartment::compute_V_rest(config);
	BOOST_REQUIRE(V.has_value());

	std::vector<double> gating_variables(3);
	Compartment::get_steady_state_gating_variables(*V, gating_variables);

	double am = (0.1 * (*V + 40)) / (1 - exp(-(*V + 40) / 10));
	double ah = 0.07 * exp(-(*V + 65) / 20);
	double an = (0.01 * (*V + 55)) / (1 - exp(-(*V + 55) / 10));
	double bm = 4.0 * exp(-(*V + 65) / 18);
	double bh = 1 / (1 + exp(-(*V + 35) / 10));
	double bn = 0.125 * exp(-(*V + 65) / 80);

	double m = am / (am + bm);
	double h = ah / (ah + bh);
	double n = an / (an + bn);

	BOOST_CHECK_SMALL(gating_variables[0] - m, 1e-12);
	BOOST_CHECK_SMALL(gating_variables[1] - h, 1e-12);
	BOOST_CHECK_SMALL(gating_variables[2] - n, 1e-12);
	BOOST_CHECK_EQUAL(gating_variables.size(), 3);
}

BOOST_AUTO_TEST_CASE(test_transition_rates) {
	HHConfig config;
	auto V = Compartment::compute_V_rest(config);
	BOOST_REQUIRE(V.has_value());

	BOOST_CHECK_SMALL(Compartment::alpha_m(*V) - (0.1 * (*V + 40.0)) / (1.0 - exp(-(*V + 40.0) / 10.0)), 1e-12);
	BOOST_CHECK_SMALL(Compartment::alpha_h(*V) - 0.07 * exp(-(*V + 65) / 20), 1e-12);
	BOOST_CHECK_SMALL(Compartment::alpha_n(*V) - (0.01 * (*V + 55.0)) / (1.0 - exp(-(*V + 55.0) / 10.0)), 1e-12);
	BOOST_CHECK_SMALL(Compartment::beta_m(*V) - 4.0 * exp(-(*V + 65) / 18), 1e-12);
	BOOST_CHECK_SMALL(Compartment::beta_h(*V) - 1.0 / (1.0 + exp(-(*V + 35) / 10)), 1e-12);
	BOOST_CHECK_SMALL(Compartment::beta_n(*V) - 0.125 * exp(-(*V + 65) / 80), 1e-12);
}