#ifndef COMPARTMENT_HPP
#define COMPARTMENT_HPP

#include "HHConfig.hpp"
#include "utils/root-finder.hpp"
#include <cassert>
#include <cmath>
#include <optional>
#include <functional>
#include <stdexcept>
#include <vector>

namespace Compartment {
	// Unit conversions to improve readability of the code and avoid magic numbers.
	constexpr double UM_TO_CM = 1e-4;
	constexpr double NA_TO_UA = 1e-3;
	constexpr double S_TO_MS = 1e3;

	// Compute the membrane area in cm^2 given the radius and length in microns.
	inline double membrane_area_cm2(double radius_um, double length_um) {
		return 2.0 * M_PI * UM_TO_CM * radius_um * UM_TO_CM * length_um;
	}

	// Injected current must be converted into current density, which in our project is uA/cm^2.
	inline double get_injected_current_density(double injected_current_nA, double radius_um, double length_um) {
		return (injected_current_nA * NA_TO_UA) / membrane_area_cm2(radius_um, length_um);
	}

	// The assignment defines the axial conductance as mS/cm^2
	inline double get_axial_conductance_density(double radius_um, double axial_resistivity, double length_um) {
		return radius_um  / (2.0 * axial_resistivity * length_um * length_um * UM_TO_CM) * S_TO_MS; // a factor UM_TO_CM has been removed from numerator and denominator, since they cancel out
	}

	// Voltage dependent transition rates for gating variables
	inline double alpha_m(double V) { return (0.1 * (V + 40.0)) / (1.0 - exp(-(V + 40.0) / 10.0)); }
	inline double alpha_h(double V) { return 0.07 * exp(-(V + 65) / 20); }
	inline double alpha_n(double V) { return (0.01 * (V + 55.0)) / (1.0 - exp(-(V + 55.0) / 10.0)); }
	inline double beta_m(double V) { return 4.0 * exp(-(V + 65) / 18); }
	inline double beta_h(double V) { return 1.0 / (1.0 + exp(-(V + 35) / 10)); }
	inline double beta_n(double V) { return 0.125 * exp(-(V + 65) / 80); }

	// Compute the total membrane current density (in uA/cm^2) for a given voltage and gating variables
	inline double get_membrane_current_density(double V, double m, double h, double n, const HHConfig &config) {
		return (config.g_Na * pow(m, 3) * h * (V - config.E_Na)) // current density of Na
				+ (config.g_K * pow(n, 4) * (V - config.E_K))    // current density of K
				+ (config.g_leak * (V - config.E_leak));         // current density of leak
	}

	// Calculate the total membrane conductance density (in mS/cm^2) for a given gating variable state
	inline double get_membrane_conductance_density(double m, double h, double n, const HHConfig &config) {
		return (config.g_Na * pow(m, 3) * h) + (config.g_K * pow(n, 4)) + config.g_leak;
	}

	// Compute the steady-state values of the gating variables m, h, n for a given voltage V and store them in the provided vector ss (which must have size 3).
	inline void get_steady_state_gating_variables(double V, std::vector<double> &ss) {
		ss[0] = alpha_m(V) / (alpha_m(V) + beta_m(V)); // m_inf 
		ss[1] = alpha_h(V) / (alpha_h(V) + beta_h(V)); // h_inf
		ss[2] = alpha_n(V) / (alpha_n(V) + beta_n(V)); // n_inf
	}

	// Compute the resting potential (i.e. the voltage at which the membrane current is zero) using Newton-Raphson, starting from an initial guess V0 (default -65 mV). This condition was based on resting potentials found in the literature. The exact value in this case is around -64.97 mV for the default parameters, which tells us that the initial guess is very close.
	std::optional<double> compute_V_rest(const HHConfig &config, double V0 = -65.0);

} // namespace Compartment

#endif // COMPARTMENT_HPP