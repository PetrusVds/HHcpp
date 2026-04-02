#include "HHConfig.hpp"
#include "model.hpp"
#include "gating-system.hpp"
#include "state.hpp"
#include "time-integrator.hpp"
#include "utils/compartment.hpp"
#include "voltage-system.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <optional>

// A very simple main program to demonstrate a mock usage of our Hodgkin-Huxley implementation. This is not meant to be a comprehensive test or benchmark, just a sanity check that the code runs and produces reasonable output. The real use-case for this project will be importing it into python to get high-performance simulations with a nice interface.
int main() {
	const int N = 10;       // number of compartments
	const double L = 20.0;  // compartment length [µm]
	const double dt = 0.01; // timestep [ms]
	const double T = 50.0;  // total simulation time [ms]
	const int steps = static_cast<int>(T / dt);
	const int print_every = 100; // print every 1 ms
	const HHConfig config;
	auto V_rest = Compartment::compute_V_rest(config);
	if (!V_rest.has_value()) {
		std::cerr << "Error: failed to compute resting potential. Check config parameters.\n";
		return 1;
	}

	State state = HHModel::make_steady_state(N, config);
	VoltageSystem voltage_sys(N, L, config);
	GatingSystem gating_sys(N);
	ForwardEuler fe;
	ImplicitEuler ie;

	std::cout << "=== Hodgkin-Huxley Parameters ===\n";
	std::cout << std::fixed << std::setprecision(6);
	std::cout << "  C        = " << config.C << " uF/cm^2\n";
	std::cout << "  g_Na     = " << config.g_Na << " mS/cm^2  (0.12 S/cm^2)\n";
	std::cout << "  g_K      = " << config.g_K << " mS/cm^2  (0.036 S/cm^2)\n";
	std::cout << "  g_leak   = " << config.g_leak << " mS/cm^2  (0.0003 S/cm^2)\n";
	std::cout << "  E_Na     = " << config.E_Na << " mV\n";
	std::cout << "  E_K      = " << config.E_K << " mV\n";
	std::cout << "  E_leak   = " << config.E_leak << " mV\n";
	std::cout << "  i_ext    = " << config.i_ext << " nA \n";
	std::cout << "  a        = " << config.a << " um \n";
	std::cout << "  r_a      = " << config.r_a << " Ohm cm\n";
	{
		double injected_current_density = Compartment::get_injected_current_density(config.i_ext, config.a, L);
		double axial_conductance_density = Compartment::get_axial_conductance_density(config.a, config.r_a, L);
		std::cout << "  i_ext_density (computed) = " << injected_current_density << " uA/cm^2\n";
		std::cout << "  g_axial       (computed) = " << axial_conductance_density << " mS/cm^2\n";
	}
	std::cout << "\n";

	// Verify steady-state gating variables at rest
	{
		std::vector<double> ss(3);
		Compartment::get_steady_state_gating_variables(*V_rest, ss);
		std::cout << "=== Steady-state gating at V_rest = " << *V_rest << " mV ===\n";
		std::cout << "  m_inf = " << ss[0] << "\n";
		std::cout << "  h_inf = " << ss[1] << "\n";
		std::cout << "  n_inf = " << ss[2] << "\n";
		double membrane_current_density = Compartment::get_membrane_current_density(*V_rest, ss[0], ss[1], ss[2], config);
		std::cout << "  I_ion at rest = " << membrane_current_density << " uA/cm^2  (should be ~0 for resting potential)\n\n";
	}

	//  Simulation
	std::cout << "=== Simulation (dt=" << dt << " ms, T=" << T << " ms, print every " << (print_every * dt) << " ms) ===\n";
	std::cout << std::setw(8) << "t [ms]" 
				<< std::setw(12) << "V[0] [mV]"
				<< std::setw(10) << "m[0]" 
				<< std::setw(10) << "h[0]"
				<< std::setw(10) << "n[0]" << "\n";
	std::cout << std::string(50, '-') << "\n";

	double t = 0.0;
	for (int step = 0; step < steps; ++step) {
		if (step % print_every == 0) {
			std::cout << std::setw(8) << t 
						<< std::setw(12) << state[0]
						<< std::setw(10) << state[N + 0] 
						<< std::setw(10) << state[2 * N + 0] 
						<< std::setw(10) << state[3 * N + 0] << "\n";
		}
		ie.step(voltage_sys, state, dt);
		ie.step(gating_sys, state, dt);
		t += dt;
	}

	std::cout << std::setw(8) << t 
				<< std::setw(12) << state[0] 
				<< std::setw(10) << state[N + 0] 
				<< std::setw(10) << state[2 * N + 0] 
				<< std::setw(10) << state[3 * N + 0] << "\n";

	std::cout << "\n=== Done ===\n";

	return 0;
}
