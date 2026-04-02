#ifndef VOLTAGE_SYSTEM_HPP
#define VOLTAGE_SYSTEM_HPP

#include "system.hpp"
#include "utils/compartment.hpp"
#include <cassert>

// The VoltageSystem is responsible for updating the voltage of each compartment based on the current gating variables, as described in the assignment. It does not update the gating variables themselves; that is the responsibility of the GatingSystem! The VoltageSystem also implements an implicit solver for the voltage update, which requires solving a tridiagonal system of equations at each timestep. This is done using the Thomas algorithm implemented in utils.hpp.
class VoltageSystem : public TridiagonalSystem {
public:
	VoltageSystem(int num_compartments, double L, const HHConfig &config) : N(num_compartments), L_um(L), config(config) {
		_init_cache(); // initialize the cache based on the initial config values
	}
	void rhs(const State &state, std::vector<double> &dydt) override;
    
    void assemble_tridiagonal_system(const State &state, double dt, std::vector<double> &lower, std::vector<double> &diag, std::vector<double> &upper, std::vector<double> &rhs) override;

	int get_size() const override { return N; } // number of voltages in the state vector
	int get_offset() const override { return 0; } // offset of the voltage variables in the state vector (they start at index 0)

private:
	int N;
	double L_um;            // compartment length [µm]
	const HHConfig &config; // reference to the configuration struct
	double axial_conductance_density;
	double injected_current_density;

	// last-seen config values
    double cached_a;
    double cached_r_a;
    double cached_i_ext;

    void _init_cache();
    void refresh_cached_axial_conductance_density();
    void refresh_cached_injected_current_density();
    void refresh_if_needed();
};

#endif // VOLTAGE_SYSTEM_HPP
