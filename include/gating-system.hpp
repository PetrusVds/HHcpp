#ifndef GATING_SYSTEM_HPP
#define GATING_SYSTEM_HPP

#include "system.hpp"
#include "utils/compartment.hpp"
#include <cmath>


// The GatingSystem is responsible for updating the gating variables (m, h, n) based on the current voltage. It does not update the voltage itself; that is the responsibility of the VoltageSystem!
class GatingSystem : public DiagonalSystem {
public:
	GatingSystem(int num_compartments) : N(num_compartments) {}

	void rhs(const State &state, std::vector<double> &dydt) override;
    
    void assemble_diagonal_system(const State &state, double dt, std::vector<double> &diag, std::vector<double> &rhs) override;

	void exponential_solve(State &state, double dt);

	int get_size() const override { return 3 * N; } // number of gating variables in the state vector (3 per compartment)
	int get_offset() const override { return N; } // offset of the gating variables in the state vector (they start after the N voltages), order: m, h, n.

private:
	int N;
};

#endif // GATING_SYSTEM_HPP
