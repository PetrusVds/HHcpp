#include "utils/compartment.hpp"

std::optional<double> Compartment::compute_V_rest(const HHConfig &config, double V0) {
    std::vector<double> ss(3); // allocate once
    return RootFinder::newton_raphson(
        [&config, &ss](double V) { // lambda for the function whose root we want to find
            get_steady_state_gating_variables(V, ss);
            return get_membrane_current_density(V, ss[0], ss[1], ss[2], config);
        }, V0);
}