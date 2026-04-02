#include "voltage-system.hpp"
#include "time-integrator.hpp"

void VoltageSystem::rhs(const State &state, std::vector<double> &dydt) {
    refresh_if_needed();

    double V, m, h, n;
    double membrane_current_density, axial_current_density;

    for (int i = 0; i < N; ++i) {
        V = state[i];
        m = state[N + i];
        h = state[2 * N + i];
        n = state[3 * N + i];

        membrane_current_density = Compartment::get_membrane_current_density(V, m, h, n, config);
        axial_current_density = 0.0;

        if (i > 0) axial_current_density += axial_conductance_density * (state[i - 1] - V);
        if (i < N - 1) axial_current_density += axial_conductance_density * (state[i + 1] - V);

        dydt[i] = (-membrane_current_density + ((i == 0) ? injected_current_density : 0.0) + axial_current_density) / config.C; 
    }
}

void VoltageSystem::assemble_tridiagonal_system(const State &state, double dt, std::vector<double> &lower, std::vector<double> &diag, std::vector<double> &upper, std::vector<double> &rhs) {
    refresh_if_needed();

    double V, m, h, n;
    double membrane_current_density, membrane_conductance_density;

    for (int i = 0; i < N; ++i) {
        V = state[i];
        m = state[N + i];
        h = state[2 * N + i];
        n = state[3 * N + i];

        membrane_current_density = Compartment::get_membrane_current_density(V, m, h, n, config);
        membrane_conductance_density = Compartment::get_membrane_conductance_density(m, h, n, config);

        if (i > 0) lower[i] = -axial_conductance_density;
        diag[i] = config.C / dt + axial_conductance_density * ((i > 0) + (i < N - 1)) + membrane_conductance_density;
        if (i < N - 1) upper[i] = -axial_conductance_density;
        rhs[i] = (config.C / dt) * V + ((i == 0) ? injected_current_density : 0.0) - membrane_current_density + membrane_conductance_density * V; 
    }
}

void VoltageSystem::_init_cache() {
    refresh_cached_axial_conductance_density();
    refresh_cached_injected_current_density();
}

void VoltageSystem::refresh_cached_axial_conductance_density() {
    axial_conductance_density = Compartment::get_axial_conductance_density(config.a, config.r_a, L_um);
    cached_a = config.a;
    cached_r_a = config.r_a;
}

void VoltageSystem::refresh_cached_injected_current_density() {
    injected_current_density = Compartment::get_injected_current_density(config.i_ext, config.a, L_um);
    cached_a = config.a;
    cached_i_ext = config.i_ext;
}

void VoltageSystem::refresh_if_needed() {
    bool radius_changed = (config.a != cached_a);
    bool axial_resistivity_changed = (config.r_a != cached_r_a);
    bool injected_current_changed = (config.i_ext != cached_i_ext);

    if (radius_changed || axial_resistivity_changed) refresh_cached_axial_conductance_density();
    if (radius_changed || injected_current_changed) refresh_cached_injected_current_density();
}

void TridiagonalSystem::accept(TimeIntegrator &integrator, State &state, double dt) {
    integrator.visit(*this, state, dt);
}