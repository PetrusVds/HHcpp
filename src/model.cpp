#include "model.hpp"
#include "gating-system.hpp"
#include "voltage-system.hpp"
#include "utils/compartment.hpp"
#include <stdexcept>

// --- Model Base Class Implementation ---

void Model::run(double t_end, double dt, bool reset_before_run) {
    sim.run(t_end, dt, reset_before_run);
}

std::vector<std::vector<double>> Model::slice(int offset) const {
    std::vector<std::vector<double>> out;
    for (const auto &snap : sim.trajectory()) {
        out.push_back({snap.begin() + offset, snap.begin() + offset + N});
    }
    return out;
}

// --- HHModel Implementation ---

HHModel::HHModel(int N, double L, std::shared_ptr<TimeIntegrator> voltage_integrator,  std::shared_ptr<TimeIntegrator> gating_integrator, std::shared_ptr<HHConfig> config) : Model(N, make_steady_state(N, *config)), config(config) {
    sim.add_subsystem(std::make_shared<VoltageSystem>(N, L, *this->config), voltage_integrator);
    sim.add_subsystem(std::make_shared<GatingSystem>(N), gating_integrator);
}

State HHModel::make_steady_state(int N, const HHConfig &config) {
    auto V_rest = Compartment::compute_V_rest(config);
    if (!V_rest.has_value()) {
        throw std::runtime_error("make_steady_state: failed to compute resting potential");
    }
    
    std::vector<double> ss(3);
    Compartment::get_steady_state_gating_variables(*V_rest, ss);

    State hs(4 * N);
    for (int i = 0; i < N; ++i) {
        hs[i] = *V_rest;
        hs[N + i] = ss[0];     // m
        hs[2 * N + i] = ss[1]; // h
        hs[3 * N + i] = ss[2]; // n
    }
    return hs;
}