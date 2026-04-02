#include "simulator.hpp"

void Simulator::add_subsystem(std::shared_ptr<System> sys, std::shared_ptr<TimeIntegrator> integrator) {
    subsystems.push_back({sys, integrator});
}

void Simulator::run(double t_end, double dt, bool reset_before_run) {
    if (reset_before_run) {
        reset(); // reset the simulator to its initial state and clear the trajectory before running
    }
    
    double t_stop = current_time + t_end;
    while (current_time < t_stop + 1e-12) {
        record(current_time);
        for (auto &[sys, integrator] : subsystems) {
            integrator->step(*sys, state, dt);
        }
        current_time += dt;
    }
}

void Simulator::record(double t) {
    timesteps.push_back(t);
    snapshots.push_back(state.y);
}

void Simulator::reset() {
    current_time = 0.0;
    state = initial_state;
    snapshots.clear();
    timesteps.clear();
}