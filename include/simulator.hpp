#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP

#include "state.hpp"
#include "system.hpp"
#include "time-integrator.hpp"
#include <memory>
#include <vector>

// The Simulator class is responsible for managing the state of the system, running the simulation loop, and recording the trajectory of the state over time. It allows you to add multiple subsystems (like VoltageSystem and GatingSystem) with their respective time integrators, and it handles the execution of these subsystems at each timestep.
class Simulator {
public:
    Simulator(State initial_state) : state(initial_state), initial_state(std::move(initial_state)) {}

    void add_subsystem(std::shared_ptr<System> sys, std::shared_ptr<TimeIntegrator> integrator);
    void run(double t_end, double dt, bool reset_before_run = true);

    const std::vector<std::vector<double>> &trajectory() const { return snapshots; }
    const std::vector<double> &times() const { return timesteps; }

private:
    double current_time = 0.0;

    void record(double t); // Records the current state and time into the trajectory. 
    void reset(); // Resets the simulator to its initial state and clears the recorded trajectory. 

    State state;
    State initial_state; // fallback to initial state on reset
    std::vector<std::pair<std::shared_ptr<System>, std::shared_ptr<TimeIntegrator>>> subsystems;
    std::vector<std::vector<double>> snapshots;
    std::vector<double> timesteps;
};

#endif // SIMULATOR_HPP