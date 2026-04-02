#ifndef MODEL_HPP
#define MODEL_HPP

#include "simulator.hpp"
#include "state.hpp"
#include "time-integrator.hpp"
#include "HHConfig.hpp"
#include <memory>
#include <vector>

// A model is a wrapper class that contains a Simulator and provides a convenient interface for running simulations and accessing the results. The HHModel class is a specific implementation of a Model that sets up the Hodgkin-Huxley system with the appropriate subsystems and initial conditions.
class Model {
public:
    virtual ~Model() = default;
    
    // Declaration only
    void run(double t_end, double dt, bool reset_before_run = true);

    const std::vector<double> &times() const { return sim.times(); }

protected:
    int N;
    Simulator sim;

    Model(int N, State initial_state) : N(N), sim(initial_state) {}

    std::vector<std::vector<double>> slice(int offset) const;
};

class HHModel : public Model {
public:
	HHModel(int N, double L, std::shared_ptr<TimeIntegrator> voltage_integrator = std::make_shared<ImplicitEuler>(), std::shared_ptr<TimeIntegrator> gating_integrator = std::make_shared<ImplicitEuler>(), std::shared_ptr<HHConfig> config = std::make_shared<HHConfig>());

	std::vector<std::vector<double>> voltages() const { return slice(0); }
	std::vector<std::vector<double>> ms() const { return slice(N); }
	std::vector<std::vector<double>> hs() const { return slice(2 * N); }
	std::vector<std::vector<double>> ns() const { return slice(3 * N); }

	// Helper function to create an initial State at the resting potential, with gating variables set to their steady-state values at that potential.
	static State make_steady_state(int N, const HHConfig &config);

private:
	std::shared_ptr<HHConfig> config;
};

#endif // MODEL_HPP