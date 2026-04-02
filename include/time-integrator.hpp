#ifndef TIME_INTEGRATOR_HPP
#define TIME_INTEGRATOR_HPP

#include "utils/sparse-system-solver.hpp"
#include "state.hpp"
#include "system.hpp"
#include <stdexcept>

// Forward declarations
class System;
class DiagonalSystem;
class TridiagonalSystem;

// An abstract base class for time integrators.
class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    virtual void step(System &system, State &state, double dt) = 0;

    virtual void visit(DiagonalSystem &system, State &state, double dt) {
        throw std::runtime_error("TimeIntegrator::visit: no visit method implemented for DiagonalSystem");
    }

    virtual void visit(TridiagonalSystem &system, State &state, double dt) {
        throw std::runtime_error("TimeIntegrator::visit: no visit method implemented for TridiagonalSystem");
    }
};

// A simple Forward Euler integrator.
class ForwardEuler : public TimeIntegrator {
public:
    using TimeIntegrator::visit; // unhide the visit methods, so we can test the throws

    void step(System &system, State &state, double dt) override;
};

// An Implicit Euler integrator that relies on the System to implement the implicit solve.
class ImplicitEuler : public TimeIntegrator {
public:
    void step(System &system, State &state, double dt) override;
    void visit(DiagonalSystem &system, State &state, double dt) override;
    void visit(TridiagonalSystem &system, State &state, double dt) override;
};

// An Exponential Euler integrator which is frequently used for gating variables in the Hodgkin-Huxley model
class ExponentialEuler : public TimeIntegrator {
public:
    using TimeIntegrator::visit; // unhide the visit methods, so we can test the throws

    void step(System &system, State &state, double dt) override;
    void visit(DiagonalSystem &system, State &state, double dt) override;
};

#endif // TIME_INTEGRATOR_HPP