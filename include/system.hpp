#ifndef SYSTEM_HPP
#define SYSTEM_HPP

#include "state.hpp"
#include "time-integrator.hpp"
#include <vector>

class TimeIntegrator; // forward declaration

// The System class is an abstract base class that defines the interface for any system of equations we want to solve. 
class System {
public:
    virtual ~System() = default;

    // The rhs() function computes the right-hand side of the system of equations, i.e. the derivatives of the state variables. The dydt vector should be filled with the computed derivatives.
    virtual void rhs(const State &state, std::vector<double> &dydt) = 0;
    virtual void accept(TimeIntegrator &integrator, State &state, double dt) = 0;

    virtual int get_size() const = 0; // Returns the number of state variables to be updated by the integrator
    virtual int get_offset() const = 0; // Returns the index offset in the state vector where this system's variables start (useful for systems that only update a subset of the state variables)
};

// Our GatingSystem is Diagonal
class DiagonalSystem : public System {
public:
    virtual void assemble_diagonal_system(const State &state, double dt, std::vector<double> &diag, std::vector<double> &rhs) = 0;
    void accept(TimeIntegrator &integrator, State &state, double dt) override;
};

// Our VoltageSystem is Tridiagonal
class TridiagonalSystem : public System {
public:
    virtual void assemble_tridiagonal_system(const State &state, double dt, std::vector<double> &lower, std::vector<double> &diag, std::vector<double> &upper, std::vector<double> &rhs) = 0;
    void accept(TimeIntegrator &integrator, State &state, double dt) override;
};

#endif // SYSTEM_HPP
