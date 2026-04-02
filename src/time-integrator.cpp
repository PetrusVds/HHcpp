#include "time-integrator.hpp"
#include "system.hpp"
#include "gating-system.hpp"

// --- Forward Euler Implementation ---
void ForwardEuler::step(System &system, State &state, double dt) {
    int size = system.get_size();
    int offset = system.get_offset();

    std::vector<double> dydt(size);
    system.rhs(state, dydt);

    for (int i = 0; i < size; ++i) {
        state[offset + i] += dt * dydt[i];
    }
}

// --- Implicit Euler Implementation ---
void ImplicitEuler::step(System &system, State &state, double dt) {
    system.accept(*this, state, dt); // calls visit for correct system type
}

void ImplicitEuler::visit(DiagonalSystem &system, State &state, double dt) {
    int size = system.get_size();
    int offset = system.get_offset();

    std::vector<double> diag(size, 0.0);
    std::vector<double> rhs(size, 0.0);
    system.assemble_diagonal_system(state, dt, diag, rhs);

    for (int i = 0; i < size; ++i) {
        state.y[offset + i] = rhs[i] / diag[i];
    }
}

void ImplicitEuler::visit(TridiagonalSystem &system, State &state, double dt) {
    int size = system.get_size();
    int offset = system.get_offset();

    std::vector<double> lower(size, 0.0);
    std::vector<double> diag(size, 0.0);
    std::vector<double> upper(size, 0.0);
    std::vector<double> rhs(size, 0.0);
    system.assemble_tridiagonal_system(state, dt, lower, diag, upper, rhs);
    SparseSystemSolver::thomas(lower, diag, upper, rhs, state.y, offset);
}

// --- Exponential Euler Implementation ---

void ExponentialEuler::step(System &system, State &state, double dt) {
    system.accept(*this, state, dt);
}

void ExponentialEuler::visit(DiagonalSystem &system, State &state, double dt) {
    GatingSystem &gating_system = dynamic_cast<GatingSystem&>(system);
    gating_system.exponential_solve(state, dt);
}