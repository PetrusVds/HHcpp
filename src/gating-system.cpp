#include "gating-system.hpp"
#include "time-integrator.hpp"

// Calculates the right-hand side of the gating variable ODEs...
void GatingSystem::rhs(const State &state, std::vector<double> &dydt) {
    double V, m, n, h;
    double am, bm, ah, bh, an, bn;

    for (int i = 0; i < N; ++i) {
        V = state[i];
        m = state[N + i];
        h = state[2 * N + i];
        n = state[3 * N + i];

        am = Compartment::alpha_m(V);
        bm = Compartment::beta_m(V);
        ah = Compartment::alpha_h(V);
        bh = Compartment::beta_h(V);
        an = Compartment::alpha_n(V);
        bn = Compartment::beta_n(V);

        dydt[i] = am * (1.0 - m) - bm * m;
        dydt[N + i] = ah * (1.0 - h) - bh * h;
        dydt[2 * N + i] = an * (1.0 - n) - bn * n;
    }
}

void GatingSystem::assemble_diagonal_system(const State &state, double dt, std::vector<double> &diag, std::vector<double> &rhs) {
    double V, m, n, h;
    double am, bm, ah, bh, an, bn;

    for (int i = 0; i < N; ++i) {
        V = state[i];

        am = Compartment::alpha_m(V);
        bm = Compartment::beta_m(V);
        ah = Compartment::alpha_h(V);
        bh = Compartment::beta_h(V);
        an = Compartment::alpha_n(V);
        bn = Compartment::beta_n(V);

        m = state[N + i];
        h = state[2 * N + i];
        n = state[3 * N + i];
        
        diag[i] = 1.0 + dt * (am + bm);
        rhs[i] = m + dt * am;
        diag[N + i] = 1.0 + dt * (ah + bh);
        rhs[N + i] = h + dt * ah;
        diag[2 * N + i] = 1.0 + dt * (an + bn);
        rhs[2 * N + i] = n + dt * an;
    }
}

// We assume voltage to be constant, which also allows us to implement an exponential update for the gating variables.
void GatingSystem::exponential_solve(State &state, double dt) {
    std::vector<double> ss(3);
    double V, m, n, h;
    double m_inf, h_inf, n_inf;
    double am, bm, ah, bh, an, bn;
    double tau_m, tau_h, tau_n;

    for (int i = 0; i < N; ++i) {
        V = state[i];
        Compartment::get_steady_state_gating_variables(V, ss);
        m_inf = ss[0];
        h_inf = ss[1];
        n_inf = ss[2];

        am = Compartment::alpha_m(V);
        bm = Compartment::beta_m(V);
        ah = Compartment::alpha_h(V);
        bh = Compartment::beta_h(V);
        an = Compartment::alpha_n(V);
        bn = Compartment::beta_n(V);

        m = state[N + i];
        h = state[2 * N + i];
        n = state[3 * N + i];

        tau_m = 1.0 / (am + bm);
        tau_h = 1.0 / (ah + bh);
        tau_n = 1.0 / (an + bn);

        state[N + i] = m_inf + (m - m_inf) * exp(-dt / tau_m);
        state[2 * N + i] = h_inf + (h - h_inf) * exp(-dt / tau_h);
        state[3 * N + i] = n_inf + (n - n_inf) * exp(-dt / tau_n);
    }
}

void DiagonalSystem::accept(TimeIntegrator &integrator, State &state, double dt) {
    integrator.visit(*this, state, dt);
}