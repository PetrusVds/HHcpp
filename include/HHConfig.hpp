#ifndef HH_CONFIG_HPP
#define HH_CONFIG_HPP

// Configuration parameters for the Hodgkin-Huxley model.
struct HHConfig {
    double C = 1.0;         // Membrane capacitance (uF/cm^2)
    double g_Na = 120.0;    // Max sodium conductance (mS/cm^2)
    double g_K = 36.0;      // Max potassium conductance (mS/cm^2)
    double g_leak = 0.3;    // Leak conductance (mS/cm^2)
    double E_Na = 50.0;     // Sodium reversal potential (mV)
    double E_K = -77.0;     // Potassium reversal potential (mV)
    double E_leak = -54.3;  // Leak reversal potential (mV)
    double i_ext = 0.2;     // Injected absolute current (nA)
    double a = 1.0;         // Compartment radius (um)
    double r_a = 100.0;     // Axial resistivity (Ohm cm)
};

#endif // HH_CONFIG_HPP