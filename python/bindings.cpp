#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "model.hpp"
#include "state.hpp"
#include "time-integrator.hpp"

namespace py = pybind11;

PYBIND11_MODULE(HHcpp, m) {
	m.doc() = "C++ implementation of the Hodgkin-Huxley model supporting the simulation of single branch neurons with multiple compartments. All parameters are homogeneous across compartments, except for the external current; this is only applied to the first compartment.";

	py::class_<HHConfig, std::shared_ptr<HHConfig>>(m, "HHConfig",
		"Configuration of homogeneous parameters for the Hodgkin-Huxley model.\n\n"
		"Example:\n"
		"    config = HHcpp.HHConfig()\n"
		"    config.g_Na = 100.0\n"
		"    config.i_ext = 0.5\n"
		"    model = HHcpp.HHModel(N=10, L=100.0, config=config)\n")
		.def(py::init<>(),"Create a config with default Hodgkin-Huxley parameters.")
		.def_readwrite("C", &HHConfig::C, "Membrane capacitance (uF/cm^2)")
		.def_readwrite("g_Na", &HHConfig::g_Na, "Max sodium conductance (mS/cm^2)")
		.def_readwrite("g_K", &HHConfig::g_K, "Max potassium conductance (mS/cm^2)")
		.def_readwrite("g_leak", &HHConfig::g_leak, "Leak conductance (mS/cm^2)")
		.def_readwrite("E_Na", &HHConfig::E_Na, "Sodium reversal potential (mV)")
		.def_readwrite("E_K", &HHConfig::E_K, "Potassium reversal potential (mV)")
		.def_readwrite("E_leak", &HHConfig::E_leak, "Leak reversal potential (mV)")
		.def_readwrite("i_ext", &HHConfig::i_ext, "Injected absolute current (nA)")
		.def_readwrite("a", &HHConfig::a, "Compartment radius (um)")
		.def_readwrite("r_a", &HHConfig::r_a, "Axial resistivity (Ohm cm)");

	py::class_<TimeIntegrator, std::shared_ptr<TimeIntegrator>>(m, "TimeIntegrator", "Abstract base class for time integrators.)");

	py::class_<ForwardEuler, TimeIntegrator, std::shared_ptr<ForwardEuler>>(m, "ForwardEuler", "Explicit forward Euler method for time integration.")
		.def(py::init<>());

	py::class_<ImplicitEuler, TimeIntegrator, std::shared_ptr<ImplicitEuler>>(m, "ImplicitEuler", "Implicit Euler method for time integration.")
		.def(py::init<>());

	py::class_<ExponentialEuler, TimeIntegrator, std::shared_ptr<ExponentialEuler>>(m, "ExponentialEuler", "Exponential Euler method for time integration. Note: This method is only applicable to the gating variable updates.")
		.def(py::init<>());

	py::class_<HHModel>(m, "HHModel",
		"Single branch Hodgkin-Huxley model.\n\n"
		"Example:\n"
		"    model = HHcpp.HHModel(N=10, L=100.0)\n\n"
		"    - Option 1: Run a single simulation\n"
		"    model.run(t_end=50.0, dt=0.01)\n"
		"    times = np.array(model.times())   # time series of time points\n"
		"    Vs = np.array(model.voltages())   # time series of voltages\n "
		"   ms = np.array(model.ms())         # time series of m gating variables\n"
		"    hs = np.array(model.hs())         # time series of h gating variables\n"
		"    ns = np.array(model.ns())         # time series of n gating variables\n\n"
		"    - Option 2: Run multiple simulations with the same model instance\n"
		"    model.run(t_end=100.0, dt=0.01)   # t_begin = 0.0, t_end = 100.0\n"
		"    # run a third simulation without resetting the state\n"
		"    model.run(t_end=150.0, dt=0.01, reset_before_run=False) # t_begin = 100.0, t_end = 250.0\n"
		"    times = np.array(model.times())   # time series of time points\n"
		"    Vs = np.array(model.voltages())   # time series of voltages\n "
		"   ms = np.array(model.ms())         # time series of m gating variables\n"
		"    hs = np.array(model.hs())         # time series of h gating variables\n"
		"    ns = np.array(model.ns())         # time series of n gating variables")
		.def(py::init<int, double, std::shared_ptr<TimeIntegrator>, std::shared_ptr<TimeIntegrator>, std::shared_ptr<HHConfig>>(),
			py::arg("N"), py::arg("L"),
			py::arg("voltage_integrator") = std::make_shared<ImplicitEuler>(),
			py::arg("gating_integrator") = std::make_shared<ImplicitEuler>(),
			py::arg("config") = std::make_shared<HHConfig>(),
			"Construct a Hodgkin-Huxley model for a single branch with N compartments and compartment size L (um). Defining custom time integrators and config is optional.")
		.def("run", &HHModel::run, py::arg("t_end"), py::arg("dt"), py::arg("reset_before_run") = true,
			"Run the simulation.\n\n"
			"Args:\n"
			"    t_end: End time of the simulation (ms).\n"
			"    dt: Time step size (ms).\n"
			"    reset_before_run: Whether to reset the simulator to its initial state before running. This allows you to run multiple simulations with the same model instance without manually resetting.")
		.def("times", &HHModel::times,
			"Return the recorded time points (ms).")
		.def("voltages", &HHModel::voltages,
			"Return the voltage trajectory.")
		.def("ms", &HHModel::ms,
			"Return the m gating variable trajectory")
		.def("hs", &HHModel::hs,
			"Return the h gating variable trajectory.")
		.def("ns", &HHModel::ns,
			"Return the n gating variable trajectory.");
}