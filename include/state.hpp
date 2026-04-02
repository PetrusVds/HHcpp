#ifndef STATE_HPP
#define STATE_HPP

#include <vector>

// The State struct is a simple wrapper around a vector of doubles, representing the state of the system at any given time. Just to make the code more readable.
struct State {
	std::vector<double> y;

	State() = default;
	explicit State(int size) : y(size) {}

	int total_size() const { return y.size(); }

	double &operator[](int i) { return y[i]; }
	const double &operator[](int i) const { return y[i]; }
};

#endif // STATE_HPP