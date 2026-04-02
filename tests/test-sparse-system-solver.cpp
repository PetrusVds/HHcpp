#define BOOST_TEST_MODULE CompartmentTests
#include <boost/test/included/unit_test.hpp>

#include <vector>
#include "utils/sparse-system-solver.hpp"

// Test with the identity matrix, which should yield the same vector as the solution
BOOST_AUTO_TEST_CASE(test_identity_matrix) {
    std::vector<double> a = {0.0, 0.0, 0.0}; // Lower
    std::vector<double> b = {1.0, 1.0, 1.0}; // Main
    std::vector<double> c = {0.0, 0.0, 0.0}; // Upper
    std::vector<double> d = {1.0, 2.0, 3.0}; // RHS
    std::vector<double> x(3);

    SparseSystemSolver::thomas(a, b, c, d, x);

    // Check if the solution is close to [1, 2, 3]
    BOOST_CHECK_SMALL(x[0] - 1.0, 1e-12);
    BOOST_CHECK_SMALL(x[1] - 2.0, 1e-12);
    BOOST_CHECK_SMALL(x[2] - 3.0, 1e-12);
}

// Test with a simple tridiagonal system that has a known solution
BOOST_AUTO_TEST_CASE(test_example_system_1) {
    // Example System:
    // 2x + -1y +  0z = 1
    // -1x + 2y + -1z = 0
    // 0x + -1y +  2z = 1
    
    std::vector<double> a = {0.0, -1.0, -1.0}; // Lower
    std::vector<double> b = {2.0,  2.0,  2.0}; // Main
    std::vector<double> c = {-1.0, -1.0, 0.0}; // Upper
    std::vector<double> d = {1.0,  0.0,  1.0}; // RHS
    std::vector<double> x(3);

    SparseSystemSolver::thomas(a, b, c, d, x);

    // Check if the solution is close to [1, 1, 1]
    BOOST_CHECK_SMALL(x[0] - 1.0, 1e-12);
    BOOST_CHECK_SMALL(x[1] - 1.0, 1e-12);
    BOOST_CHECK_SMALL(x[2] - 1.0, 1e-12);
}

BOOST_AUTO_TEST_CASE(test_example_system_2) {
    // Example System:
    // 4x + -1y +  0z = 3
    // -1x + 4y + -1z = 5
    // 0x + -1y +  4z = 6
    
    std::vector<double> a = {0.0, -1.0, -1.0}; // Lower
    std::vector<double> b = {4.0,  4.0,  4.0}; // Main
    std::vector<double> c = {-1.0, -1.0, 0.0}; // Upper
    std::vector<double> d = {3.0,  5.0,  6.0}; // RHS
    std::vector<double> x(3);

    SparseSystemSolver::thomas(a, b, c, d, x);

    // I approximated my calculations by hand, so will relax the tolerance a bit here
    BOOST_CHECK_SMALL(x[0] - 1.268, 1e-3);
    BOOST_CHECK_SMALL(x[1] - 2.071, 1e-3);
    BOOST_CHECK_SMALL(x[2] - 2.018, 1e-3);
}

BOOST_AUTO_TEST_CASE(test_solution_written_with_offset) {
    std::vector<double> a = {0.0, -1.0, -1.0};
    std::vector<double> b = {2.0,  2.0,  2.0};
    std::vector<double> c = {-1.0, -1.0, 0.0};
    std::vector<double> d = {1.0,  0.0,  1.0};
    std::vector<double> x = {-5.0, -5.0, -5.0, -5.0, -5.0};

    SparseSystemSolver::thomas(a, b, c, d, x, 2);

    BOOST_CHECK_SMALL(x[2] - 1.0, 1e-12);
    BOOST_CHECK_SMALL(x[3] - 1.0, 1e-12);
    BOOST_CHECK_SMALL(x[4] - 1.0, 1e-12);
    BOOST_CHECK_SMALL(x[0] + 5.0, 1e-12);
    BOOST_CHECK_SMALL(x[1] + 5.0, 1e-12);
}