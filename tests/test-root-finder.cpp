#define BOOST_TEST_MODULE RootFinderTests
#include <boost/test/included/unit_test.hpp>

#include "utils/root-finder.hpp"

// Test with a simple quadratic function f(x) = (x - 2)^2, which has a root at x = 2
BOOST_AUTO_TEST_CASE(test_convergence_to_root_for_quadratic) {
    auto root = RootFinder::newton_raphson(
        [](double x) { 
            return pow(x - 2, 2); 
        }, 0.0);
    BOOST_REQUIRE(root.has_value());
    BOOST_CHECK_CLOSE(*root, 2.0, 1e-3); // bit relaxed tolerance since the root is a double root, which can cause slow convergence
}

// Test with a cubic function f(x) = (x - 1)(x - 2)(x - 3), which has roots at x = 1, 2, and 3
BOOST_AUTO_TEST_CASE(test_convergence_to_root_for_trigonometric) {
    auto f = [](double x) { 
        return (x - 1) * (x - 2) * (x - 3); 
    };

    auto root1 = RootFinder::newton_raphson(f, 0.5);
    auto root2 = RootFinder::newton_raphson(f, 1.75);
    auto root3 = RootFinder::newton_raphson(f, 2.25);
    auto root4 = RootFinder::newton_raphson(f, 3.5);
    auto root5 = RootFinder::newton_raphson(f, 2.5);
    auto root6 = RootFinder::newton_raphson(f, 3.0);

    BOOST_REQUIRE(root1.has_value());
    BOOST_REQUIRE(root2.has_value());
    BOOST_REQUIRE(root3.has_value());
    BOOST_REQUIRE(root4.has_value());
    BOOST_REQUIRE(root5.has_value());
    BOOST_REQUIRE(root6.has_value());

    BOOST_CHECK_CLOSE(*root1, 1.0, 1e-4);
    BOOST_CHECK_CLOSE(*root2, 2.0, 1e-4);
    BOOST_CHECK_CLOSE(*root3, 2.0, 1e-4);
    BOOST_CHECK_CLOSE(*root4, 3.0, 1e-4);
    BOOST_CHECK_CLOSE(*root5, 1.0, 1e-4); // edge case, choice of initial condition is bad
    BOOST_CHECK_CLOSE(*root6, 3.0, 1e-4);
}

BOOST_AUTO_TEST_CASE(newton_raphson_zero_derivative_test) {
    auto result = RootFinder::newton_raphson(
        [](double x) { 
            return 5.0; // constant function, derivative is zero everywhere
        }, 1.0);
    
    BOOST_CHECK(!result.has_value());
}

BOOST_AUTO_TEST_CASE(newton_raphson_convergence_fail_test) {
    auto result = RootFinder::newton_raphson(
        [](double x) { 
            return (x * x) - 1.0; // root at x = 1 and x = -1, but we will limit max_iter to 1 to force failure
        }, 100.0, 1e-10, 1);
    
    BOOST_CHECK(!result.has_value());
}