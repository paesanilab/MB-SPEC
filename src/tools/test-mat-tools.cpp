#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#include "mat_tools.h"

namespace {

constexpr double tolerance = 1e-12;

void require_close(double actual, double expected, const char* description)
{
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance) {
        std::cerr << description << ": expected " << expected
                  << ", got " << actual << '\n';
        throw std::runtime_error(description);
    }
}

void test_inversion()
{
    // Nonsymmetric, with a zero first pivot to exercise row interchange.
    const std::array<double, 9> original = {0, 2, 1, 1, 3, 0, 0, 0, 4};
    const std::array<double, 9> expected = {-1.5, 1, 0.375,
                                          0.5, 0, -0.125,
                                          0, 0, 0.25};
    auto inverse = original;
    tools::invert(3, inverse.data());
    for (size_t i = 0; i < inverse.size(); ++i) {
        require_close(inverse[i], expected[i], "matrix inverse");
    }

    // Calculate the residual independently of the production matrix helpers.
    for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
            double product = 0;
            for (size_t k = 0; k < 3; ++k) {
                product += original[3*i + k] * inverse[3*k + j];
            }
            require_close(product, i == j ? 1 : 0, "inverse identity residual");
        }
    }
}

void test_singular_matrix()
{
    double singular[] = {1, 2, 2, 4};
    try {
        tools::invert(2, singular);
    } catch (const std::runtime_error&) {
        return;
    }
    throw std::runtime_error("singular inversion did not throw");
}

void test_eigenvalues()
{
    // The 2x2 block has eigenvalues 1 and 3; the last diagonal adds 5.
    std::array<double, 9> matrix = {2, 1, 0, 1, 2, 0, 0, 0, 5};
    const auto original = matrix;
    std::array<double, 3> eigenvalues{};
    if (tools::eigenvalues(3, matrix.data(), eigenvalues.data()) != 0) {
        throw std::runtime_error("eigensolver failed");
    }
    const std::array<double, 3> expected = {1, 3, 5};
    for (size_t i = 0; i < expected.size(); ++i) {
        require_close(eigenvalues[i], expected[i], "symmetric eigenvalue");
    }
    if (matrix != original) {
        throw std::runtime_error("eigenvalues modified its input matrix");
    }
}

} // namespace

int main()
{
    try {
        test_inversion();
        test_singular_matrix();
        test_eigenvalues();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
