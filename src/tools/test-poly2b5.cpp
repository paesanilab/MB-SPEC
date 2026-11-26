#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "poly2b5.h"
#include "test-poly2b5-data.h"

namespace {

constexpr std::array<size_t, 6> sizes = {1826, 1826, 1826, 1826, 1052, 1052};
constexpr std::array<const char*, 6> names = {"pHa1", "pHa2", "pHb1", "pHb2", "pOa", "pOb"};
constexpr double guard = 1234567.25;

void make_input(int sample, std::array<double, 15>& h, std::array<double, 15>& o)
{
    for (int i = 0; i < 15; ++i) {
        h[i] = (sample == 0 ? (3*i)%7 + 1 : (5*i)%13 - 6) / 8.0;
        o[i] = (sample == 0 ? (7*i)%13 - 6 : (2*i)%7 + 1) / 8.0;
    }
}

void test_case(int sample, const std::array<std::int32_t, 9408>* expected)
{
    std::array<double, 15> h{}, o{};
    if (expected) {
        make_input(sample, h, o);
    }
    const auto original_h = h;
    const auto original_o = o;
    std::array<std::vector<double>, 6> output;
    for (size_t block = 0; block < output.size(); ++block) {
        output[block].assign(sizes[block] + 2, std::numeric_limits<double>::quiet_NaN());
        output[block].front() = output[block].back() = guard;
    }

    poly::cpt_poly2b5(h.data(), o.data(), output[0].data()+1, output[1].data()+1,
                     output[2].data()+1, output[3].data()+1,
                     output[4].data()+1, output[5].data()+1);

    size_t index = 0;
    for (size_t block = 0; block < output.size(); ++block) {
        if (output[block].front() != guard || output[block].back() != guard) {
            throw std::runtime_error(std::string(names[block]) + " wrote outside its output array");
        }
        for (size_t i = 0; i < sizes[block]; ++i, ++index) {
            // These reference inputs are binary fractions. All intermediate
            // products and sums are exact, so no numerical tolerance is needed.
            const double reference = expected ? (*expected)[index] / 32768.0 : 0;
            const double actual = output[block][i+1];
            if (!std::isfinite(actual) || actual != reference) {
                throw std::runtime_error(std::string(names[block]) + "[" + std::to_string(i)
                                         + "] differs from the original polynomial");
            }
        }
    }
    if (h != original_h || o != original_o) {
        throw std::runtime_error("polynomial modified its input arrays");
    }
}

void test_overlapping_outputs()
{
    std::array<double, 15> h, o;
    make_input(0, h, o);
    std::array<double, 1828> output;
    output.fill(std::numeric_limits<double>::quiet_NaN());
    output.front() = output.back() = guard;
    double* shared = output.data()+1;
    poly::cpt_poly2b5(h.data(), o.data(), shared, shared, shared, shared, shared, shared);
    // Preserve the original write order even when output arrays overlap:
    // pOb writes last up to 1051; pHb2 writes last for the remaining entries.
    for (size_t i = 0; i < 1826; ++i) {
        const size_t index = i < 1052 ? 8356+i : 5478+i;
        if (shared[i] != expected_0[index] / 32768.0) {
            throw std::runtime_error("overlapping outputs differ from the original polynomial");
        }
    }
    if (output.front() != guard || output.back() != guard) {
        throw std::runtime_error("overlapping outputs wrote outside their array");
    }
}

} // namespace

int main()
{
    try {
        test_case(0, nullptr);
        test_case(0, &expected_0);
        test_case(1, &expected_1);
        test_overlapping_outputs();
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
