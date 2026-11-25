#include "numeric_input.h"

#include <fstream>
#include <stdexcept>

namespace mbspec {

ResponseTable read_response(const std::string& filename, double cutoff, int time_digits)
{
    std::ifstream input(filename);
    if (!input)
        throw std::runtime_error("could not open response file '" + filename + "'");

    ResponseTable response;
    size_t line_number = 0;
    size_t samples = 0;
    double previous = 0.0;
    double previous_roundoff = 0.0;
    double interval_roundoff = 0.0;
    for (std::string line; std::getline(input, line);) {
        ++line_number;
        std::istringstream row(line);
        if ((row >> std::ws).eof())
            continue;

        double time, value;
        const auto location = filename + ':' + std::to_string(line_number) + ": ";
        if (!read_finite(row, time) || !read_finite(row, value) || !(row >> std::ws).eof())
            throw std::runtime_error(location + "expected exactly two finite numbers (time and response)");

        const double roundoff = printed_roundoff(time, time_digits);
        if (samples != 0) {
            const double interval = time - previous;
            if (!std::isfinite(interval) || interval <= 0.0)
                throw std::runtime_error(location + "response times must be finite and strictly increasing");
            const double error = roundoff + previous_roundoff;
            if (samples == 1) {
                response.interval = interval;
                interval_roundoff = error;
            } else if (std::abs(interval - response.interval) > error + interval_roundoff) {
                throw std::runtime_error(location + "response times must be evenly spaced");
            }
        }
        ++samples;
        previous = time;
        previous_roundoff = roundoff;
        if (time < cutoff)
            response.values.push_back(value);
    }
    if (input.bad())
        throw std::runtime_error("could not read response file '" + filename + "'");
    return response;
}

} // namespace mbspec
