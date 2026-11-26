#ifndef MB_SPEC_NUMERIC_INPUT_H
#define MB_SPEC_NUMERIC_INPUT_H

#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace mbspec {

inline bool read_finite(std::istream& input, double& value)
{
    std::string field;
    if (!(input >> field))
        return false;
    std::istringstream number(field);
    return bool(number >> value) && number.eof() && std::isfinite(value);
}

inline double printed_roundoff(double value, int digits)
{
    if (value == 0.0)
        return 0.0;
    return 0.5 * std::pow(10.0, std::floor(std::log10(std::abs(value))) + 1 - digits)
         + 4 * std::numeric_limits<double>::epsilon() * std::abs(value);
}

struct ResponseTable {
    std::vector<double> values;
    double interval = 0.0;
};

ResponseTable read_response(const std::string& filename, double cutoff, int time_digits);

} // namespace mbspec

#endif
