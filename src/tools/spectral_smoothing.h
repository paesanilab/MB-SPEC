#ifndef SPECTRAL_SMOOTHING_H
#define SPECTRAL_SMOOTHING_H

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace spectral {

// Applies the Gaussian smoothing convention used by the spectroscopy programs:
// a normalized kernel truncated at 3 sigma, with local edge renormalization.
inline bool ApplyGaussianSmoothing(std::vector<double>& spectrum,
                                   const double wavenumber_spacing_cm1,
                                   const double fwhm_cm1)
{
    const double sigma_cm1 =
        fwhm_cm1 / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    const double sigma = sigma_cm1 / wavenumber_spacing_cm1;
    const size_t radius = static_cast<size_t>(std::ceil(3.0 * sigma));
    if (radius >= spectrum.size())
        return false;

    const size_t kernel_size = 2 * radius + 1;
    const double center = static_cast<double>(radius);
    std::vector<double> kernel(kernel_size);
    for (size_t k = 0; k < kernel_size; ++k) {
        const double x = (static_cast<double>(k) - center) / sigma;
        kernel[k] = std::exp(-0.5 * x * x);
    }
    const double kernel_sum = std::accumulate(kernel.begin(), kernel.end(), 0.0);
    for (double& value : kernel)
        value /= kernel_sum;

    std::vector<double> smoothed(spectrum.size());
    smoothed[0] = 0.0;
    for (size_t i = 1; i < spectrum.size(); ++i) {
        const size_t first = (i > radius) ? i - radius : 1;
        const size_t last = std::min(spectrum.size() - 1, i + radius);
        double weighted_sum = 0.0;
        double weight_sum = 0.0;
        for (size_t j = first; j <= last; ++j) {
            const double weight = kernel[radius + j - i];
            weighted_sum += weight * spectrum[j];
            weight_sum += weight;
        }
        smoothed[i] = weighted_sum / weight_sum;
    }

    spectrum.swap(smoothed);
    return true;
}

} // namespace spectral

#endif // SPECTRAL_SMOOTHING_H
