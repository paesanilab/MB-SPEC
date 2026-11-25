#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>

#include "cli.h"
#include "system.h"
#include "constants.h"
#include "spectral_smoothing.h"

namespace {

const double pi2 = 2.0*M_PI;
const double speed_of_light_cm_ps =
    constants::c0*constants::cm_per_m*constants::ps_to_s; // cm/ps

} // namespace

int main(int argc, char** argv)
{
    //----------------------------------------------------------------------------//
    // Program setup and JSON configuration
    //----------------------------------------------------------------------------//

    if (const auto status = mbspec::check_arguments(argc, argv, 3, "<response.dat> <mb-spec.json>"))
        return *status;

    System sys;
    sys.SetUpFromJson("SFG", argv[2]);

    const double temperature = sys.GetTemperature();
    const double tau = sys.GetTau();
    const double int_max = sys.GetIntMax();
    const double smoothing_fwhm_cm1 = sys.GetSmoothingFwhmCm1();
    const bool has_smoothing_fwhm_cm1 = sys.HasSmoothingFwhmCm1();
    const double maximum_wavenumber_cm1 =
        sys.GetSpectrumMaxWavenumberCm1();
    const size_t num_wavenumbers = sys.GetSpectrumNumWavenumbers();
    const bool subtract_initial_response = sys.GetSubtractInitialResponse();
    const bool has_subtract_initial_response =
        sys.HasSubtractInitialResponse();
    const double transform_time_ps = sys.GetSpectrumTransformTimePs();
    const bool has_transform_time_ps = sys.HasSpectrumTransformTimePs();

    if (temperature <= 0.0) {
        std::cerr << "SFG.temperature must be positive (K)" << std::endl;
        return EXIT_FAILURE;
    }
    if (tau < 0.0) {
        std::cerr << "SFG.tau must be nonnegative (ps)" << std::endl;
        return EXIT_FAILURE;
    }
    if (int_max <= 0.0) {
        std::cerr << "SFG.int_max must be positive (ps)" << std::endl;
        return EXIT_FAILURE;
    }
    if (maximum_wavenumber_cm1 <= 0.0) {
        std::cerr << "SFG.spectrum_max_wavenumber_cm1 must be positive "
                  << "(cm^-1)" << std::endl;
        return EXIT_FAILURE;
    }
    if (num_wavenumbers < 2) {
        std::cerr << "SFG.spectrum_num_wavenumbers must be at least 2"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (!has_subtract_initial_response) {
        std::cerr << "SFG.subtract_initial_response must be specified"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (!has_transform_time_ps || transform_time_ps <= 0.0) {
        std::cerr << "SFG.spectrum_transform_time_ps must be positive"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (transform_time_ps > int_max) {
        std::cerr << "SFG.spectrum_transform_time_ps must not exceed SFG.int_max"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (has_smoothing_fwhm_cm1 && smoothing_fwhm_cm1 < 0.0) {
        std::cerr << "SFG.smoothing_fwhm_cm1 must be nonnegative"
                  << std::endl;
        return EXIT_FAILURE;
    }

    //----------------------------------------------------------------------------//
    // Response input and time-grid validation
    //----------------------------------------------------------------------------//

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "could not open '" << argv[1] << "' for reading"
                  << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<double> time;
    std::vector<double> corr;
    std::string line;
    size_t line_number = 0;
    while (std::getline(infile, line)) {
        ++line_number;
        std::istringstream line_stream(line);
        double t = 0.0;
        double response = 0.0;
        if (!(line_stream >> t >> response)) {
            std::cerr << "response file line " << line_number
                      << " must begin with time and response values"
                      << std::endl;
            return EXIT_FAILURE;
        }
        time.push_back(t);
        corr.push_back((t > int_max) ? 0.0 : response);
    }

    const size_t n_data = time.size();
    if (n_data < 2) {
        std::cerr << "response file must contain at least two time points"
                  << std::endl;
        return EXIT_FAILURE;
    }

    const double dt = time[1] - time[0];
    if (dt <= 0.0) {
        std::cerr << "response time values must be strictly increasing"
                  << std::endl;
        return EXIT_FAILURE;
    }
    for (size_t i = 2; i < n_data; ++i) {
        const double interval = time[i] - time[i - 1];
        if (interval <= 0.0) {
            std::cerr << "response time values must be strictly increasing"
                      << std::endl;
            return EXIT_FAILURE;
        }
        if (std::abs(interval - dt) > 1.0e-10*dt) {
            std::cerr << "response time values must be evenly spaced"
                      << std::endl;
            return EXIT_FAILURE;
        }
    }

    if (subtract_initial_response) {
        const double corr_0 = corr[0];
        for (size_t i = 0; i < n_data; ++i) {
            corr[i] -= corr_0;
        }
    }

    //----------------------------------------------------------------------------//
    // Time-domain window
    //----------------------------------------------------------------------------//

    if (tau > 0.0) {
        for (size_t i = 0; i < n_data; ++i) {
            corr[i] *= std::exp(-time[i]/tau);
        }
    }

    const double available_time_ps = n_data*dt;
    if (transform_time_ps > available_time_ps + 1.0e-10*dt) {
        std::cerr << "SFG.spectrum_transform_time_ps exceeds the available "
                  << "response duration (" << available_time_ps << " ps)"
                  << std::endl;
        return EXIT_FAILURE;
    }
    const size_t n_transform = static_cast<size_t>(
        std::floor(transform_time_ps/dt + 1.0e-10));
    if (n_transform < 2) {
        std::cerr << "SFG.spectrum_transform_time_ps must include at least "
                  << "two response intervals"
                  << std::endl;
        return EXIT_FAILURE;
    }

    const double wavenumber_factor = 1.0/speed_of_light_cm_ps;
    const double dnu =
        (maximum_wavenumber_cm1/wavenumber_factor)/(num_wavenumbers - 1);
    const double wavenumber_spacing_cm1 = wavenumber_factor*dnu;

    //----------------------------------------------------------------------------//
    // Half-response Fourier transform on the fixed SFG frequency grid
    //----------------------------------------------------------------------------//

    std::vector<double> freq;
    std::vector<double> intensity_im;
    std::vector<double> intensity_re;
    freq.reserve(num_wavenumbers);
    intensity_im.reserve(num_wavenumbers);
    intensity_re.reserve(num_wavenumbers);

    for (size_t k = 0; k < num_wavenumbers; ++k) {
        const double nu = k*dnu;
        double accum_im = 0.5*corr[0];
        double accum_re = 0.0;
        for (size_t i = 1; i < n_transform - 1; ++i) {
            accum_im += std::cos(pi2*nu*time[i])*corr[i];
            accum_re += std::sin(pi2*nu*time[i])*corr[i];
        }

        const double wavenumber = wavenumber_factor*nu;
        freq.push_back(wavenumber);
        intensity_im.push_back(dt*wavenumber*accum_im);
        intensity_re.push_back(-dt*wavenumber*accum_re);
    }

    const bool use_smoothing =
        has_smoothing_fwhm_cm1 && smoothing_fwhm_cm1 > 0.0;

    //----------------------------------------------------------------------------//
    // Optional Gaussian spectral smoothing
    //----------------------------------------------------------------------------//

    if (use_smoothing) {
        if (!spectral::ApplyGaussianSmoothing(
                intensity_im, wavenumber_spacing_cm1, smoothing_fwhm_cm1) ||
            !spectral::ApplyGaussianSmoothing(
                intensity_re, wavenumber_spacing_cm1, smoothing_fwhm_cm1)) {
            std::cerr << "SFG.smoothing_fwhm_cm1 is too broad for "
                      << "the available spectral grid" << std::endl;
            return EXIT_FAILURE;
        }
    }

    //----------------------------------------------------------------------------//
    // Spectrum output
    //----------------------------------------------------------------------------//

    for (size_t i = 0; i < num_wavenumbers; ++i) {
        const double magnitude2 =
            intensity_re[i]*intensity_re[i] + intensity_im[i]*intensity_im[i];
        std::cout << std::setw(15) << freq[i]
                  << std::setw(15) << intensity_im[i]
                  << std::setw(15) << intensity_re[i]
                  << std::setw(15) << magnitude2 << std::endl;
    }

    return 0;
}
