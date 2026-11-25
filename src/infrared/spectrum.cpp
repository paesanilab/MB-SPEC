#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>

#include <cassert>
#include <cstdlib>
#include <sstream>

#include "cli.h"
#include "numeric_input.h"
#include "system.h"
#include "constants.h"
#include "spectral_smoothing.h"

namespace {

//----------------------------------------------------------------------------//
// Physical constants and numerical settings
//----------------------------------------------------------------------------//

const double light_velocity = constants::c0; // m/s
const double speed_of_light_cm_s = constants::c0*constants::cm_per_m; // cm/s

} // namespace

int main(int argc, char** argv)
{
    //----------------------------------------------------------------------------//
    // Program setup and JSON configuration
    //----------------------------------------------------------------------------//

    if (const auto status = mbspec::check_arguments(argc, argv, 3, "<response.dat> <mb-spec.json>"))
        return *status;

    System sys;
    sys.SetUpFromJson("INFRARED", argv[2]);

    double T = sys.GetTemperature();
    double volume = sys.GetBoxVolume();
    double tau = sys.GetTau();
    double smoothing_fwhm_cm1 = sys.GetSmoothingFwhmCm1();
    bool has_smoothing_fwhm_cm1 = sys.HasSmoothingFwhmCm1();
    double max_int = sys.GetIntMax();
    const double transform_time_ps = sys.GetSpectrumTransformTimePs();
    const bool has_transform_time_ps = sys.HasSpectrumTransformTimePs();

    {
        if (T <= 0.0) {
            std::cerr << "INFRARED.temperature must be positive (K)" << std::endl;
            return EXIT_FAILURE;
        }
    }
    if (!has_transform_time_ps || transform_time_ps <= 0.0 ||
        transform_time_ps > max_int) {
        std::cerr << "INFRARED.spectrum_transform_time_ps must satisfy "
                  << "0 < spectrum_transform_time_ps <= int_max" << std::endl;
        return EXIT_FAILURE;
    }
    {
        if (tau < 0.0) {
            std::cerr << "INFRARED.tau must be nonnegative (ps)" << std::endl;
            return EXIT_FAILURE;
        }
    }
    {
        if (volume <= 0.0) {
            std::cerr << "INFRARED.box must define a positive simulation-cell volume (A^3)"
                      << std::endl;
            return EXIT_FAILURE;
        }
        volume *= constants::A3_to_m3; // A^3 -> m^3
    }
    {
        if (!std::is_same<decltype(tau), double>::value) {
            std::cerr << "could not convert '" << tau
                      << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    }
    {
        if (!std::is_same<decltype(max_int), double>::value) {
            std::cerr << "could not convert '" << max_int
                      << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    }

    //----------------------------------------------------------------------------//
    // Response-function input and time-grid validation
    //----------------------------------------------------------------------------//

    mbspec::ResponseTable response;
    try {
        response = mbspec::read_response(argv[1], max_int, 9);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    auto& corr = response.values;
    const size_t n_response = corr.size();
    if (n_response < 2) {
        std::cerr << "response file '" << argv[1] << "' must contain at least 2 time points below INFRARED.int_max\n";
        return EXIT_FAILURE;
    }
    const double dt_ps = response.interval;
    const double available_time_ps = n_response*dt_ps;
    if (transform_time_ps > available_time_ps + 1.0e-10*dt_ps) {
        std::cerr << "INFRARED.spectrum_transform_time_ps exceeds the available "
                  << "response duration (" << available_time_ps << " ps)" << std::endl;
        return EXIT_FAILURE;
    }
    const size_t Ndata = static_cast<size_t>(
        std::floor(transform_time_ps/dt_ps + 1.0e-10));
    if (Ndata < 2) {
        std::cerr << "INFRARED.spectrum_transform_time_ps must include at least "
                  << "two response intervals" << std::endl;
        return EXIT_FAILURE;
    }
    std::cerr << "Ndata : " << Ndata << std::endl;

    //----------------------------------------------------------------------------//
    // Correlation truncation and exponential window
    //----------------------------------------------------------------------------//

    double arg = 2.0 * M_PI / (double(2 * Ndata - 1));

    double kT = constants::kB * T;
    double beta = 1.0 / kT;

    double dt = dt_ps * constants::ps_to_s;
    double dnu = 1.0 / (double(2 * Ndata - 1)) / dt;
    const double wavenumber_spacing_cm1 = dnu / speed_of_light_cm_s;

    if (has_smoothing_fwhm_cm1 && smoothing_fwhm_cm1 < 0.0) {
        std::cerr << "INFRARED.smoothing_fwhm_cm1 must be nonnegative"
                  << std::endl;
        return EXIT_FAILURE;
    }

    if (tau != 0)
        for (size_t t = 0; t < Ndata; ++t)
            corr[t] *= std::exp(-(t * dt) / (tau * constants::ps_to_s));

    //----------------------------------------------------------------------------//
    // Fourier transform of the dipole TCF and IR absorption prefactor
    //
    // n(omega) alpha(omega)
    //   = beta omega^2 / (6 epsilon_0 c V)
    //     * integral dt exp(-i omega t) C_mu_mu(t)
    //
    // No additional quantum correction factor is applied here.
    // The interpretation of C_mu_mu(t) as classical or Kubo-transformed
    // depends on the dynamics used to generate the correlation function.
    //----------------------------------------------------------------------------//

    std::vector<double> freq;
    std::vector<double> n_alpha_spectrum;

    for (size_t inu = 0; inu < Ndata; ++inu) {

        double frequency_hz = dnu * inu;
        double omega = 2.0 * M_PI * frequency_hz;

        double cosine_transform = 0.0;

        for (size_t itime = Ndata - 1; itime > 0; --itime) {

            double nu_time = (double(itime)) * (double(inu)) * arg;
            cosine_transform += std::cos(nu_time) * corr[itime];
        }

        for (size_t itime = 0; itime < Ndata; ++itime) {

            double nu_time = (double(itime)) * (double(inu)) * arg;
            cosine_transform += std::cos(nu_time) * corr[itime];
        }

        double ir_prefactor =
            beta * omega * omega /
            (6.0 * constants::epsilon_vacuum * light_velocity * volume);

        double n_alpha =
            cosine_transform * ir_prefactor * dt
            * constants::m_per_cm;

        freq.push_back(frequency_hz / speed_of_light_cm_s);
        n_alpha_spectrum.push_back(n_alpha);
    }

    //----------------------------------------------------------------------------//
    // Optional Gaussian spectral smoothing
    //----------------------------------------------------------------------------//

    const bool use_smoothing =
        has_smoothing_fwhm_cm1 && smoothing_fwhm_cm1 > 0.0;

    if (use_smoothing) {
        if (!spectral::ApplyGaussianSmoothing(
                n_alpha_spectrum, wavenumber_spacing_cm1, smoothing_fwhm_cm1)) {
            std::cerr << "INFRARED.smoothing_fwhm_cm1 is too broad for "
                      << "the available spectral grid" << std::endl;
            return EXIT_FAILURE;
        }

        //----------------------------------------------------------------------------//
        // Smoothed spectrum output
        //----------------------------------------------------------------------------//

        for (size_t n = 0; n < n_alpha_spectrum.size(); ++n) {
            std::cout << freq[n] << '\t' << n_alpha_spectrum[n] << std::endl;
        }
    }

    else {

        //----------------------------------------------------------------------------//
        // Spectrum output
        //----------------------------------------------------------------------------//

        for (size_t n = 0; n < n_alpha_spectrum.size(); ++n) {
            std::cout << freq[n] << '\t' << n_alpha_spectrum[n] << std::endl;
        }
    }

    return 0;
}
