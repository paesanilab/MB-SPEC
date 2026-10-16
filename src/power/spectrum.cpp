#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <random>

#include <cstdlib>
#include <sstream>

#include "system.h"
#include "constants.h"
#include "spectral_smoothing.h"

namespace {

const double conv_freq = constants::c0*constants::cm_per_m; // cm/s

} // namespace

int main(int argc, char** argv)
{
    // -------------------------------------------------------------------------
    // Command-line input and configuration
    // -------------------------------------------------------------------------

    if (argc != 3) {
        std::cerr << "usage: " << argv[0]
                  << " <response.dat> <mb-spec.json>" << std::endl;
        return EXIT_FAILURE;
    }

    System sys;
    sys.SetUpFromJson("POWER", argv[2]);

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Error: cannot open file " << argv[1] << "\n";
        return EXIT_FAILURE;
    }


    // -------------------------------------------------------------------------
    // POWER spectrum settings and validation
    // -------------------------------------------------------------------------
    double tau = sys.GetTau();
    const double smoothing_fwhm_cm1 = sys.GetSmoothingFwhmCm1();
    const bool has_smoothing_fwhm_cm1 = sys.HasSmoothingFwhmCm1();
    double max_int = sys.GetIntMax();
    const double transform_time_ps = sys.GetSpectrumTransformTimePs();
    const bool has_transform_time_ps = sys.HasSpectrumTransformTimePs();

    {
        if (!std::is_same<decltype(tau), double>::value) {
            std::cerr << "could not convert '" << tau
            << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    }
    if (!has_transform_time_ps || transform_time_ps <= 0.0) {
        std::cerr << "POWER.spectrum_transform_time_ps must be positive"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (transform_time_ps > max_int) {
        std::cerr << "POWER.spectrum_transform_time_ps must not exceed POWER.int_max"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (tau < 0.0) {
        std::cerr << "POWER.tau must be nonnegative" << std::endl;
        return EXIT_FAILURE;
    }
    {
        if (!std::is_same<decltype(max_int), double>::value) {
            std::cerr << "could not convert '" << max_int
            << "' to real number" << std::endl;
            return EXIT_FAILURE;
        }
    }

    // -------------------------------------------------------------------------
    // Response input and correlation-time selection
    // -------------------------------------------------------------------------
    std::vector<double> time;
    std::vector<double> corr;

    size_t Ndata = 0;

    double t(0);
    double c(0);

    while (infile >> t >> c) {
        if (t >= max_int)
            break;
        ++Ndata;
        time.push_back(t);
        corr.push_back(c);
    }

    infile.close();

    std::cerr << "Ndata : " << Ndata << std::endl;


    if (Ndata < 2) {
        std::cerr << "response file must contain at least two points below POWER.int_max"
                  << std::endl;
        return EXIT_FAILURE;
    }
    const double dt_ps = time[1] - time[0];
    if (dt_ps <= 0.0) {
        std::cerr << "response time values must be strictly increasing"
                  << std::endl;
        return EXIT_FAILURE;
    }
    const double available_time_ps = Ndata*dt_ps;
    if (transform_time_ps > available_time_ps + 1.0e-10*dt_ps) {
        std::cerr << "POWER.spectrum_transform_time_ps exceeds the available "
                  << "response duration (" << available_time_ps << " ps)"
                  << std::endl;
        return EXIT_FAILURE;
    }
    const size_t n_transform = static_cast<size_t>(
        std::floor(transform_time_ps/dt_ps + 1.0e-10));
    if (n_transform < 2) {
        std::cerr << "POWER.spectrum_transform_time_ps must include at least "
                  << "two response intervals" << std::endl;
        return EXIT_FAILURE;
    }

    // -------------------------------------------------------------------------
    // Frequency grid and classical normalization
    // -------------------------------------------------------------------------
    double arg = 2.0 * M_PI / (double(2 * n_transform - 1));

    double dt = dt_ps * constants::ps_to_s;
    double dnu = 1.0 / (double(2 * n_transform - 1)) / dt;
    const double wavenumber_spacing_cm1 = dnu / conv_freq;

    if (corr[0] <= 0.0) {
        std::cerr << "POWER response must have a positive C(0) for normalization"
                  << std::endl;
        return EXIT_FAILURE;
    }
    const double velocity_normalization = 2.0 * conv_freq / corr[0];

    if (has_smoothing_fwhm_cm1 && smoothing_fwhm_cm1 < 0.0) {
        std::cerr << "POWER.smoothing_fwhm_cm1 must be nonnegative" << std::endl;
        return EXIT_FAILURE;
    }

    // -------------------------------------------------------------------------
    // Optional exponential time-domain window
    // -------------------------------------------------------------------------

    if(tau != 0)
        for (size_t t = 0; t < n_transform; ++t)
            corr[t] *= std::exp(-(t*dt)/(tau*constants::ps_to_s));

    // -------------------------------------------------------------------------
    // Cosine transform
    // -------------------------------------------------------------------------
    std::vector<double> freq;
    std::vector<double> abs;

    for (size_t inu = 0; inu < n_transform; ++inu){

        double nu = dnu * inu;

	double absorbance = 0.0;

	for (size_t itime = n_transform - 1; itime > 0; --itime){

	    double nu_time = (double(itime)) * (double(inu)) * arg;
	    absorbance += std::cos(nu_time) * corr[itime];

	}

	for (size_t itime = 0; itime < n_transform; ++itime){

	    double nu_time = (double(itime)) * (double(inu)) * arg;
	    absorbance += std::cos(nu_time) * corr[itime];

	}

        absorbance *= velocity_normalization;
        freq.push_back(nu/conv_freq);
	abs.push_back(absorbance);
	//std::cout << (nu/conv_freq) << '\t'<< absorbance << std::endl;

    }
  

    // -------------------------------------------------------------------------
    // Optional Gaussian spectral smoothing
    // -------------------------------------------------------------------------
    const bool use_smoothing =
        has_smoothing_fwhm_cm1 && smoothing_fwhm_cm1 > 0.0;
    std::vector<double> output = abs;
    if (use_smoothing) {
        if (!spectral::ApplyGaussianSmoothing(
                output, wavenumber_spacing_cm1, smoothing_fwhm_cm1)) {
            std::cerr << "POWER.smoothing_fwhm_cm1 is too broad for the available spectral grid" << std::endl;
            return EXIT_FAILURE;
        }
    }

    for (size_t n = 0; n < output.size(); ++n)
        std::cout << freq[n] << '\t' << output[n] << std::endl;
    return 0;
}    
