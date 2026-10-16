#include <cmath>
#include <vector>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>

#include <cassert>
#include <cstdlib>
#include <sstream>
#include "system.h"
#include "constants.h"
#include "spectral_smoothing.h"

#define RAMAN yes
#define DESYM yes

namespace {

const double pi2 = 2.0*M_PI;

const double speed_of_light_cm_ps =
    constants::c0*constants::cm_per_m*constants::ps_to_s; // cm/ps

//const double epsilon_vacuum = 8.8541878e-12;
const double light_velocity = constants::c0; // m/s

//const double T = 298.15; // K

} // namespace

int main(int argc, char** argv)
{
    //----------------------------------------------------------------------------//
    // Program setup and JSON configuration
    //----------------------------------------------------------------------------//

    if (argc != 3) {
        std::cerr << "usage: " << argv[0]
                  << " <response.dat> <mb-spec.json>" << std::endl;
        return EXIT_FAILURE;
    }

    System sys;
    sys.SetUpFromJson("RAMAN", argv[2]);

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Error: cannot open file " << argv[1] << "\n";
        return EXIT_FAILURE;
    }


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
            std::cerr << "RAMAN.temperature must be positive (K)" << std::endl;
            return EXIT_FAILURE;
        }
    }
    if (!has_transform_time_ps || transform_time_ps <= 0.0 ||
        transform_time_ps > max_int) {
        std::cerr << "RAMAN.spectrum_transform_time_ps must satisfy "
                  << "0 < spectrum_transform_time_ps <= int_max" << std::endl;
        return EXIT_FAILURE;
    }
    {
        if (tau < 0.0) {
            std::cerr << "RAMAN.tau must be nonnegative (ps)" << std::endl;
            return EXIT_FAILURE;
        }
    }
    {
        if (volume <= 0.0) {
            std::cerr << "RAMAN.box must define a positive simulation-cell volume (A^3)"
                      << std::endl;
            return EXIT_FAILURE;
        }
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
    // Response-function input and optional time cutoff
    //----------------------------------------------------------------------------//

    std::vector<double> time;
    std::vector<double> corr;

    double t(0);
    double correlation_value(0);

    while (infile >> t >> correlation_value) {
        if (t >= max_int)
            break;
        time.push_back(t);
        corr.push_back(correlation_value);
    }

	infile.close();

    const size_t n_response = time.size();
    if (n_response < 4) {
        std::cerr << "response file must contain at least four time points" << std::endl;
        return EXIT_FAILURE;
    }

    const double dt = time[1] - time[0]; // ps
    if (dt <= 0.0) {
        std::cerr << "response time values must be strictly increasing" << std::endl;
        return EXIT_FAILURE;
    }
    for (size_t i = 2; i < n_response; ++i) {
        if (time[i] <= time[i - 1]) {
            std::cerr << "response time values must be strictly increasing" << std::endl;
            return EXIT_FAILURE;
        }
    }

    const double available_time_ps = n_response*dt;
    if (transform_time_ps > available_time_ps + 1.0e-10*dt) {
        std::cerr << "RAMAN.spectrum_transform_time_ps exceeds the available "
                  << "response duration (" << available_time_ps << " ps)" << std::endl;
        return EXIT_FAILURE;
    }
    const size_t nt = static_cast<size_t>(
        std::floor(transform_time_ps/dt + 1.0e-10));
    if (nt < 2) {
        std::cerr << "RAMAN.spectrum_transform_time_ps must include at least "
                  << "two response intervals" << std::endl;
        return EXIT_FAILURE;
    }
	std::cerr << "Ndata : " << nt << std::endl;

    const double wavenumber_factor = 1.0/speed_of_light_cm_ps;

    if (has_smoothing_fwhm_cm1 && smoothing_fwhm_cm1 < 0.0) {
        std::cerr << "RAMAN.smoothing_fwhm_cm1 must be nonnegative" << std::endl;
        return EXIT_FAILURE;
    }

    //----------------------------------------------------------------------------//
    // Time-domain window
    //----------------------------------------------------------------------------//

    if(tau != 0)
	for (size_t t = 0; t < nt; ++t)
	    corr[t] *= std::exp(-(t*dt)/tau);
    //----------------------------------------------------------------------------//
    // Cosine transform and Raman frequency factors
    //----------------------------------------------------------------------------//

    const double dnu = 1.0/(double(2*nt - 1)*dt); // ps^-1

    std::vector<double> freq;
    std::vector<double> intensity;

    for (size_t k = 0; k < nt; ++k) {
	const double nu = k*dnu; // ps-1

	// assume autocorrelation has faded to 0 by upper integration limit
	double accum = 0.5*corr[0];
	for (size_t t = 1; t < (nt - 1); ++t)
	    accum += std::cos(pi2*nu*t*dt)*corr[t];

	// changing even integral from \int_{-\infty}^{\infty} 
	// to 2*\int_0^{\infty}
	accum = 2.0*accum*dt; //ps

	freq.push_back(wavenumber_factor*nu);
	intensity.push_back(accum);
	//intensity.push_back(accum*factor_quant/100.0);

    }

#ifdef RAMAN
    for (size_t k = 0; k < freq.size(); ++k)
	//intensity[k] *= 1.0/std::pow((20500.0-freq[k]),4);
	intensity[k] *= freq[k]*light_velocity/std::pow((20500.0-freq[k]),4);
#endif

#ifdef DESYM
    for (size_t k = 0; k < freq.size(); ++k)
	intensity[k] *= std::tanh( 0.5*(1.0/T/constants::kB)*constants::hbar_Js
				   * pi2*freq[k]*constants::cm_per_m*light_velocity );
#endif

    if (sys.GetBoseEinsteinCorrection()) {
        // R_BE = R/omega; the finite-correlation limit at omega=0 is zero.
        intensity[0] = 0.0;
        for (size_t k = 1; k < freq.size(); ++k) {
            const double omega = pi2*freq[k]*constants::cm_per_m*light_velocity; // rad/s
            intensity[k] /= omega;
        }
    }

    //----------------------------------------------------------------------------//
    // Optional Gaussian spectral smoothing
    //----------------------------------------------------------------------------//

    const bool use_smoothing =
        has_smoothing_fwhm_cm1 && smoothing_fwhm_cm1 > 0.0;
    if (use_smoothing) {
        const double wavenumber_spacing_cm1 = wavenumber_factor*dnu;
        if (!spectral::ApplyGaussianSmoothing(
                intensity, wavenumber_spacing_cm1, smoothing_fwhm_cm1)) {
            std::cerr << "RAMAN.smoothing_fwhm_cm1 is too broad for "
                      << "the available spectral grid" << std::endl;
            return EXIT_FAILURE;
        }
    }

    //----------------------------------------------------------------------------//
    // Spectrum output
    //----------------------------------------------------------------------------//

    for (size_t k = 0; k < freq.size(); ++k)
	std::cout << std::setw(15) << freq[k]
		  << std::setw(15) << intensity[k] << std::endl;

}
