#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <sstream>

#include "cli.h"
#include "system.h"

namespace {

void do_autocorr(const char* filename, const double int_max,
                 const double trajectory_duration_ps, const bool depol,
                 const bool subtract_mean_polarizability,
                 std::vector<double>& autocorrelation, double& dt)
{
    //----------------------------------------------------------------------------//
    // Polarizability-tensor input
    //----------------------------------------------------------------------------//

    std::vector<double> bxx;
    std::vector<double> bxy;
    std::vector<double> bxz;
    std::vector<double> byy;
    std::vector<double> byz;
    std::vector<double> bzz;
    std::vector<double> aiso;

    std::ifstream ifs(filename);

    if (!ifs) {
        std::cerr << " ** Error ** : could not open '"
                  << filename << "' for reading"
                  << std::endl;
        exit(1);
    }

//    std::string line;
//    getline(ifs,line);
//    std::istringstream iss(line);
//    iss >> vol;

    double previous_time = 0.0;
    bool has_previous_time = false;
    std::string line;
    size_t line_number = 0;
    while (std::getline(ifs, line)) {
        ++line_number;
        std::istringstream line_stream(line);
        double t, my_axx, my_ayy, my_azz, my_axy, my_axz, my_ayz;
        std::string extra_field;
        if (!(line_stream >> t >> my_axx >> my_axy >> my_axz
                          >> my_ayy >> my_ayz >> my_azz) ||
            line_stream >> extra_field) {
            std::cerr << " ** Error ** : '" << filename << "' line "
                      << line_number << " must contain exactly seven numeric columns: "
                      << "time alpha_xx alpha_xy alpha_xz alpha_yy alpha_yz alpha_zz"
                      << std::endl;
            exit(1);
        }
        if (trajectory_duration_ps > 0.0 && t > trajectory_duration_ps + 1.0e-12)
            break;

        if (has_previous_time) {
            const double interval = t - previous_time;
            if (interval <= 0.0) {
                std::cerr << " ** Error ** : time values in '" << filename
                          << "' must be strictly increasing" << std::endl;
                exit(1);
            }
            if (bxx.size() == 1) {
                dt = interval;
            } else if (std::abs(interval - dt) > 1.0e-10 * dt) {
                std::cerr << " ** Error ** : time values in '" << filename
                          << "' must be evenly spaced" << std::endl;
                exit(1);
            }
        }
        previous_time = t;
        has_previous_time = true;

        const double a_iso = (my_axx + my_ayy + my_azz)/3.0;
        bxx.push_back(my_axx - a_iso);
        byy.push_back(my_ayy - a_iso);
        bzz.push_back(my_azz - a_iso);
        bxy.push_back(my_axy);
        bxz.push_back(my_axz);
        byz.push_back(my_ayz);
        aiso.push_back(a_iso);
    }
    int this_file_nsamples = (int) bxx.size();
    if (this_file_nsamples < 2 || dt <= 0.0) {
        std::cerr << " ** Error ** : '" << filename
                  << "' must contain at least two time-ordered tensor samples"
                  << std::endl;
        exit(1);
    }

    if (subtract_mean_polarizability) {
        double mean_aiso = 0.0;
        double mean_bxx = 0.0;
        double mean_bxy = 0.0;
        double mean_bxz = 0.0;
        double mean_byy = 0.0;
        double mean_byz = 0.0;
        double mean_bzz = 0.0;

        for (int i = 0; i < this_file_nsamples; ++i) {
            mean_aiso += aiso[i];
            mean_bxx += bxx[i];
            mean_bxy += bxy[i];
            mean_bxz += bxz[i];
            mean_byy += byy[i];
            mean_byz += byz[i];
            mean_bzz += bzz[i];
        }

        const double inverse_samples = 1.0/this_file_nsamples;
        mean_aiso *= inverse_samples;
        mean_bxx *= inverse_samples;
        mean_bxy *= inverse_samples;
        mean_bxz *= inverse_samples;
        mean_byy *= inverse_samples;
        mean_byz *= inverse_samples;
        mean_bzz *= inverse_samples;

        for (int i = 0; i < this_file_nsamples; ++i) {
            aiso[i] -= mean_aiso;
            bxx[i] -= mean_bxx;
            bxy[i] -= mean_bxy;
            bxz[i] -= mean_bxz;
            byy[i] -= mean_byy;
            byz[i] -= mean_byz;
            bzz[i] -= mean_bzz;
        }
    }

    const int int_time_ind = static_cast<int>(int_max / dt);
    if (int_time_ind < 1) {
        std::cerr << " ** Error ** : RAMAN.int_max must be at least one "
                  << "sampling interval (" << dt << ")" << std::endl;
        exit(1);
    }
    const int n_lags = std::min(this_file_nsamples, int_time_ind);

    //----------------------------------------------------------------------------//
    // Isotropic or depolarized autocorrelation
    //----------------------------------------------------------------------------//

    std::vector<double> corr(n_lags, 0.0);
    std::vector<double> ntime(n_lags, 0.0);

    for(int i = 0; i < this_file_nsamples; ++i){
	for(int j = i; j < std::min(this_file_nsamples, i + n_lags); ++j){
	    int time = (j-i);

        if (depol) {

            // xx = xx*xx + xy*yx + xz*zx
            double xx = bxx[i]*bxx[j] + bxy[i]*bxy[j] + bxz[i]*bxz[j];
            
            // yy = yx*xy + yy*yy + yz*zy
            double yy = bxy[i]*bxy[j] + byy[i]*byy[j] + byz[i]*byz[j];

            // zz = zx*xz + zy*yz + zz*zz
            double zz = bxz[i]*bxz[j] + byz[i]*byz[j] + bzz[i]*bzz[j];

            corr[time] += xx + yy + zz;
        }
        else {
            corr[time] += aiso[i]*aiso[j];
        }
	    ++ntime[time];
	}
    }

    for(size_t i = 0; i < corr.size(); ++i)
	corr[i] /= ntime[i];

    std::cerr << "* loaded " << n_lags
              << " samples from '" << filename << "'"
              << std::endl;

    autocorrelation.swap(corr);
}

} // namespace

int main(int argc, char **argv)
{
    //----------------------------------------------------------------------------//
    // Program setup and JSON configuration
    //----------------------------------------------------------------------------//

    const auto exit_code = mbspec::check_arguments(argc, argv, {"polarizability_file", "mb-spec.json"});

    if (exit_code.has_value()) {
        return *exit_code;
    }

    const auto polarizability_path = argv[1];
    const auto config_path = argv[2];

    System sys;
    sys.SetUpFromJson("RAMAN", config_path);

    double int_max(0);
    {
        int_max = sys.GetIntMax();
        if (int_max <= 0.0) {
            std::cerr << "RAMAN.int_max must be positive" << std::endl;
            return EXIT_FAILURE;
        }
    }
    const double trajectory_duration_ps = sys.GetTrajectoryDurationPs();
    if (trajectory_duration_ps < 0.0 ||
        (trajectory_duration_ps > 0.0 && int_max > trajectory_duration_ps)) {
        std::cerr << "RAMAN requires 0 < int_max <= trajectory_duration_ps "
                  << "when trajectory_duration_ps is nonzero" << std::endl;
        return EXIT_FAILURE;
    }

    // depol = true
    bool depol = sys.GetRamanResponse();
    const bool subtract_mean_polarizability =
        sys.GetSubtractMeanPolarizability();

    //----------------------------------------------------------------------------//
    // Response calculation and output
    //----------------------------------------------------------------------------//

    std::vector<double> autocorrelation;
    double dt = 0.0;
    do_autocorr(polarizability_path, int_max, trajectory_duration_ps, depol,
                subtract_mean_polarizability,
                autocorrelation, dt);

    std::cout << std::setprecision(8) << std::scientific;

    for (size_t k = 0; k < autocorrelation.size(); ++k) {
        std::cout << std::setw(18) << dt*k
                  << std::setw(18) << autocorrelation[k]
                  << '\n';
    }
}
