#include <cmath>
#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <sstream>

#include "cli.h"
#include "system.h"

namespace {

const double debye = 3.33564e-30; // Coulomb*m

//----------------------------------------------------------------------------//
// Dipole autocorrelation for one trajectory
//----------------------------------------------------------------------------//

void do_autocorr(const char* filename, const double int_max,
                 const double trajectory_duration_ps,
                 const bool subtract_mean_dipole,
                 std::vector<double>& autocorrelation, double& dt)
{
    //----------------------------------------------------------------------------//
    // Dipole input and time-grid validation
    //----------------------------------------------------------------------------//

    std::vector<double> v3x;
    std::vector<double> v3y;
    std::vector<double> v3z;
    double previous_time = 0.0;
    double this_dt = 0.0;

    std::ifstream ifs(filename);
    if (!ifs) {
        std::cerr << " ** Error ** : could not open '"
                  << filename << "' for reading"
                  << std::endl;
        exit(1);
    }

    // Input rows are exactly: time [ps], mu_x [D], mu_y [D], mu_z [D].
    std::string line;
    size_t line_number = 0;
    while (std::getline(ifs, line)) {
        ++line_number;
        std::istringstream line_stream(line);
        double t, dx, dy, dz;
        std::string extra_field;
        if (!(line_stream >> t >> dx >> dy >> dz) || line_stream >> extra_field) {
            std::cerr << " ** Error ** : '" << filename << "' line " << line_number
                      << " must contain exactly four numeric columns: time mu_x mu_y mu_z"
                      << std::endl;
            exit(1);
        }
        if (trajectory_duration_ps > 0.0 && t > trajectory_duration_ps + 1.0e-12)
            break;

        if (!v3x.empty()) {
            const double interval = t - previous_time;
            if (interval <= 0.0) {
                std::cerr << " ** Error ** : time values in '" << filename
                          << "' must be strictly increasing" << std::endl;
                exit(1);
            }
            if (v3x.size() == 1) {
                this_dt = interval;
            } else if (std::abs(interval - this_dt) > 1.0e-10 * this_dt) {
                std::cerr << " ** Error ** : time values in '" << filename
                          << "' must be evenly spaced" << std::endl;
                exit(1);
            }
        }
        previous_time = t;
        // Convert the total dipole from Debye to SI before correlating.
        v3x.push_back(dx*debye);
        v3y.push_back(dy*debye);
        v3z.push_back(dz*debye);
    }

    int this_file_nsamples = (int) v3x.size();
    if (this_file_nsamples < 2 || this_dt <= 0.0) {
        std::cerr << " ** Error ** : '" << filename
                  << "' must contain at least two time-ordered dipole samples"
                  << std::endl;
        exit(1);
    }

    dt = this_dt;

    if (subtract_mean_dipole) {
        const double mean_x = std::accumulate(v3x.begin(), v3x.end(), 0.0)
                            / this_file_nsamples;
        const double mean_y = std::accumulate(v3y.begin(), v3y.end(), 0.0)
                            / this_file_nsamples;
        const double mean_z = std::accumulate(v3z.begin(), v3z.end(), 0.0)
                            / this_file_nsamples;

        for (int i = 0; i < this_file_nsamples; ++i) {
            v3x[i] -= mean_x;
            v3y[i] -= mean_y;
            v3z[i] -= mean_z;
        }
    }

    //----------------------------------------------------------------------------//
    // Autocorrelation evaluation
    //----------------------------------------------------------------------------//

    // int_max is an exclusive upper limit: the final lag is below int_max.
    const int int_time_ind = static_cast<int>(int_max / dt);
    if (int_time_ind < 1) {
        std::cerr << " ** Error ** : INFRARED.int_max must be at least one "
                  << "sampling interval (" << dt << ")" << std::endl;
        exit(1);
    }

    const int n_lags = std::min(this_file_nsamples, int_time_ind);

    std::vector<double> corr_x(n_lags, 0.0);
    std::vector<double> corr_y(n_lags, 0.0);
    std::vector<double> corr_z(n_lags, 0.0);
    std::vector<double> ntime(n_lags, 0.0);

    // Estimate C(k dt) = <mu(i) . mu(i+k)> using every valid time origin.
    // Each lag is normalized by its own number of available time origins.
    for(int i = 0; i < this_file_nsamples; ++i){
	for(int j = i; j < std::min(this_file_nsamples, i + n_lags); ++j){
	    int time = (j-i);

	    corr_x[time] += v3x[i]*v3x[j];
	    corr_y[time] += v3y[i]*v3y[j];
	    corr_z[time] += v3z[i]*v3z[j];
	    ++ntime[time];
	}
    }

    for(int i = 0; i < n_lags; ++i){
	corr_x[i] /= ntime[i];
    	corr_y[i] /= ntime[i];
    	corr_z[i] /= ntime[i];
    }

    std::cerr << "* loaded " << n_lags
              << " samples from '" << filename << "'"
              << std::endl;

    autocorrelation.resize(n_lags);
    for (int k = 0; k < n_lags; ++k) {
        autocorrelation[k] = corr_x[k] + corr_y[k] + corr_z[k];
    }
}

} // namespace

int main(int argc, char **argv)
{

    //----------------------------------------------------------------------------//
    // Program setup and response output
    //----------------------------------------------------------------------------//

    const auto exit_code = mbspec::check_arguments(argc, argv, {"dipole_file", "mb-spec.json"});

    if (exit_code.has_value()) {
        return *exit_code;
    }

    const auto dipole_path = argv[1];
    const auto config_path = argv[2];

    System sys;
    sys.SetUpFromJson("INFRARED", config_path);

    double int_max(0);
    {
    	int_max = sys.GetIntMax();
	if (int_max <= 0.0) {
        std::cerr << "INFRARED.int_max must be positive" << std::endl;
    	    return EXIT_FAILURE;
    	}
    }
    const double trajectory_duration_ps = sys.GetTrajectoryDurationPs();
    if (trajectory_duration_ps < 0.0 ||
        (trajectory_duration_ps > 0.0 && int_max > trajectory_duration_ps)) {
        std::cerr << "INFRARED requires 0 < int_max <= trajectory_duration_ps "
                  << "when trajectory_duration_ps is nonzero" << std::endl;
        return EXIT_FAILURE;
    }

    const bool subtract_mean_dipole = sys.GetSubtractMeanDipole();
    std::vector<double> autocorrelation;
    double dt = 0.0;

    do_autocorr(dipole_path, int_max, trajectory_duration_ps,
                subtract_mean_dipole, autocorrelation, dt);

    // Output columns are lag time [ps] and C(t) [(C m)^2].
    std::cout << std::setprecision(8) << std::scientific;

    for (size_t k = 0; k < autocorrelation.size(); ++k) {
        std::cout << std::setw(18) << dt*k
                  << std::setw(18) << autocorrelation[k]
                  << '\n';
    }
}
