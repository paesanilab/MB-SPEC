#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <omp.h>

#include "system.h"
#include "read-lammpstrj.h"

namespace {

// -------------------------------------------------------------------------
// Native LAMMPS velocity-frame reader
// -------------------------------------------------------------------------

bool read_lammpstrj_velocity_frame(std::istream& input, size_t skip_initial_atoms,
                                   double total_velocity[3])
{
    std::string line;
    if (!std::getline(input, line))
        return false;
    if (line != "ITEM: TIMESTEP")
        throw std::runtime_error("expected 'ITEM: TIMESTEP' in LAMMPS trajectory");

    if (!std::getline(input, line))
        throw std::runtime_error("missing LAMMPS timestep value");
    if (!std::getline(input, line) || line != "ITEM: NUMBER OF ATOMS")
        throw std::runtime_error("expected 'ITEM: NUMBER OF ATOMS' in LAMMPS trajectory");
    if (!std::getline(input, line))
        throw std::runtime_error("missing LAMMPS atom count");
    const size_t natoms = static_cast<size_t>(std::stoull(line));
    if (skip_initial_atoms >= natoms)
        throw std::runtime_error("POWER.skip_initial_atoms must be smaller than the atom count");

    if (!std::getline(input, line) || line.rfind("ITEM: BOX BOUNDS", 0) != 0)
        throw std::runtime_error("expected 'ITEM: BOX BOUNDS' in LAMMPS trajectory");
    for (size_t i = 0; i < 3; ++i) {
        if (!std::getline(input, line))
            throw std::runtime_error("incomplete LAMMPS box bounds");
    }

    // ---------------------------------------------------------------------
    // Velocity-column discovery from the LAMMPS atom header
    // ---------------------------------------------------------------------

    if (!std::getline(input, line) || line.rfind("ITEM: ATOMS", 0) != 0)
        throw std::runtime_error("expected 'ITEM: ATOMS' in LAMMPS trajectory");
    std::istringstream header(line.substr(std::string("ITEM: ATOMS").size()));
    std::vector<std::string> columns;
    for (std::string column; header >> column; )
        columns.push_back(column);

    const auto find_column = [&columns](const std::string& name) {
        const auto found = std::find(columns.begin(), columns.end(), name);
        return found == columns.end() ? columns.size() :
            static_cast<size_t>(found - columns.begin());
    };
    size_t vx_column = find_column("vx");
    size_t vy_column = find_column("vy");
    size_t vz_column = find_column("vz");
    if (vx_column == columns.size() || vy_column == columns.size() ||
        vz_column == columns.size()) {
        throw std::runtime_error(
            "LAMMPS trajectory must provide vx, vy, and vz velocity columns");
    }

    // ---------------------------------------------------------------------
    // Selected-atom total velocity for one frame
    // ---------------------------------------------------------------------

    std::fill(total_velocity, total_velocity + 3, 0.0);
    const size_t last_velocity_column = std::max(vx_column,
        std::max(vy_column, vz_column));
    for (size_t atom = 0; atom < natoms; ++atom) {
        if (!std::getline(input, line))
            throw std::runtime_error("incomplete LAMMPS atom section");
        std::istringstream atom_line(line);
        std::vector<double> values;
        for (double value = 0.0; atom_line >> value; )
            values.push_back(value);
        if (values.size() <= last_velocity_column)
            throw std::runtime_error("LAMMPS atom record does not contain all velocity columns");
        if (atom >= skip_initial_atoms) {
            total_velocity[0] += values[vx_column];
            total_velocity[1] += values[vy_column];
            total_velocity[2] += values[vz_column];
        }
    }
    return true;
}

// -------------------------------------------------------------------------
// Collective velocity autocorrelation
// -------------------------------------------------------------------------

void do_autocorr(const std::vector<double>& data, double int_max)
{
    const size_t stride = 4;
    const size_t nframes = data.size()/stride;
    if (nframes < 2)
        throw std::runtime_error("not enough frames for a velocity autocorrelation");

    const double dt = data[stride] - data[0];
    if (dt <= 0.0)
        throw std::runtime_error("trajectory times must be strictly increasing");
    size_t int_time_ind = static_cast<size_t>(std::floor(int_max/dt));
    if (int_time_ind < 1)
        throw std::runtime_error("POWER.int_max must include at least one sampling interval");
    int_time_ind = std::min(int_time_ind, nframes);

    // ---------------------------------------------------------------------
    // Correlation storage and thread-local accumulation
    // ---------------------------------------------------------------------

    std::vector<double> corr_vel(3*int_time_ind, 0.0);
    std::vector<size_t> nav(int_time_ind, 0);

#pragma omp parallel
    {
        std::vector<double> corr_private(3*int_time_ind, 0.0);
        std::vector<size_t> nav_private(int_time_ind, 0);

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < nframes; ++i) {
            const size_t last = std::min(nframes, i + int_time_ind);
            for (size_t j = i; j < last; ++j) {
                const size_t lag = j - i;
                ++nav_private[lag];
                for (size_t component = 0; component < 3; ++component) {
                    corr_private[3*lag + component] +=
                        data[stride*i + 1 + component]*
                        data[stride*j + 1 + component];
                }
            }
        }

#pragma omp critical
        {
            for (size_t lag = 0; lag < int_time_ind; ++lag) {
                nav[lag] += nav_private[lag];
                for (size_t component = 0; component < 3; ++component)
                    corr_vel[3*lag + component] +=
                        corr_private[3*lag + component];
            }
        }
    }

    // ---------------------------------------------------------------------
    // Time-origin normalization and response output
    // ---------------------------------------------------------------------

    for (size_t lag = 0; lag < int_time_ind; ++lag) {
        double corr_total = 0.0;
        for (size_t component = 0; component < 3; ++component)
            corr_total += corr_vel[3*lag + component]/nav[lag];
        std::cout << dt*lag << ' ' << corr_total << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    // -------------------------------------------------------------------------
    // Command-line input and POWER configuration
    // -------------------------------------------------------------------------

    if (argc != 3) {
        std::cerr << "usage: " << argv[0]
                  << " <velocity_trajectory> <mb-spec.json>" << std::endl;
        return EXIT_FAILURE;
    }

    System sys;
    sys.SetUpFromJson("POWER", argv[2]);
    const double dump_frequency = sys.GetDumpFrequency();
    const double int_max = sys.GetIntMax();
    const size_t skip_initial_atoms = sys.GetSkipInitialAtoms();
    const double trajectory_duration_ps = sys.GetTrajectoryDurationPs();
    if (dump_frequency <= 0.0 || int_max <= 0.0) {
        std::cerr << "POWER.dump_frequency and POWER.int_max must be positive"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (trajectory_duration_ps < 0.0) {
        std::cerr << "POWER.trajectory_duration_ps must be nonnegative"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (trajectory_duration_ps > 0.0 && int_max > trajectory_duration_ps) {
        std::cerr << "POWER.int_max must not exceed POWER.trajectory_duration_ps"
                  << std::endl;
        return EXIT_FAILURE;
    }

    // -------------------------------------------------------------------------
    // Trajectory input and per-frame total velocity
    // -------------------------------------------------------------------------

    const std::string filename = argv[1];
    const size_t suffix = filename.find_last_of('.');
    const std::string format = suffix == std::string::npos ? "" :
        filename.substr(suffix + 1);
    std::ifstream input(filename);
    if (!input) {
        std::cerr << "could not open '" << filename << "' for reading"
                  << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<double> total_velocities;
    try {
        if (format == "xyz") {
            size_t lineno = 0;
            size_t imcon = sys.GetImcon();
            std::vector<double> box = sys.GetBoxVector();
            for (size_t frame_number = 0; ; ++frame_number) {
                const double time = dump_frequency*frame_number;
                if (trajectory_duration_ps > 0.0 &&
                    time > trajectory_duration_ps + 1.0e-12) {
                    break;
                }
                opt::lammpstrj_frame frame;
                opt::read_xyz_frame(lineno, imcon, input, frame, time, box);
                if (frame.natm == 0)
                    break;
                if (skip_initial_atoms >= frame.natm)
                    throw std::runtime_error(
                        "POWER.skip_initial_atoms must be smaller than the atom count");
                double total_velocity[3] = {0.0, 0.0, 0.0};
                for (size_t atom = skip_initial_atoms; atom < frame.natm; ++atom) {
                    for (size_t component = 0; component < 3; ++component)
                        total_velocity[component] += frame.xyz[3*atom + component];
                }
                total_velocities.push_back(time);
                total_velocities.insert(total_velocities.end(), total_velocity,
                                        total_velocity + 3);
            }
        } else {
            for (size_t frame_number = 0; ; ++frame_number) {
                const double time = dump_frequency*frame_number;
                if (trajectory_duration_ps > 0.0 &&
                    time > trajectory_duration_ps + 1.0e-12) {
                    break;
                }
                double total_velocity[3];
                if (!read_lammpstrj_velocity_frame(input, skip_initial_atoms,
                                                   total_velocity)) {
                    break;
                }
                total_velocities.push_back(time);
                total_velocities.insert(total_velocities.end(), total_velocity,
                                        total_velocity + 3);
            }
        }
    } catch (const std::exception& error) {
        std::cerr << "power trajectory input error: " << error.what()
                  << std::endl;
        return EXIT_FAILURE;
    }

    // -------------------------------------------------------------------------
    // Optional removal of the trajectory-average total velocity
    // -------------------------------------------------------------------------

    const size_t nframes = total_velocities.size()/4;
    if (nframes < 2) {
        std::cerr << "trajectory must contain at least two usable frames"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (sys.GetSubtractMeanVelocity()) {
        double mean_total_velocity[3] = {0.0, 0.0, 0.0};
        for (size_t frame = 0; frame < nframes; ++frame) {
            for (size_t component = 0; component < 3; ++component)
                mean_total_velocity[component] +=
                    total_velocities[4*frame + 1 + component];
        }
        for (size_t component = 0; component < 3; ++component)
            mean_total_velocity[component] /= nframes;
        for (size_t frame = 0; frame < nframes; ++frame) {
            for (size_t component = 0; component < 3; ++component)
                total_velocities[4*frame + 1 + component] -=
                    mean_total_velocity[component];
        }
    }

    // -------------------------------------------------------------------------
    // Collective velocity autocorrelation
    // -------------------------------------------------------------------------

    try {
        do_autocorr(total_velocities, int_max);
    } catch (const std::exception& error) {
        std::cerr << "power response error: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
