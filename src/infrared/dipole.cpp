#ifdef HAVE_CONFIG_H
#include "cli.h"
#include "config.h"
#endif // HAVE_CONFIG_H

#include <omp.h>

//----------------------------------------------------------------------------//
// Standard-library dependencies
//----------------------------------------------------------------------------//

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <algorithm>
#include <array>
#include <atomic>

//----------------------------------------------------------------------------//
// MB-SPEC dependencies and water dipole models
//----------------------------------------------------------------------------//

#include "read-lammpstrj.h"
#include "constants.h"

#include "ttm4-dip.h"
#include "ltp2011.h"
#include "system.h"

#include "dip_poly2b5.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

//----------------------------------------------------------------------------//
// Water dipole models
//----------------------------------------------------------------------------//

static h2o::ltp2011 dip1b_model;

//----------------------------------------------------------------------------//
// Periodic-boundary geometry helpers
//----------------------------------------------------------------------------//

inline double distance3(size_t imcon, const double* cell, const double* rcell,
		       const double* r1, const double* r2)
{
    assert(imcon == 0 || imcon == 1 || imcon == 2 || imcon == 3);

    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    if (imcon == 1 || imcon == 2) {

        dx -= cell[0]*std::floor(dx/cell[0] + 0.5);
        dy -= cell[4]*std::floor(dy/cell[4] + 0.5);
        dz -= cell[8]*std::floor(dz/cell[8] + 0.5);

    } else if (imcon == 3) {

        const double ssx = rcell[0]*dx + rcell[3]*dy + rcell[6]*dz;
        const double ssy = rcell[1]*dx + rcell[4]*dy + rcell[7]*dz;
        const double ssz = rcell[2]*dx + rcell[5]*dy + rcell[8]*dz;

        const double xss = ssx - std::floor(ssx + 0.5);
        const double yss = ssy - std::floor(ssy + 0.5);
        const double zss = ssz - std::floor(ssz + 0.5);

        dx = cell[0]*xss + cell[3]*yss + cell[6]*zss;
        dy = cell[1]*xss + cell[4]*yss + cell[7]*zss;
        dz = cell[2]*xss + cell[5]*yss + cell[8]*zss;

    }

    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

//----------------------------------------------------------------------------//

inline void image(const size_t imcon, const double* cell, const double* rcell,
           const double* r1, const double* r2, double* crd_out)
{
    assert(imcon == 0 || imcon == 1 || imcon == 2 || imcon == 3);

    double dx = r1[0] - r2[0];
    double dy = r1[1] - r2[1];
    double dz = r1[2] - r2[2];

    if (imcon == 1 || imcon == 2) {

       double facx = std::floor(dx/cell[0] + 0.5);
       double facy = std::floor(dy/cell[4] + 0.5);
       double facz = std::floor(dz/cell[8] + 0.5);

       crd_out[0] = r2[0] + cell[0]*facx;
       crd_out[1] = r2[1] + cell[4]*facy;
       crd_out[2] = r2[2] + cell[8]*facz;

    } else if (imcon == 3) {

       const double ssx = rcell[0]*dx + rcell[3]*dy + rcell[6]*dz;
       const double ssy = rcell[1]*dx + rcell[4]*dy + rcell[7]*dz;
       const double ssz = rcell[2]*dx + rcell[5]*dy + rcell[8]*dz;

       const double xss = ssx - std::floor(ssx + 0.5);
       const double yss = ssy - std::floor(ssy + 0.5);
       const double zss = ssz - std::floor(ssz + 0.5);

       dx = cell[0]*xss + cell[3]*yss + cell[6]*zss;
       dy = cell[1]*xss + cell[4]*yss + cell[7]*zss;
       dz = cell[2]*xss + cell[5]*yss + cell[8]*zss;

       crd_out[0] = r1[0] - dx;
       crd_out[1] = r1[1] - dy;
       crd_out[2] = r1[2] - dz;

    }
}

//----------------------------------------------------------------------------//

inline void image_dimer(const size_t imcon, const double* cell, 
                const double* rcell, const double* r1, const double* r2,
		      	double* dimer)
{
    std::fill(dimer, dimer + 18, 0.0);

    // take r1 oxygen as central atom, image its two hydrogens
     
    for(size_t i = 0; i < 3; ++i)
	    dimer[i] = r1[i];
    
    for(size_t a = 1; a < 3; ++a)
	    image(imcon, cell, rcell, r1, r1 + 3*a, dimer + 3*a);

    // Now image the oxygen of the second molecule

    image(imcon, cell, rcell, r1, r2, dimer + 9);

    // Finally, reattach hydrogens of molec two to oxygen of molec 2

    for(size_t a = 1; a < 3; ++a)
	    image(imcon, cell, rcell, dimer + 9, r2 + 3*a, dimer + 9 + 3*a);
}

//----------------------------------------------------------------------------//

double f_switch2(const double& r, const double r2i, const double r2f)
{
    if (r > r2f)
	    return 0.0;
    else if (r > r2i) {
	    const double x = (r - r2i)/(r2f - r2i);
	    return 1 + x*x*(2*x - 3);
    } else
	    return 1.0;
}

//----------------------------------------------------------------------------//
// End of local helpers
//----------------------------------------------------------------------------//

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

  //----------------------------------------------------------------------------//
  // Program setup and JSON input
  //----------------------------------------------------------------------------//

  const auto exit_code = mbspec::check_arguments(argc, argv, {"trajectory", "mb-spec.json"});

  if (exit_code.has_value()) {
    return *exit_code;
  }

  const std::string trajectory_path = argv[1];
  const auto config_path = argv[2];

  System sys;
  sys.SetUpFromJson("INFRARED", config_path);

  std::vector<double> box(9);

  box = sys.GetBoxVector();

  const double r2i = sys.GetCutoff1();
  const double r2f = sys.GetCutoff2();

  //----------------------------------------------------------------------------//
  // Trajectory and input-format setup
  //----------------------------------------------------------------------------//

  std::cout.setf(std::ios_base::showpoint);
  std::cout.precision(9);

  size_t lineno(0);
  size_t frameno(0);

  std::ifstream ifs_cmd(trajectory_path);
  if (!ifs_cmd) {
    std::cerr << "could not open '" << trajectory_path << "' for reading" << std::endl;
    return EXIT_FAILURE;
  }

  const size_t pos = trajectory_path.rfind('.');
  std::string format = (pos == std::string::npos) ? "" : trajectory_path.substr(pos + 1);
  std::transform(format.begin(), format.end(), format.begin(), ::tolower);

  size_t imcon(0);
  {
    imcon = sys.GetImcon();
  }

  if (imcon == 0){
    imcon = 2;
  } else if (imcon > 3) {
    std::cerr << "Cannot calculate system of imcon = " << imcon << std::endl;
    return EXIT_FAILURE;
  } else if (imcon > 0 && imcon < 3) {
    if (imcon == 1) {
      std::cerr << "imcon = " << imcon << " (cubic)" << std::endl;
    } else if (imcon == 2) {
      std::cerr << "imcon = " << imcon << " (orthorhombic)" << std::endl;
    }
  } else if (imcon == 3) {
    std::cerr << "imcon = " << imcon << " (monoclinic)" << std::endl;
  }

  const bool dms = sys.GetDMS();

  double dump(0);
  dump = sys.GetDumpFrequency();
  const size_t skip_initial_atoms = sys.GetSkipInitialAtoms();
  const double trajectory_duration_ps = sys.GetTrajectoryDurationPs();
  if (trajectory_duration_ps < 0.0) {
    std::cerr << "INFRARED.trajectory_duration_ps must be nonnegative" << std::endl;
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------//
  // IR dipole-model option validation
  //----------------------------------------------------------------------------//

  if (dms) {
    std::cerr << "Using LTP 1B DMS" << std::endl;
  } else {
    std::cerr << "Reading from " << trajectory_path << std::endl;
  }


  // Toggles names: 1B, NB (1BNB), MB to lowercase
  std::string mb = sys.GetDipoleMb();
  std::transform(mb.begin(), mb.end(), mb.begin(), ::tolower);

  const bool add_poly_2b = mb == "1b+2b" || mb == "1b+2b+nb";
  const bool replace_induced_2b = !dms && mb == "1b+2b+nb";
  const bool add_permanent = !dms && mb != "nb";
  const bool add_induced = !dms &&
      (mb == "nb" || mb == "1b+nb" || mb == "1b+2b+nb");

  const bool valid_mode = dms
      ? (mb == "1b" || mb == "1b+2b")
      : (mb == "1b" || mb == "nb" || mb == "1b+nb" ||
         mb == "1b+2b" || mb == "1b+2b+nb");
  if (!valid_mode) {
    if (dms) {
      std::cerr << "With INFRARED.DMS=true, INFRARED.dipole_mb must be '1b' or '1b+2b'"
                << std::endl;
    } else {
      std::cerr << "With INFRARED.DMS=false, INFRARED.dipole_mb must be '1b', 'nb', "
                << "'1b+nb', '1b+2b', or '1b+2b+nb'" << std::endl;
    }
    return EXIT_FAILURE;
  }

  const bool write_molecular_properties = sys.GetWriteMolecularProperties();
  if (write_molecular_properties && (!dms || mb != "1b")) {
    std::cerr << "INFRARED.write_molecular_properties requires DMS=true and "
              << "INFRARED.dipole_mb='1b'" << std::endl;
    return EXIT_FAILURE;
  }

  std::ofstream molecular_output;
  if (write_molecular_properties) {
    molecular_output.open("molecular_dipoles.dat");
    if (!molecular_output) {
      std::cerr << "could not open 'molecular_dipoles.dat' for writing"
                << std::endl;
      return EXIT_FAILURE;
    }
    molecular_output.setf(std::ios_base::showpoint);
    molecular_output.precision(9);
    molecular_output << "# molecule time_ps mu_x_D mu_y_D mu_z_D" << std::endl;
  }

  // duh can't read dipoles if there r no dipoles
  if (!dms && format.compare("xyz") == 0) {
    std::cerr << ".lammpstrj format required to read dipoles" << std::endl;
    return EXIT_FAILURE;
  }

  //----------------------------------------------------------------------------//
  // Frame-batch reading and OpenMP dipole evaluation
  //----------------------------------------------------------------------------//

  while(!ifs_cmd.eof() ) {
    std::vector<opt::lammpstrj_frame> frames;

    const int batch_size = std::max(1, omp_get_max_threads());
    while((int) frames.size() < batch_size && !ifs_cmd.eof()){
      if (trajectory_duration_ps > 0.0 &&
          dump*frameno > trajectory_duration_ps + 1.0e-12) break;
      opt::lammpstrj_frame this_frame;
      try {
        if(format.compare("xyz") == 0) {
          opt::read_xyz_frame(lineno, imcon, ifs_cmd, this_frame, dump * frameno, box);
        } else {
          opt::read_lammpstrj_frame(lineno, imcon, ifs_cmd, this_frame,
                                    dump * frameno, !dms);
        }
      } catch (const std::exception& error) {
        std::cerr << "IR trajectory input error in '" << trajectory_path
                  << "': " << error.what() << std::endl;
        return EXIT_FAILURE;
      }

      frameno++;

      if (format != "xyz" && ifs_cmd.eof()) break;
      if (this_frame.natm == 0) break;

      frames.push_back(this_frame);

    }

    int n_frames = frames.size();

    if (n_frames == 0) break;


    std::vector<std::array<double, 3>> frame_dipoles(n_frames);
    std::vector<double> frame_times(n_frames);
    std::vector<std::vector<std::array<double, 3>>> molecular_dipoles;
    if (write_molecular_properties) {
      for (int frame_index = 0; frame_index < n_frames; ++frame_index) {
        if (skip_initial_atoms >= frames[frame_index].natm) {
          std::cerr << "INFRARED.skip_initial_atoms must be smaller than the atom count"
                    << std::endl;
          return EXIT_FAILURE;
        }
      }
      molecular_dipoles.resize(n_frames);
      for (int frame_index = 0; frame_index < n_frames; ++frame_index) {
        molecular_dipoles[frame_index].resize(
            (frames[frame_index].natm - skip_initial_atoms)/3);
      }
    }
    std::atomic<bool> invalid_water(false);
    std::atomic<bool> invalid_selection(false);

#pragma omp parallel
    {
      // ttm4_dip has mutable work storage, so each OpenMP thread needs its own instance.
      ttm::ttm4_dip dip2b_ind_model;

#pragma omp for
      for (int frame_index = 0; frame_index < n_frames; ++frame_index) {
      const opt::lammpstrj_frame& frame = frames[frame_index];
      if (skip_initial_atoms >= frame.natm) {
        invalid_selection = true;
        continue;
      }
      const size_t selected_natm = frame.natm - skip_initial_atoms;
      const size_t nmolec = selected_natm/3;
      if ((add_poly_2b || dms) &&
          selected_natm % 3 != 0) {
        invalid_water = true;
        continue;
      }

      double total_dip[3] = {0.0, 0.0, 0.0};
      double dip2b_corr[3] = {0.0, 0.0, 0.0};
      double dip2b_ind[3] = {0.0, 0.0, 0.0};

      // Explicit short-range MB-mu two-body polynomial correction.
      if (add_poly_2b) {
        double dimer[18];
        for (size_t i = 0; i < nmolec; ++i) {
          for (size_t j = i + 1; j < nmolec; ++j) {
            const size_t i9 = 3*skip_initial_atoms + 9*i;
            const size_t j9 = 3*skip_initial_atoms + 9*j;
            const double rOO = distance3(imcon, frame.cell, frame.rcell,
                                         frame.xyz + i9, frame.xyz + j9);

            if (rOO < r2f) {
              const double s = f_switch2(rOO, r2i, r2f);
              image_dimer(imcon, frame.cell, frame.rcell,
                          frame.xyz + i9, frame.xyz + j9, dimer);

              double tmp_dip2b[3];
              h2o::dip_poly2b5::eval(dimer, dimer + 9, tmp_dip2b);
              for (size_t b = 0; b < 3; ++b) {
                dip2b_corr[b] += s*tmp_dip2b[b];
              }

              if (replace_induced_2b) {
                // Replace the induced short-range 2B term already present in NB.
                dip2b_ind_model(2, dimer, tmp_dip2b);
                for (size_t b = 0; b < 3; ++b) {
                  dip2b_ind[b] += s*tmp_dip2b[b];
                }
              }
            }
          }
        }
      }

      for(size_t i = skip_initial_atoms; i < frame.natm; ++i) {
        const size_t i3 = 3*i;
        const size_t i4 = 4*i;

        if (!dms) {
          for (size_t j = 0; j < 3; ++j) {
            if (add_permanent) {
              total_dip[j] += frame.dipmol[i4 + j];
            }
            if (add_induced) {
              total_dip[j] += frame.dipind[i4 + j];
            }
          }
        } else if ((i - skip_initial_atoms) % 3 == 0) {
          double monomer[9];
          double dip_mol[3];
          std::copy(frame.xyz + i3, frame.xyz + i3 + 3, monomer);
          for(size_t a = 1; a < 3; ++a) {
            image(imcon, frame.cell, frame.rcell, frame.xyz + i3,
                  frame.xyz + i3 + 3*a, monomer + 3*a);
          }
          dip1b_model.cartesian(1, monomer, dip_mol);
          if (write_molecular_properties) {
            molecular_dipoles[frame_index][(i - skip_initial_atoms)/3] = {
                dip_mol[0], dip_mol[1], dip_mol[2]};
          }
          for (size_t component = 0; component < 3; ++component) {
            total_dip[component] += dip_mol[component];
          }
        }
      }

      if (add_poly_2b) {
        for (size_t component = 0; component < 3; ++component) {
          total_dip[component] += dip2b_corr[component];
          if (replace_induced_2b) {
            // 1B+2B+NB: add the polynomial and remove the replaced induction.
            total_dip[component] -= dip2b_ind[component];
          }
        }
      }

      frame_times[frame_index] = frame.time;
      frame_dipoles[frame_index] = {total_dip[0], total_dip[1], total_dip[2]};
      }
    }

    if (invalid_water) {
      std::cerr << "MB-mu calculations can only be used with water trajectories"
                << std::endl;
      return EXIT_FAILURE;
    }
    if (invalid_selection) {
      std::cerr << "INFRARED.skip_initial_atoms must be smaller than the atom count"
                << std::endl;
      return EXIT_FAILURE;
    }

    //----------------------------------------------------------------------------//
    // Ordered dipole-series output
    //----------------------------------------------------------------------------//

    for (int i = 0; i < n_frames; ++i) {
      std::cout << frame_times[i] << ' '
                << std::setw(14) << frame_dipoles[i][0] << ' '
                << std::setw(14) << frame_dipoles[i][1] << ' '
                << std::setw(14) << frame_dipoles[i][2] << std::endl;
      if (write_molecular_properties) {
        for (size_t molecule = 0; molecule < molecular_dipoles[i].size();
             ++molecule) {
          const std::array<double, 3>& dipole = molecular_dipoles[i][molecule];
          molecular_output << molecule + 1 << ' ' << frame_times[i] << ' '
                           << dipole[0] << ' ' << dipole[1] << ' '
                           << dipole[2] << std::endl;
        }
      }
    }

  }

}
