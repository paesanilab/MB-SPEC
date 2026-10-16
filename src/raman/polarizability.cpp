#ifdef HAVE_CONFIG_H
#include "config.h"
#endif // HAVE_CONFIG_H

#include <omp.h>

#include <cmath>
#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <array>
#include <atomic>

#include "read-lammpstrj.h"
#include "constants.h"
#include "system.h"

#include "avila.h"

#include "alpha-el-poly2b3.h"

#include "alpha-ttm4.h"

#include "ttm4-periodic-es.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////


static x2o::avila pol_model;

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

////////////////////////////////////////////////////////////////////////////////

} // namespace

////////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{

    //----------------------------------------------------------------------------//
    // Program setup and JSON configuration
    //----------------------------------------------------------------------------//

    if (argc < 3) {
        std::cerr << "usage: " << argv[0] << " <file_name> <mb-spec.json>" << std::endl;
        return EXIT_FAILURE;
    }

    System sys;
    sys.SetUpFromJson("RAMAN", argv[argc - 1]);

    std::vector<double> box(9);

    box = sys.GetBoxVector();

    const double r2i = sys.GetCutoff1();
    const double r2f = sys.GetCutoff2();

    if (r2i < 0.0 || r2f <= r2i) {
        std::cerr << "RAMAN.cutoff_1 and RAMAN.cutoff_2 must satisfy "
                  << "0 <= cutoff_1 < cutoff_2" << std::endl;
        return EXIT_FAILURE;
    }


    std::string fname = argv[1];
    std::ifstream ifs_cmd(fname);
    if (!ifs_cmd) {
        std::cerr << "could not open '" << fname << "' for reading" << std::endl;
        return EXIT_FAILURE;
    }
    const size_t pos = fname.rfind('.');
    std::string format = (pos == std::string::npos) ? "" : fname.substr(pos + 1);
    std::transform(format.begin(), format.end(), format.begin(), ::tolower);


    // Open files for Writing

    std::cout.setf(std::ios_base::showpoint);
    std::cout.precision(9);

    std::cerr.setf(std::ios_base::scientific);
    std::cerr.precision(12);

    size_t lineno(0);
    size_t frameno(0);
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
        } else if (imcon == 3) {
            std::cerr << "imcon = " << imcon << " (monoclinic)" << std::endl;
        }
    } else if (imcon == 3) {
        std::cerr << "imcon = " << imcon << " (monoclinic)" << std::endl;
    }

    double dump(0);
    dump = sys.GetDumpFrequency();
    const size_t skip_initial_atoms = sys.GetSkipInitialAtoms();
    const double trajectory_duration_ps = sys.GetTrajectoryDurationPs();
    if (trajectory_duration_ps < 0.0) {
        std::cerr << "RAMAN.trajectory_duration_ps must be nonnegative" << std::endl;
        return EXIT_FAILURE;
    }

    // Toggles names: 1B, NB (1BNB), MB to lowercase
    std::string mb = sys.GetPolarizabilityMb();
    std::transform(mb.begin(), mb.end(), mb.begin(), ::tolower);

    const bool add_1b = mb != "nb";
    const bool add_nb = mb == "nb" || mb == "1b+nb" ||
                        mb == "1b+2b+nb";
    const bool add_2b = mb == "1b+2b" || mb == "1b+2b+nb";
    const bool add_pair_induction = mb == "1b+2b";
    const bool valid_mode = mb == "1b" || mb == "nb" ||
                            mb == "1b+nb" || mb == "1b+2b" ||
                            mb == "1b+2b+nb";
    if (!valid_mode) {
        std::cerr << "RAMAN.polarizability_mb must be '1b', 'nb', '1b+nb', '1b+2b', or "
                  << "'1b+2b+nb'" << std::endl;
        return EXIT_FAILURE;
    }
    if (add_nb && imcon == 3) {
        std::cerr << "RAMAN.imcon=3 (triclinic) is unsupported for polarizability_mb='"
                  << mb << "': the internal NB model uses only orthorhombic cell lengths. "
                  << "Use an orthorhombic cell, select '1b' or '1b+2b', or supply "
                  << "externally calculated polarizabilities to response."
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (add_nb) {
        std::cerr << "WARNING: RAMAN NB-containing modes use the internal "
                  << "TTM4 periodic model without an Ewald treatment. "
                  << "They are not validated for periodic bulk simulations."
                  << std::endl;
    }

    const bool write_molecular_properties = sys.GetWriteMolecularProperties();
    if (write_molecular_properties && mb != "1b") {
        std::cerr << "RAMAN.write_molecular_properties requires "
                  << "RAMAN.polarizability_mb='1b'"
                  << std::endl;
        return EXIT_FAILURE;
    }

    std::ofstream molecular_output;
    if (write_molecular_properties) {
        molecular_output.open("molecular_polarizabilities.dat");
        if (!molecular_output) {
            std::cerr << "could not open 'molecular_polarizabilities.dat' "
                      << "for writing" << std::endl;
            return EXIT_FAILURE;
        }
        molecular_output.setf(std::ios_base::showpoint);
        molecular_output.precision(9);
        molecular_output << "# molecule time_ps alpha_xx_A3 alpha_xy_A3 "
                         << "alpha_xz_A3 alpha_yx_A3 alpha_yy_A3 "
                         << "alpha_yz_A3 alpha_zx_A3 alpha_zy_A3 "
                         << "alpha_zz_A3" << std::endl;
    }


    //----------------------------------------------------------------------------//
    // Frame-batch reading and OpenMP tensor evaluation
    //----------------------------------------------------------------------------//

    while (!ifs_cmd.eof()) {

        std::vector<opt::lammpstrj_frame> frames;
        const int batch_size = std::max(1, omp_get_max_threads());

        while ((int) frames.size() < batch_size && !ifs_cmd.eof()) {
            if (trajectory_duration_ps > 0.0 &&
                dump*frameno > trajectory_duration_ps + 1.0e-12) break;
            opt::lammpstrj_frame this_frame;
            try {
                if(format.compare("xyz") == 0) {
                    opt::read_xyz_frame(lineno, imcon, ifs_cmd, this_frame, dump * frameno, box);
                } else {
                    opt::read_lammpstrj_frame(lineno, imcon, ifs_cmd, this_frame,
                                              dump * frameno, false);
                }
            } catch (const std::exception& error) {
                std::cerr << "Raman trajectory input error in '" << argv[1]
                          << "': " << error.what() << std::endl;
                return EXIT_FAILURE;
            }
            if(format != "xyz" && ifs_cmd.eof()) break;
            if(this_frame.natm == 0) break;
            frames.push_back(this_frame);
            ++frameno;
        }

        const int n_frames = frames.size();
        if (n_frames == 0) break;

        std::vector<std::array<double, 9>> frame_polarizabilities(n_frames);
        std::vector<std::vector<std::array<double, 9>>> molecular_polarizabilities;
        if (write_molecular_properties) {
            for (int frame_index = 0; frame_index < n_frames; ++frame_index) {
                if (skip_initial_atoms >= frames[frame_index].natm) {
                    std::cerr << "RAMAN.skip_initial_atoms must be smaller than the atom count"
                              << std::endl;
                    return EXIT_FAILURE;
                }
            }
            molecular_polarizabilities.resize(n_frames);
            for (int frame_index = 0; frame_index < n_frames; ++frame_index) {
                molecular_polarizabilities[frame_index].resize(
                    (frames[frame_index].natm - skip_initial_atoms)/3);
            }
        }
        std::atomic<bool> invalid_water(false);
        std::atomic<bool> invalid_selection(false);

#pragma omp parallel for
        for (int frame_index = 0; frame_index < n_frames; ++frame_index) {
            const opt::lammpstrj_frame& frame = frames[frame_index];
            if (skip_initial_atoms >= frame.natm) {
                invalid_selection = true;
                continue;
            }
            if ((frame.natm - skip_initial_atoms) % 3 != 0) {
                invalid_water = true;
                continue;
            }

            const size_t nmolec = (frame.natm - skip_initial_atoms)/3;
            const size_t coordinate_offset = 3*skip_initial_atoms;
            double total_pol[9] = {0.0};
            double pol1b[9] = {0.0};
            double pol2b_corr[9] = {0.0};
            double polMb[9] = {0.0};
            double dimer[18];

            std::vector<std::array<double, 9>> monomer_ttm_polarizabilities;
            if (add_pair_induction) {
                monomer_ttm_polarizabilities.resize(nmolec);
                for (size_t i = 0; i < nmolec; ++i) {
                    const size_t i9 = coordinate_offset + 9*i;
                    double monomer[9];
                    std::copy(frame.xyz + i9, frame.xyz + i9 + 3, monomer);
                    for (size_t a = 1; a < 3; ++a)
                        image(imcon, frame.cell, frame.rcell, frame.xyz + i9,
                              frame.xyz + i9 + 3*a, monomer + 3*a);
                    ttm::alpha_tot(1, monomer, monomer_ttm_polarizabilities[i].data());
                }
            }

            for (size_t i = 0; i < nmolec; ++i) {
                size_t i9 = coordinate_offset + 9*i;

                /////
                // 1B contribution: avila
                /////

            // One-body Avila and removable TTM4 monomer contributions.
                {
      	            double monomer[9];
      	            std::copy(frame.xyz + i9, frame.xyz + i9 + 3, monomer);
      	            for(size_t a = 1; a < 3; ++a)
      	                image(imcon, frame.cell, frame.rcell, frame.xyz + i9,
      		                  frame.xyz + i9 + 3*a, monomer + 3*a);

      	            // Avila 1B polarizability
      	            double pol_mol[9];
      	            std::fill(pol_mol, pol_mol + 9, 0.0);
      	            pol_model.alpha(monomer, pol_mol);

                    if (write_molecular_properties) {
                        std::copy(pol_mol, pol_mol + 9,
                                  molecular_polarizabilities[frame_index][i].begin());
                    }
      	            for (size_t a = 0; a < 9; ++a)
      	                pol1b[a] += pol_mol[a];

      	            // TTM4 1B polarizability
                    if (add_nb) {
                        double pol_ind_1b[9];
      	                std::fill(pol_ind_1b, pol_ind_1b + 9, 0.0);
      	                ttm::alpha_tot(1, monomer, pol_ind_1b);
      	                for(size_t a = 0; a < 9; ++a)
      	                    polMb[a] -= pol_ind_1b[a];
                    }
                }
                
                /////
                // 1B+2B switches both the polynomial and isolated-pair induction.
                // The full NB mode already includes induction, so adds only the polynomial.
                /////

                if (add_2b) {
                    for (size_t j = i+1; j < nmolec; ++j) {
                        const size_t j9 = coordinate_offset + 9*j;

    	                double rOO = distance3(imcon, frame.cell, frame.rcell, 
      		 	                               frame.xyz + i9, frame.xyz + j9);

      	                if (rOO < r2f) {

                            const double s = f_switch2(rOO, r2i, r2f);
      	                    image_dimer(imcon, frame.cell, frame.rcell, 
      			                frame.xyz + i9, frame.xyz + j9, dimer);

      	                    double corr[9];
      	                    std::fill(corr, corr + 9, 0.0);
      	                    h2o::alpha_el_poly2b3::eval(dimer, dimer + 9, corr);

                            if (add_pair_induction) {
                                double dimer_ttm_polarizability[9];
                                ttm::alpha_tot(2, dimer, dimer_ttm_polarizability);
                                for (size_t b = 0; b < 9; ++b)
                                    corr[b] += dimer_ttm_polarizability[b]
                                        - monomer_ttm_polarizabilities[i][b]
                                        - monomer_ttm_polarizabilities[j][b];
                            }

      	                    for (size_t b = 0; b < 9; ++b) {
      		                pol2b_corr[b] += s*corr[b];
      	                    }
      	                }    
                    }
                }
            }

            /////
            // MB contribution: alpha^{ttm4,NB}(1, \dots, N) 
            //                  - \sum_{i  }^N \alpha^{ttm4,1B}(i)
            //
            //    Note: alpha^{ttm4,NB} contains a monomer polarizability
            //          Therefore, the 1B polarizability
            //          must be removed to obtain the many-body polarizability
            //
            //          The 1B polarizability has already been
            //          removed above in the polMb term
            /////
 	
            if (add_nb) {

	        double ind[9];
	        std::fill(ind, ind + 9, 0.0);

	        ttm::ttm4_periodic_es ind_model;
	        ind_model(nmolec, frame.xyz + coordinate_offset, frame.cell);
	        ind_model.ptensor(nmolec, ind);

	        for(size_t a = 0; a < 9; ++a)
		    polMb[a] += ind[a];
            }

	    for(size_t a = 0; a < 9; ++a) {
	        total_pol[a] = pol2b_corr[a] + polMb[a];
                if (add_1b) {
                    total_pol[a] += pol1b[a];
                }
            }

            std::copy(total_pol, total_pol + 9,
                      frame_polarizabilities[frame_index].begin());
        }

        //----------------------------------------------------------------------------//
        // Ordered tensor output
        //----------------------------------------------------------------------------//

        if (invalid_water) {
            std::cerr << "MB-alpha can only be calculated for water trajectories"
                      << std::endl;
            return EXIT_FAILURE;
        }
        if (invalid_selection) {
            std::cerr << "RAMAN.skip_initial_atoms must be smaller than the atom count"
                      << std::endl;
            return EXIT_FAILURE;
        }

        for (int i = 0; i < n_frames; ++i) {
            const std::array<double, 9>& pol = frame_polarizabilities[i];
            std::cout << frames[i].time << ' '
                      << std::setw(14) << pol[0] << ' '
                      << std::setw(14) << pol[1] << ' '
                      << std::setw(14) << pol[2] << ' '
                      << std::setw(14) << pol[4] << ' '
                      << std::setw(14) << pol[5] << ' '
                      << std::setw(14) << pol[8] << std::endl;
            if (write_molecular_properties) {
                for (size_t molecule = 0;
                     molecule < molecular_polarizabilities[i].size(); ++molecule) {
                    const std::array<double, 9>& pol_mol =
                        molecular_polarizabilities[i][molecule];
                    molecular_output << molecule + 1 << ' ' << frames[i].time;
                    for (size_t component = 0; component < 9; ++component) {
                        molecular_output << ' ' << pol_mol[component];
                    }
                    molecular_output << std::endl;
                }
            }
        }
    }

    return 0;
}
