//////////////////////////////////////////////////////////////////////////
//  Calculates SFG response function for given polarization condition.  //
//  Options for polarizations combinations are: zzz, zxx, xzx, and xxz. //
//////////////////////////////////////////////////////////////////////////
#ifdef HAVE_CONFIG_H
#include "cli.h"
#include "config.h"
#endif // HAVE_CONFIG_H

#include <omp.h>

#include <cassert>
#include <cstdlib>

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include <cmath>
#include <algorithm>
#include <array>

#include "read-lammpstrj.h"
#include "system.h"
#include "pbc-tools.h"
#include "struct-sfg-pols.h"

#include "ps.h"
#include "avila.h"
#include "dip_poly2b5.h"
#include "alpha-el-poly2b3.h"
#include "alpha-ttm4.h"
#include "ttm4-dip.h"

////////////////////////////////////////////////////////////////////////////////

namespace {

////////////////////////////////////////////////////////////////////////////////

std::vector<opt::lammpstrj_frame> frames;
static x2o::avila pol_model;

std::vector<sfg::eprop> eprops;
std::vector<sfg::connect> connects;

//----------------------------------------------------------------------------//

inline double g_scale(const double* crd, const double zc1, const double zc2,
                      const int surface)
{
    if ((surface > 0 && crd[2] < 0.0) ||
        (surface < 0 && crd[2] > 0.0)) {
        return 0.0;
    }

    double z = std::fabs(crd[2]);
    if(z < zc1)
	return 0.0;
    else if(z > zc2)
	return copysign(1.0, crd[2]);
    else{
	const double c = std::cos(M_PI*(z-zc2)/(2*(zc1-zc2)));
	return copysign(c*c, crd[2]); // multiply switch by sign of z (crd[2])
    }
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

       dx -= cell[0]*std::floor(dx/cell[0] + 0.5);
       dy -= cell[4]*std::floor(dy/cell[4] + 0.5);
       dz -= cell[8]*std::floor(dz/cell[8] + 0.5);

       crd_out[0] = r1[0] - dx;
       crd_out[1] = r1[1] - dy;
       crd_out[2] = r1[2] - dz;

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

inline void image_monomer(const size_t imcon, const double* cell, const double* rcell,
       const double* r1, double* mon)
{
    std::fill(mon, mon + 9, 0.0);

    // take r1 oxygen as central atom, image its two hydrogens

    for(size_t i = 0; i < 3; ++i)
        mon[i] = r1[i];

    for(size_t a = 1; a < 3; ++a)
        image(imcon, cell, rcell, r1, r1 + 3*a, mon + 3*a);
}

inline void image_dimer(const size_t imcon, const double* cell,
                        const double* rcell, const double* r1,
                        const double* r2, double* dimer)
{
    image_monomer(imcon, cell, rcell, r1, dimer);
    image(imcon, cell, rcell, r1, r2, dimer + 9);
    for (size_t a = 1; a < 3; ++a) {
        image(imcon, cell, rcell, dimer + 9, r2 + 3*a,
              dimer + 9 + 3*a);
    }
}

inline double f_switch2(const double r, const double r2i, const double r2f)
{
    if (r > r2f) return 0.0;
    if (r > r2i) {
        const double x = (r - r2i)/(r2f - r2i);
        return 1.0 + x*x*(2.0*x - 3.0);
    }
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

    const auto exit_code = mbspec::check_arguments(argc, argv, {"trajectory.lammpstrj", "mb-spec.json"});

    if (exit_code.has_value()) {
        return *exit_code;
    }

    const auto trajectory_path = argv[1];
    const auto config_path = argv[2];

    System sys;
    sys.SetUpFromJson("SFG", config_path);

    const double zc1 = sys.GetZc1();

    const double zc2 = sys.GetZc2();

    if (!std::isfinite(zc1) || !std::isfinite(zc2) || zc1 < 0.0 || zc2 <= zc1) {
        std::cerr << "SFG.zc1 and SFG.zc2 must both be supplied as finite numbers "
                  << "in angstrom satisfying 0 <= zc1 < zc2; they define the "
                  << "interface weights. Set both keys in the SFG JSON section."
                  << std::endl;
        return EXIT_FAILURE;
    }
    
    const double int_max = sys.GetIntMax();

    if (int_max <= 0.0) {
        std::cerr << "SFG.int_max must be positive" << std::endl;
        return EXIT_FAILURE;
    }
    const double trajectory_duration_ps = sys.GetTrajectoryDurationPs();
    const size_t skip_initial_atoms = sys.GetSkipInitialAtoms();
    if (trajectory_duration_ps < 0.0) {
        std::cerr << "SFG.trajectory_duration_ps must be nonnegative" << std::endl;
        return EXIT_FAILURE;
    }
    if (trajectory_duration_ps > 0.0 && int_max > trajectory_duration_ps) {
        std::cerr << "SFG.int_max must not exceed SFG.trajectory_duration_ps"
                  << std::endl;
        return EXIT_FAILURE;
    }

    const double rOO_max = sys.GetCrossCorrelationCutoff();
    const size_t neighbor_sampling_stride = sys.GetCrossCorrelationSamplingStride();
    const bool include_autocorrelation = sys.GetIncludeAutocorrelation();
    const bool include_cross_correlation = sys.GetIncludeCrossCorrelation();
    const bool center_slab = sys.GetCenterSlab();
    const bool shift_slab_by_half_cell = sys.GetShiftSlabByHalfCell();
    const bool time_dependent_switch = sys.GetTimeDependentSwitch();
    const bool use_ps_dipole = sys.GetUsePsDipole();
    std::string surface = sys.GetSurface();
    std::transform(surface.begin(), surface.end(), surface.begin(), ::tolower);
    int surface_selection = 0;
    if (surface == "positive") {
        surface_selection = 1;
    } else if (surface == "negative") {
        surface_selection = -1;
    } else if (surface != "both") {
        std::cerr << "SFG.surface must be positive, negative, or both"
                  << std::endl;
        return EXIT_FAILURE;
    }

    if (rOO_max <= 0.0) {
        std::cerr << "SFG.cross_correlation_cutoff must be positive"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (neighbor_sampling_stride == 0) {
        std::cerr << "SFG.cross_correlation_sampling_stride must be positive"
                  << std::endl;
        return EXIT_FAILURE;
    }
    if (!include_autocorrelation && !include_cross_correlation) {
        std::cerr << "SFG must include autocorrelation and/or cross correlation"
                  << std::endl;
        return EXIT_FAILURE;
    }

    //----------------------------------------------------------------------------//
    // Polarization components and molecular-property models
    //----------------------------------------------------------------------------//

    std::string ijk = sys.GetIjk();
    std::transform(ijk.begin(), ijk.end(), ijk.begin(), ::tolower);

    int dipk(0);
    int dip2k(0);
    int polij(0);
    int pol2ij(0);

    std::cerr << "ijk = " << ijk << std::endl;

    // Define possible polarization combinations

    if (ijk.compare("zzz") == 0) {
        dipk = 2;		// z component of dipole vector
        dip2k = 2;		// z component of dipole vector
        polij = 8;		// zz element of polarizabilty tensor
        pol2ij = 8;		// zz element of polarizabilty tensor
    } else if (ijk.compare("zxx") == 0) { //|| (ijk.compare("zyy") == 0) {
        dipk = 0;		// x component of dipole vector
        dip2k = 1;		// y component of dipole vector
        polij = 6;		// zx element of polarizability tensor
        pol2ij = 7;		// zy element of polarizability tensor
    } else if (ijk.compare("xzx") == 0) { //|| (ijk.compare("yzy") == 0) {
        dipk = 0;		// x component of dipole vector
        dip2k = 1;		// y component of dipole vector
        polij = 2;		// xz element of polarizability tensor
        pol2ij = 5;		// yz element of polarizability tensor
    } else if (ijk.compare("xxz") == 0) { //|| (ijk.compare("yyz") == 0) {
        dipk = 2;		// z component of dipole vector
        dip2k = 2;		// z component of dipole vector
        polij = 0;		// xx element of polarizability tensor
        pol2ij = 4;		// yy element of polarizability tensor
    } else {
	std::cerr << "not a valid polarization combination" << std::endl;
	return EXIT_FAILURE;
    }

    std::string dipole_mb = sys.GetDipoleMb();
    std::string polarizability_mb = sys.GetPolarizabilityMb();
    std::transform(dipole_mb.begin(), dipole_mb.end(), dipole_mb.begin(), ::tolower);
    std::transform(polarizability_mb.begin(), polarizability_mb.end(),
                   polarizability_mb.begin(), ::tolower);
    const bool valid_dipole_mb = dipole_mb == "1b" || dipole_mb == "nb" ||
        dipole_mb == "1b+nb" || dipole_mb == "1b+2b" ||
        dipole_mb == "1b+2b+nb";
    const bool valid_polarizability_mb = polarizability_mb == "1b" ||
        polarizability_mb == "1b+2b";
    if (!valid_dipole_mb || !valid_polarizability_mb) {
        std::cerr << "SFG.dipole_mb must be '1b', 'nb', '1b+nb', '1b+2b', "
                  << "or '1b+2b+nb'; SFG.polarizability_mb must be '1b' "
                  << "or '1b+2b'" << std::endl;
        return EXIT_FAILURE;
    }
    const bool add_dipole_1b = dipole_mb != "nb";
    const bool add_dipole_nb = dipole_mb == "nb" || dipole_mb == "1b+nb" || dipole_mb == "1b+2b+nb";
    const bool add_dipole_2b = dipole_mb == "1b+2b" || dipole_mb == "1b+2b+nb";
    const bool replace_induced_2b = dipole_mb == "1b+2b+nb";
    const bool add_polarizability_2b = polarizability_mb == "1b+2b";
    const double r2i = sys.GetCutoff1();
    const double r2f = sys.GetCutoff2();
    if (r2i < 0.0 || r2f <= r2i) {
        std::cerr << "SFG.cutoff_1 and SFG.cutoff_2 must satisfy "
                  << "0 <= cutoff_1 < cutoff_2" << std::endl;
        return EXIT_FAILURE;
    }

    //----------------------------------------------------------------------------//
    // Periodic-boundary convention and trajectory sampling
    //----------------------------------------------------------------------------//

    size_t imcon(0);
    {
	imcon = sys.GetImcon();
    if (imcon == 0){
      imcon = 2;
    } else if (imcon < 0 || imcon > 3) {
      std::cerr << "Cannot calculate system of imcon = " << imcon << std::endl;
      return EXIT_FAILURE;
    } else if (imcon > 0 && imcon <= 3) {
      if (imcon == 1) {
        std::cerr << "imcon = " << imcon << " (cubic)" << std::endl;
      } else if (imcon == 2) {
        std::cerr << "imcon = " << imcon << " (orthorhombic)" << std::endl;
      } else if (imcon == 3) {
        std::cerr << "imcon = " << imcon << " (triclinic)" << std::endl;
      }
    } else {
      std::cerr << "imcon value of " << imcon << " not a valid integer" << std::endl;
      return EXIT_FAILURE;
    }
	}
    
    const double dump = sys.GetDumpFrequency();
    if (dump <= 0.0) {
        std::cerr << "SFG.dump_frequency must be positive (ps)" << std::endl;
        return EXIT_FAILURE;
    }

    //----------------------------------------------------------------------------//
    // Trajectory loading, validation, and optional slab centering
    //----------------------------------------------------------------------------//

    try {
        opt::read_lammpstrj(trajectory_path, frames, imcon, dump,
                            trajectory_duration_ps);
        std::cerr << " # Frames : " << frames.size() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    if (frames.size() < 2) {
        std::cerr << "SFG response requires at least two trajectory frames"
                  << std::endl;
        return EXIT_FAILURE;
    }

    //----------------------------------------------------------------------------//
    // Frame validation and molecule mass assignment
    //----------------------------------------------------------------------------//

    for (size_t n = 0; n < frames.size(); ++n) {
        if (skip_initial_atoms >= frames[n].natm) {
            std::cerr << "SFG.skip_initial_atoms must be smaller than the atom count"
                      << std::endl;
            return EXIT_FAILURE;
        }
        if ((frames[n].natm - skip_initial_atoms) % 3 != 0) {
            std::cerr << "SFG response requires contiguous O--H--H water triplets"
                      << std::endl;
            return EXIT_FAILURE;
        }
        for (size_t i = skip_initial_atoms; i < frames[n].natm; i += 3) {
            frames[n].mass[i] = 15.999;
            frames[n].mass[i + 1] = 1.008;
            frames[n].mass[i + 2] = 1.008;
        }
    }

    double dt = (frames[1].time - frames[0].time);
    int int_time_ind = (int) (int_max / dt);
    if (dt <= 0.0 || int_time_ind < 1) {
        std::cerr << "SFG.int_max must include at least one sampling interval"
                  << std::endl;
        return EXIT_FAILURE;
    }

    std::cerr << dt << "   " << int_time_ind << std::endl;

    //----------------------------------------------------------------------------//
    // Optional per-frame slab centering
    //----------------------------------------------------------------------------//

    if (center_slab) {
        #pragma omp parallel for
        for (size_t n = 0; n < frames.size(); ++n) {
            if (imcon < 3) {
                tools::center_slab_at_origin(
                    frames[n].natm - skip_initial_atoms, frames[n].cell,
                    frames[n].mass + skip_initial_atoms,
                    frames[n].xyz + 3*skip_initial_atoms,
                    shift_slab_by_half_cell);
            } else if (imcon == 3) {
                tools::center_slab_at_origin3(
                    frames[n].natm - skip_initial_atoms, imcon, frames[n].cell,
                    frames[n].rcell, frames[n].mass + skip_initial_atoms,
                    frames[n].xyz + 3*skip_initial_atoms,
                    shift_slab_by_half_cell);
            }
        }
    }

    //----------------------------------------------------------------------------//
    // Per-frame molecular properties and cross-correlation partners
    //----------------------------------------------------------------------------//

    for(size_t n = 0; n < frames.size(); ++n) {

        const size_t nmolec = (frames[n].natm - skip_initial_atoms)/3;
        const size_t coordinate_offset = 3*skip_initial_atoms;
        const size_t dipole_offset = 4*skip_initial_atoms;
        //------------------------------------------------------------------------//
        // Molecular one-body properties and interface weights
        //------------------------------------------------------------------------//

        { // calculate electrostatic properties
            sfg::eprop eprop_this_frame;
            eprop_this_frame.nmolec  = nmolec;
            eprop_this_frame.dip_k   = new double[eprop_this_frame.nmolec];
            eprop_this_frame.dip2_k  = new double[eprop_this_frame.nmolec];
            eprop_this_frame.pol_ij  = new double[eprop_this_frame.nmolec];
            eprop_this_frame.pol2_ij = new double[eprop_this_frame.nmolec];
            eprop_this_frame.gscale  = new double[eprop_this_frame.nmolec];
            std::vector<std::array<double, 9>> monomer_ttm_polarizabilities;
            if (add_polarizability_2b)
                monomer_ttm_polarizabilities.resize(nmolec);
 
            #pragma omp parallel for
            for(size_t i = 0; i < nmolec; ++i){
                const size_t i9 = coordinate_offset + 9*i;
 
                double mon[9];
                image_monomer(imcon, frames[n].cell, frames[n].rcell, frames[n].xyz + i9, mon);
 
                // k-component of the dipole
                double moldip_k = frames[n].dipmol[dipole_offset + 12*i + dipk];
                double moldip2_k = frames[n].dipmol[dipole_offset + 12*i + dip2k];
                if (use_ps_dipole) {
                    double ps_dipole[3];
                    ps::dipole(mon, ps_dipole);
                    moldip_k = ps_dipole[dipk];
                    moldip2_k = ps_dipole[dip2k];
                }
 
                eprop_this_frame.dip_k[i] =
                    (add_dipole_1b ? moldip_k : 0.0) +
                    (add_dipole_nb ? frames[n].dipind[dipole_offset + 12*i + dipk] : 0.0);
                eprop_this_frame.dip2_k[i] =
                    (add_dipole_1b ? moldip2_k : 0.0) +
                    (add_dipole_nb ? frames[n].dipind[dipole_offset + 12*i + dip2k] : 0.0);

                // ij-component of the polarizability
                double pol_mol[9];
      	        std::fill(pol_mol, pol_mol + 9, 0.0);
                pol_model.alpha(mon, pol_mol);
                eprop_this_frame.pol_ij[i] = pol_mol[polij];
                eprop_this_frame.pol2_ij[i] = pol_mol[pol2ij];
                if (add_polarizability_2b)
                    ttm::alpha_tot(1, mon, monomer_ttm_polarizabilities[i].data());
     
     		
                // scaling that accounts for the double-interface
                double g_i = g_scale(frames[n].xyz + i9, zc1, zc2,
                                     surface_selection);
                eprop_this_frame.gscale[i] = g_i;
            }

            //----------------------------------------------------------------------//
            // Switched and interface-partitioned two-body corrections
            //----------------------------------------------------------------------//

            if (add_dipole_2b || add_polarizability_2b) {
                std::vector<double> dipole_2b_k(nmolec, 0.0);
                std::vector<double> dipole_2b2_k(nmolec, 0.0);
                std::vector<double> polarizability_2b(nmolec, 0.0);
                std::vector<double> polarizability2_2b(nmolec, 0.0);

                #pragma omp parallel
                {
                    ttm::ttm4_dip induced_dipole_2b;
                    #pragma omp for
                    for (size_t i = 0; i < nmolec; ++i) {
                        const size_t i9 = coordinate_offset + 9*i;
                        for (size_t j = i + 1; j < nmolec; ++j) {
                            const size_t j9 = coordinate_offset + 9*j;
                            const double rOO = tools::distance3(
                                imcon, frames[n].cell, frames[n].rcell,
                                frames[n].xyz + i9, frames[n].xyz + j9);
                            if (rOO >= r2f) continue;

                            const double abs_g_i = std::abs(eprop_this_frame.gscale[i]);
                            const double abs_g_j = std::abs(eprop_this_frame.gscale[j]);
                            const double abs_g_sum = abs_g_i + abs_g_j;
                            if (abs_g_sum == 0.0) continue;

                            double dimer[18];
                            image_dimer(imcon, frames[n].cell, frames[n].rcell,
                                        frames[n].xyz + i9, frames[n].xyz + j9,
                                        dimer);
                            const double pair_weight = f_switch2(rOO, r2i, r2f);
                            const double partition_i = abs_g_i/abs_g_sum;
                            const double partition_j = abs_g_j/abs_g_sum;

                            if (add_dipole_2b) {
                                double correction[3];
                                h2o::dip_poly2b5::eval(dimer, dimer + 9,
                                                       correction);
                                if (replace_induced_2b) {
                                    double induced[3];
                                    induced_dipole_2b(2, dimer, induced);
                                    for (size_t k = 0; k < 3; ++k)
                                        correction[k] -= induced[k];
                                }
                                #pragma omp atomic
                                dipole_2b_k[i] += pair_weight*partition_i*correction[dipk];
                                #pragma omp atomic
                                dipole_2b_k[j] += pair_weight*partition_j*correction[dipk];
                                #pragma omp atomic
                                dipole_2b2_k[i] += pair_weight*partition_i*correction[dip2k];
                                #pragma omp atomic
                                dipole_2b2_k[j] += pair_weight*partition_j*correction[dip2k];
                            }
                            if (add_polarizability_2b) {
                                double correction[9];
                                h2o::alpha_el_poly2b3::eval(dimer, dimer + 9,
                                                             correction);
                                // The polynomial corrects the TTM pair increment.
                                // Switch and partition both terms together.
                                double dimer_ttm_polarizability[9];
                                ttm::alpha_tot(2, dimer, dimer_ttm_polarizability);
                                for (size_t k = 0; k < 9; ++k)
                                    correction[k] += dimer_ttm_polarizability[k]
                                        - monomer_ttm_polarizabilities[i][k]
                                        - monomer_ttm_polarizabilities[j][k];
                                #pragma omp atomic
                                polarizability_2b[i] += pair_weight*partition_i*correction[polij];
                                #pragma omp atomic
                                polarizability_2b[j] += pair_weight*partition_j*correction[polij];
                                #pragma omp atomic
                                polarizability2_2b[i] += pair_weight*partition_i*correction[pol2ij];
                                #pragma omp atomic
                                polarizability2_2b[j] += pair_weight*partition_j*correction[pol2ij];
                            }
                        }
                    }
                }

                #pragma omp parallel for
                for (size_t i = 0; i < nmolec; ++i) {
                    eprop_this_frame.dip_k[i] += dipole_2b_k[i];
                    eprop_this_frame.dip2_k[i] += dipole_2b2_k[i];
                    eprop_this_frame.pol_ij[i] += polarizability_2b[i];
                    eprop_this_frame.pol2_ij[i] += polarizability2_2b[i];
                }
            }
 
            eprops.push_back(eprop_this_frame);
 
        } // end electrostatics calculation
 
        //------------------------------------------------------------------------//
        // Cross-correlation neighbor lists
        //------------------------------------------------------------------------//

        { // determine connectivity
 
            sfg::connect connect_this_frame;
            connect_this_frame.nmolec = nmolec;
     
            for (size_t i = 0; i < nmolec; ++i) { 
                if (eprops[n].gscale[i] != 0.0) {
                    if (include_cross_correlation) {
                        const size_t i9 = 9*i;
                        std::vector<int> mol;
     
                        #pragma omp parallel for
                        for(size_t j = 0; j < nmolec; ++j){
                            double avg_rij(0);
                            size_t count(0);
                            if (j != i){
                            // Assume that int_time_ind is short enough to 
                            // prohibit diffusion between solvation shells
                            // and that the average_rij cannot be < rOO_max
                            // unless rij at the initial time is 
                            // less than 1.5*rOO_max
                            if (tools::distance3(imcon, frames[n].cell,
                                frames[n].rcell,
                                frames[n].xyz + coordinate_offset + i9,
                                frames[n].xyz + coordinate_offset + j*9)
                                < 2.0*rOO_max){
     
                                // calculate average rij distance over the next int_time_ind frames
                                for(size_t a = n;
                                    a < std::min(frames.size(), n + int_time_ind + 1);
                                    a += neighbor_sampling_stride){
         
                                    #pragma omp atomic
                                    avg_rij += tools::distance3(imcon,
                                                 frames[n].cell,
                                                 frames[n].rcell,
                                                 frames[n].xyz + coordinate_offset + i9,
                                                 frames[a].xyz + coordinate_offset + j*9);
         
                                    #pragma omp atomic
                                    ++count;
                                }
          
                                if ((avg_rij/count) < rOO_max) {
                                    #pragma omp critical
                                    {
                                    mol.push_back((int)j);
                                    }
                                }
     		            }
     		        }
                    }
     
                    connect_this_frame.partners.push_back(mol);
                    }

                    connect_this_frame.central.push_back(i);
     	        }
            }
     
      	    connects.push_back(connect_this_frame);
        }
    }
 
    if (eprops.size() != frames.size()) {
        std::cerr << "eprops.size() != frames.size()" << std::endl;
        exit(1);
    }
 
    std::cerr << "Deallocating frames" << std::endl;
    //frames.erase(frames.begin(), frames.end());
    std::vector<opt::lammpstrj_frame>().swap(frames);
 
    //----------------------------------------------------------------------------//
    // SFG response accumulation and output
    //----------------------------------------------------------------------------//

    std::vector<double> corr(int_time_ind, 0.0);
    std::vector<int> ntime(int_time_ind, 0);
 
    std::cerr << "Beginning correlation function" << std::endl;
 
    #pragma omp parallel for
    for(size_t m = 0; m < eprops.size(); ++m){ 	// dipole, time = 0
        for(size_t n = m; 						//polarizability, time = 't'
 	      n < std::min( eprops.size(), m + int_time_ind); ++ n){ 
 
            int time = n - m;
 
            double* dip_k_0  = eprops[m].dip_k;
            double* dip2_k_0  = eprops[m].dip2_k;
            double* pol_ij_t = eprops[n].pol_ij;
            double* pol2_ij_t = eprops[n].pol2_ij;

            for(size_t i = 0; i < connects[m].central.size(); ++i){
                const size_t i_ind = connects[m].central[i];
                const double g_i0 = eprops[m].gscale[i_ind];
 
            if (include_autocorrelation) {
                 const double g3 = time_dependent_switch
                     ? g_i0*eprops[n].gscale[i_ind]*eprops[n].gscale[i_ind]
                     : g_i0*g_i0*g_i0;
 
                 #pragma omp atomic 
                 corr[time] += g3*0.5*(dip_k_0[i_ind]*pol_ij_t[i_ind]+dip2_k_0[i_ind]*pol2_ij_t[i_ind]);
            }
 
            if (include_cross_correlation) {
                for(size_t j = 0; j < connects[m].partners[i].size(); ++j){
                    const size_t j_ind = connects[m].partners[i][j];
                    const double g_j0 = eprops[n].gscale[j_ind];
                    {
                    #pragma omp atomic
                    corr[time] += g_i0*g_j0*g_j0*0.5*(dip_k_0[i_ind]*pol_ij_t[j_ind]+dip2_k_0[i_ind]*pol2_ij_t[j_ind]);
                    }
                 }
            }
            }
 
 	   #pragma omp atomic
 	   ++ntime[time];
        }
    }
 
    //----------------------------------------------------------------------------//
    // Time-origin normalization and response output
    //----------------------------------------------------------------------------//

    std::cout << std::setprecision(10) << std::scientific;
 
    for(int i = 0; i < int_time_ind; ++i)
        std::cout << std::setw(20) << dt*i 
        << std::setw(20) << corr[i]/((double)ntime[i])
        << std::setw(20) << corr[i]
        << std::setw(20) << ((double)ntime[i])
        << std::endl;
}
