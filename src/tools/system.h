#ifndef SYSTEM_H
#define SYSTEM_H

#include <vector>
#include <tuple>
#include <string>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <memory>
#include <unordered_map>
#include <cmath>
#include "json/json.h"

class System {

public:
    System();
    ~System();

    nlohmann::json GetJsonConfig();

    void SetUpFromJson(std::string func, char *json_file = nullptr);
    void SetUpFromJson(std::string json_text, std::string func);
    void SetUpFromJson(nlohmann::json j, std::string func);

    std::string GetJsonText();

    std::vector<double> GetBoxVector();
    std::vector<double> GetBoxMaxMin();
    size_t GetImcon();
    double GetDumpFrequency();
    double GetIntMax();
    double GetTemperature();
    double GetTau();
    double GetSmoothingFwhmCm1();
    bool HasSmoothingFwhmCm1();
    double GetSpectrumMaxWavenumberCm1();
    size_t GetSpectrumNumWavenumbers();
    bool GetSubtractInitialResponse();
    bool HasSubtractInitialResponse();
    double GetSpectrumTransformTimePs();
    bool HasSpectrumTransformTimePs();
    bool GetDMS();
    bool GetSubtractMeanDipole();
    bool GetSubtractMeanPolarizability();
    bool GetBoseEinsteinCorrection();
    bool GetWriteMolecularProperties();
    bool GetSubtractMeanVelocity();
    size_t GetSkipInitialAtoms();
    double GetTrajectoryDurationPs();
    double GetCutoff1();
    double GetCutoff2();
    double GetBoxVolume();
    std::string GetIjk();
    double GetZc1();
    double GetZc2();
    double GetCrossCorrelationCutoff();
    size_t GetCrossCorrelationSamplingStride();
    bool GetIncludeAutocorrelation();
    bool GetIncludeCrossCorrelation();
    bool GetCenterSlab();
    bool GetShiftSlabByHalfCell();
    bool GetTimeDependentSwitch();
    bool GetUsePsDipole();
    std::string GetDipoleMb();
    std::string GetPolarizabilityMb();
    std::string GetSurface();
    bool GetRamanResponse();

    
private:

    size_t imcon_;
    double dump_frequency_;
    double int_max_;
    double temperature_;
    double tau_;
    double smoothing_fwhm_cm1_;
    bool has_smoothing_fwhm_cm1_;
    double spectrum_max_wavenumber_cm1_;
    size_t spectrum_num_wavenumbers_;
    bool subtract_initial_response_;
    bool has_subtract_initial_response_;
    double spectrum_transform_time_ps_;
    bool has_spectrum_transform_time_ps_;
    bool dms_;
    bool subtract_mean_dipole_;
    bool subtract_mean_polarizability_;
    bool bose_einstein_correction_;
    bool write_molecular_properties_;
    bool subtract_mean_velocity_;
    size_t skip_initial_atoms_;
    double trajectory_duration_ps_;

    double cutoff1_;
    double cutoff2_;
    std::string ijk_;
    double zc1_;
    double zc2_;  
    double cross_correlation_cutoff_;
    size_t cross_correlation_sampling_stride_;
    bool include_autocorrelation_;
    bool include_cross_correlation_;
    bool center_slab_;
    bool shift_slab_by_half_cell_;
    bool time_dependent_switch_;
    bool use_ps_dipole_;
    std::string dipole_mb_;
    std::string polarizability_mb_;
    std::string surface_;

    bool initialized_;
    bool monomer_json_read_;
    bool mpi_initialized_;
    std::string raman_response_;

    std::vector<double> box_;
    std::vector<double> box_vector_;
    std::vector<double> box_maxmin_;

    nlohmann::json mbx_j_;

    // // ?
    // bool simcell_periodic_;
    // int tag_counter_;

    // ?
    std::vector<double> BoxVectorToMaxMin(const std::vector<double>& box, size_t imcon);
    std::vector<double> BoxMaxMinToVector(const std::vector<double>& box, size_t imcon);
};

#endif
