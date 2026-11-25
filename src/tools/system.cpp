#include "system.h"

#include <iomanip>
#include <initializer_list>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string_view>

#define EXIT_FAILURE 1

namespace {

template <typename T>
bool ReadJsonOption(nlohmann::json& options, const std::string& section,
                    const char* key, T& destination,
                    std::initializer_list<std::string_view> workflows = {}) {
    if (workflows.size() != 0 && std::find(workflows.begin(), workflows.end(), section) == workflows.end())
        return false;

    const auto found = options.find(key);
    if (found == options.end())
        return false;

    const auto& value = *found;
    const nlohmann::json default_value = T{};
    bool valid = value.type() == default_value.type();
    std::string expected = default_value.type_name();

    if (default_value.is_number_float()) {
        valid = value.is_number();
        expected = "a number";
    } else if (default_value.is_number_unsigned()) {
        valid = false;
        if (value.is_number_unsigned()) {
            valid = value.get<nlohmann::json::number_unsigned_t>() <=
                    std::numeric_limits<size_t>::max();
        } else if (value.is_number_integer()) {
            const auto integer = value.get<nlohmann::json::number_integer_t>();
            valid = integer >= 0 &&
                static_cast<nlohmann::json::number_unsigned_t>(integer) <=
                    std::numeric_limits<size_t>::max();
        }
        expected = "a nonnegative integer within the size_t range";
    } else if (default_value.is_array()) {
        valid = value.is_array() &&
            std::all_of(value.begin(), value.end(),
                        [](const nlohmann::json& element) { return element.is_number(); });
        expected = "an array of numbers";
    }

    if (!valid) {
        throw std::runtime_error(section + "." + key + " must be " + expected +
                                 " (received " + value.type_name() + ")");
    }
    destination = value.get<T>();
    options.erase(found);
    return true;
}

} // namespace

System::System() {
    initialized_ = false;
    monomer_json_read_ = false;
    mpi_initialized_ = false;
    // simcell_periodic_ = false;
    // tag_counter_ = 1;

    std::cerr << std::setprecision(20);

    cutoff1_ = 5.5;
    cutoff2_ = 7.5;

    imcon_ = 2;
    box_ = {0, 0, 0, 0, 0, 0};
    dump_frequency_ = 0.0;
    int_max_ = 0.0;
    temperature_ = 0.0;
    tau_ = 0.0;
    smoothing_fwhm_cm1_ = 0.0;
    has_smoothing_fwhm_cm1_ = false;
    spectrum_max_wavenumber_cm1_ = 0.0;
    spectrum_num_wavenumbers_ = 0;
    subtract_initial_response_ = false;
    has_subtract_initial_response_ = false;
    spectrum_transform_time_ps_ = 0.0;
    has_spectrum_transform_time_ps_ = false;
    dms_ = false;
    subtract_mean_dipole_ = false;
    subtract_mean_polarizability_ = false;
    bose_einstein_correction_ = false;
    write_molecular_properties_ = false;
    subtract_mean_velocity_ = false;
    skip_initial_atoms_ = 0;
    trajectory_duration_ps_ = 0.0;
    // Missing interface bounds must not act as physical default distances.
    zc1_ = std::numeric_limits<double>::quiet_NaN();
    zc2_ = std::numeric_limits<double>::quiet_NaN();
    cross_correlation_cutoff_ = 4.0;
    cross_correlation_sampling_stride_ = 250;
    include_autocorrelation_ = true;
    include_cross_correlation_ = true;
    center_slab_ = false;
    shift_slab_by_half_cell_ = false;
    time_dependent_switch_ = false;
    use_ps_dipole_ = false;
    dipole_mb_ = "";
    polarizability_mb_ = "";
    surface_ = "both";

    mbx_j_ = {};
}

System::~System() {}

nlohmann::json System::GetJsonConfig() {
    return mbx_j_;
}

void System::SetUpFromJson(nlohmann::json j, std::string func) {

    if (!j.is_object())
        throw std::runtime_error("JSON configuration must be an object");
    if (!j.contains(func))
        throw std::runtime_error("missing JSON workflow section " + func);
    auto& options = j.at(func);
    if (!options.is_object())
        throw std::runtime_error(func + " must be a JSON object");

    ReadJsonOption(options, func, "imcon", imcon_);

    mbx_j_["MB-SPEC"]["imcon"] = imcon_;

    ReadJsonOption(options, func, "dump_frequency", dump_frequency_);

    mbx_j_["MB-SPEC"]["dump_frequency"] = dump_frequency_;

    ReadJsonOption(options, func, "DMS", dms_, {"INFRARED"});

    mbx_j_["MB-SPEC"]["DMS"] = dms_;

    ReadJsonOption(options, func, "subtract_mean_dipole", subtract_mean_dipole_, {"INFRARED"});

    mbx_j_["MB-SPEC"]["subtract_mean_dipole"] = subtract_mean_dipole_;

    ReadJsonOption(options, func, "subtract_mean_polarizability",
                   subtract_mean_polarizability_, {"RAMAN"});

    mbx_j_["MB-SPEC"]["subtract_mean_polarizability"] =
        subtract_mean_polarizability_;

    ReadJsonOption(options, func, "write_molecular_properties", write_molecular_properties_, {"INFRARED", "RAMAN"});

    mbx_j_["MB-SPEC"]["write_molecular_properties"] =
        write_molecular_properties_;

    ReadJsonOption(options, func, "subtract_mean_velocity", subtract_mean_velocity_, {"POWER"});

    mbx_j_["MB-SPEC"]["subtract_mean_velocity"] = subtract_mean_velocity_;

    ReadJsonOption(options, func, "skip_initial_atoms", skip_initial_atoms_);

    mbx_j_["MB-SPEC"]["skip_initial_atoms"] = skip_initial_atoms_;

    ReadJsonOption(options, func, "trajectory_duration_ps", trajectory_duration_ps_);

    mbx_j_["MB-SPEC"]["trajectory_duration_ps"] = trajectory_duration_ps_;

    ReadJsonOption(options, func, "int_max", int_max_);

    mbx_j_["MB-SPEC"]["int_max"] = int_max_;

    ReadJsonOption(options, func, "temperature", temperature_, {"INFRARED", "RAMAN", "SFG"});

    mbx_j_["MB-SPEC"]["temperature"] = temperature_;

    ReadJsonOption(options, func, "tau", tau_);

    mbx_j_["MB-SPEC"]["tau"] = tau_;

    if (ReadJsonOption(options, func, "smoothing_fwhm_cm1", smoothing_fwhm_cm1_)) {
        has_smoothing_fwhm_cm1_ = true;
    }

    if (has_smoothing_fwhm_cm1_) {
        mbx_j_["MB-SPEC"]["smoothing_fwhm_cm1"] = smoothing_fwhm_cm1_;
    }

    ReadJsonOption(options, func, "spectrum_max_wavenumber_cm1",
                   spectrum_max_wavenumber_cm1_, {"SFG"});

    mbx_j_["MB-SPEC"]["spectrum_max_wavenumber_cm1"] =
        spectrum_max_wavenumber_cm1_;

    ReadJsonOption(options, func, "spectrum_num_wavenumbers", spectrum_num_wavenumbers_, {"SFG"});

    mbx_j_["MB-SPEC"]["spectrum_num_wavenumbers"] =
        spectrum_num_wavenumbers_;

    if (ReadJsonOption(options, func, "subtract_initial_response",
                       subtract_initial_response_, {"SFG"})) {
        has_subtract_initial_response_ = true;
    }

    if (has_subtract_initial_response_) {
        mbx_j_["MB-SPEC"]["subtract_initial_response"] =
            subtract_initial_response_;
    }

    if (ReadJsonOption(options, func, "spectrum_transform_time_ps",
                       spectrum_transform_time_ps_)) {
        has_spectrum_transform_time_ps_ = true;
    }

    if (has_spectrum_transform_time_ps_) {
        mbx_j_["MB-SPEC"]["spectrum_transform_time_ps"] =
            spectrum_transform_time_ps_;
    }

    ReadJsonOption(options, func, "cutoff_1", cutoff1_, {"INFRARED", "RAMAN", "SFG"});

    mbx_j_["MB-SPEC"]["cutoff_1"] = cutoff1_;

    ReadJsonOption(options, func, "cutoff_2", cutoff2_, {"INFRARED", "RAMAN", "SFG"});

    mbx_j_["MB-SPEC"]["cutoff_2"] = cutoff2_;

    ReadJsonOption(options, func, "ijk", ijk_, {"SFG"});

    mbx_j_["MB-SPEC"]["ijk"] = ijk_;

    ReadJsonOption(options, func, "raman_response", raman_response_, {"RAMAN"});

    mbx_j_["MB-SPEC"]["raman_response"] = raman_response_;

    ReadJsonOption(options, func, "bose_einstein_correction", bose_einstein_correction_, {"RAMAN"});

    mbx_j_["MB-SPEC"]["bose_einstein_correction"] = bose_einstein_correction_;

    ReadJsonOption(options, func, "zc1", zc1_, {"SFG"});

    mbx_j_["MB-SPEC"]["zc1"] = zc1_;

    ReadJsonOption(options, func, "zc2", zc2_, {"SFG"});

    mbx_j_["MB-SPEC"]["zc2"] = zc2_;

    ReadJsonOption(options, func, "cross_correlation_cutoff", cross_correlation_cutoff_, {"SFG"});

    mbx_j_["MB-SPEC"]["cross_correlation_cutoff"] =
        cross_correlation_cutoff_;

    ReadJsonOption(options, func, "cross_correlation_sampling_stride",
                   cross_correlation_sampling_stride_, {"SFG"});

    mbx_j_["MB-SPEC"]["cross_correlation_sampling_stride"] =
        cross_correlation_sampling_stride_;

    ReadJsonOption(options, func, "include_autocorrelation", include_autocorrelation_, {"SFG"});

    mbx_j_["MB-SPEC"]["include_autocorrelation"] =
        include_autocorrelation_;

    ReadJsonOption(options, func, "include_cross_correlation", include_cross_correlation_, {"SFG"});

    mbx_j_["MB-SPEC"]["include_cross_correlation"] =
        include_cross_correlation_;

    ReadJsonOption(options, func, "center_slab", center_slab_, {"SFG"});

    mbx_j_["MB-SPEC"]["center_slab"] = center_slab_;

    ReadJsonOption(options, func, "shift_slab_by_half_cell", shift_slab_by_half_cell_, {"SFG"});

    mbx_j_["MB-SPEC"]["shift_slab_by_half_cell"] =
        shift_slab_by_half_cell_;

    ReadJsonOption(options, func, "time_dependent_switch", time_dependent_switch_, {"SFG"});

    mbx_j_["MB-SPEC"]["time_dependent_switch"] =
        time_dependent_switch_;

    ReadJsonOption(options, func, "use_ps_dipole", use_ps_dipole_, {"SFG"});

    mbx_j_["MB-SPEC"]["use_ps_dipole"] = use_ps_dipole_;

    ReadJsonOption(options, func, "dipole_mb", dipole_mb_, {"INFRARED", "SFG"});
    mbx_j_["MB-SPEC"]["dipole_mb"] = dipole_mb_;

    ReadJsonOption(options, func, "polarizability_mb", polarizability_mb_, {"RAMAN", "SFG"});
    mbx_j_["MB-SPEC"]["polarizability_mb"] = polarizability_mb_;

    ReadJsonOption(options, func, "surface", surface_, {"SFG"});
    mbx_j_["MB-SPEC"]["surface"] = surface_;

    ReadJsonOption(options, func, "box", box_);
    if (!options.empty())
        throw std::runtime_error("unknown option " + func + "." + options.begin().key());

    // if (box_.size() == 9) {
    //     box_maxmin_ = BoxVectorToMaxMin(box_, imcon_);
    // } 
    // else if (box_.size() == 6) {
    //     box_maxmin_ = box_;
    //     box_ = BoxMaxMinToVector(box_maxmin_, imcon_);
    // }


    mbx_j_["MB-SPEC"]["box"] = box_;
}

void System::SetUpFromJson(std::string json_text, std::string func) {
    nlohmann::json j = nlohmann::json::parse(json_text);
    SetUpFromJson(j, func);
}

void System::SetUpFromJson(std::string func, char *json_file) {
    try {
        nlohmann::json j = {{func, nlohmann::json::object()}};
        if (json_file != nullptr) {
            std::ifstream ifjson(json_file);
            if (!ifjson)
                throw std::runtime_error("could not open JSON file '" +
                                         std::string(json_file) + "'");
            j = nlohmann::json::parse(ifjson);
        }
        SetUpFromJson(j, func);
    } catch (const std::exception& error) {
        std::cerr << "Configuration error";
        if (json_file != nullptr)
            std::cerr << " in '" << json_file << "'";
        std::cerr << ": " << error.what() << std::endl;
        std::exit(EXIT_FAILURE);
    } 
}

std::string System::GetJsonText() {
    std::ostringstream ss;
    ss << mbx_j_;
    return ss.str();
}

std::vector<double> System::GetBoxVector() {
    if (box_.size() == 9) {
        return box_;
    } else if (box_.size() == 6) {
        return BoxMaxMinToVector(box_, imcon_);
    } else {
        throw std::runtime_error("box must be in {xx, xy, xz, yx, yy, yz, zx, zy, zz} or {xmin, xmax, ymin, ymax, zmin, zmax} format");
    }
}

std::vector<double> System::GetBoxMaxMin() {
    if (box_.size() == 6) {
        return box_;
    } else if (box_.size() == 9) {
        return BoxVectorToMaxMin(box_, imcon_);
    } else {
        throw std::runtime_error("box must be in {xx, xy, xz, yx, yy, yz, zx, zy, zz} or {xmin, xmax, ymin, ymax, zmin, zmax} format");
    }
}

size_t System::GetImcon() {
    return imcon_;
}

double System::GetDumpFrequency () {
    return dump_frequency_;
}

double System::GetIntMax() {
    return int_max_;
}

double System::GetTemperature() {
    return temperature_;
}

double System::GetTau() {
    return tau_;
}

double System::GetSmoothingFwhmCm1() {
    return smoothing_fwhm_cm1_;
}

bool System::HasSmoothingFwhmCm1() {
    return has_smoothing_fwhm_cm1_;
}

double System::GetSpectrumMaxWavenumberCm1() {
    return spectrum_max_wavenumber_cm1_;
}

size_t System::GetSpectrumNumWavenumbers() {
    return spectrum_num_wavenumbers_;
}

bool System::GetSubtractInitialResponse() {
    return subtract_initial_response_;
}

bool System::HasSubtractInitialResponse() {
    return has_subtract_initial_response_;
}

double System::GetSpectrumTransformTimePs() {
    return spectrum_transform_time_ps_;
}

bool System::HasSpectrumTransformTimePs() {
    return has_spectrum_transform_time_ps_;
}

bool System::GetDMS() {
    return dms_;
}

bool System::GetSubtractMeanDipole() {
    return subtract_mean_dipole_;
}

bool System::GetSubtractMeanPolarizability() {
    return subtract_mean_polarizability_;
}

bool System::GetBoseEinsteinCorrection() {
    return bose_einstein_correction_;
}

bool System::GetWriteMolecularProperties() {
    return write_molecular_properties_;
}

bool System::GetSubtractMeanVelocity() {
    return subtract_mean_velocity_;
}

size_t System::GetSkipInitialAtoms() {
    return skip_initial_atoms_;
}

double System::GetTrajectoryDurationPs() {
    return trajectory_duration_ps_;
}

double System::GetCutoff1() {
    return cutoff1_;
}

double System::GetCutoff2() {
    return cutoff2_;
}


// TODO: change box so xx xy xz yx yy yz zx zy zz then can calc for monoclinic
double System::GetBoxVolume() {
    if (box_.size() == 9) {
        if (imcon_ < 3) {
            return box_[0] * box_[4] * box_[8];
        }
        else {
            double cross_x = box_[4]*box_[8] - box_[5]*box_[7];
            double cross_y = -box_[1]*box_[8] + box_[7]*box_[2];
            double cross_z = box_[1]*box_[5] - box_[4]*box_[2];

            return box_[0]*cross_x + box_[3]*cross_y + box_[6]*cross_z;
        }
    }
    if (box_.size() == 6) {
        if (imcon_ < 3) {
            box_vector_ = BoxMaxMinToVector(box_, imcon_);
            double vmm = (box_[1] - box_[0]) * (box_[3] - box_[2]) * (box_[5] - box_[4]);
            double vv = box_vector_[0] * box_vector_[4] * box_vector_[8];
            if (vmm == vv) {
                return vmm;
            } else {
                throw std::runtime_error("Cannot calculate box volume using Min/Max coordinates and imcon = " + std::to_string(imcon_));
            }
        } else {
            throw std::runtime_error("Cannot calculate box volume using Min/Max coordinates and imcon = " + std::to_string(imcon_));
        }
    }

    throw std::runtime_error("Cannot calculate box volume: expected a 6- or 9-element box");
}

std::string System::GetIjk() {
        return ijk_;
}

double System::GetZc1() {
    return zc1_;
}

double System::GetZc2() {
    return zc2_;
}

double System::GetCrossCorrelationCutoff() {
    return cross_correlation_cutoff_;
}

size_t System::GetCrossCorrelationSamplingStride() {
    return cross_correlation_sampling_stride_;
}

bool System::GetIncludeAutocorrelation() {
    return include_autocorrelation_;
}

bool System::GetIncludeCrossCorrelation() {
    return include_cross_correlation_;
}

bool System::GetCenterSlab() {
    return center_slab_;
}

bool System::GetShiftSlabByHalfCell() {
    return shift_slab_by_half_cell_;
}

bool System::GetTimeDependentSwitch() {
    return time_dependent_switch_;
}

bool System::GetUsePsDipole() {
    return use_ps_dipole_;
}

std::string System::GetDipoleMb() {
    return dipole_mb_;
}

std::string System::GetPolarizabilityMb() {
    return polarizability_mb_;
}

std::string System::GetSurface() {
    return surface_;
}

// true if depol
// TODO: .tolower() catch all 'depol' 'dep' 'd' t/f 'iso' etc
bool System::GetRamanResponse() {
    if (raman_response_.compare("depolarized") == 0) {
        return true;
    }
    return false;
}

std::vector<double> System::BoxVectorToMaxMin(const std::vector<double>& box, size_t imcon) {
    if (imcon < 3) {
        std::vector<double> box_out = {0, box[0], 0, box[4], 0, box[8]};
        return box_out;
    } else {
        throw std::runtime_error("Cannot calculate box using vector coordinates and imcon = " + std::to_string(imcon));
    }
}

std::vector<double> System::BoxMaxMinToVector(const std::vector<double>& box, size_t imcon) {
    if (imcon < 3) {
        std::vector<double> box_out = {box[1] - box[0], 0, 0, 0, box[3] - box[2], 0, 0, 0, box[5] - box[4]};
        return box_out;
    } else {
        throw std::runtime_error("Cannot calculate box using Min/Max coordinates and imcon = " + std::to_string(imcon));
    }
    
}
