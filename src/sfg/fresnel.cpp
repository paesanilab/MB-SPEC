#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cli.h"
#include "numeric_input.h"
#include "json/json.h"
#include "constants.h"

namespace {

using Complex = std::complex<double>;

struct ComponentSpectrum {
    std::vector<double> wavenumber;
    std::vector<double> real;
    std::vector<double> imaginary;
};

struct RefractiveIndexTable {
    std::vector<double> wavelength_um;
    std::vector<double> value;
};

template <typename T>
T ReadOption(nlohmann::json& config, const char* key) {
    const auto found = config.find(key);
    if (found == config.end())
        throw std::runtime_error(std::string("missing Fresnel option ") + key);
    try {
        const auto value = found->get<T>();
        config.erase(found);
        return value;
    } catch (const nlohmann::json::exception& error) {
        throw std::runtime_error(std::string("Fresnel option ") + key + ": " + error.what());
    }
}

ComponentSpectrum ReadComponentSpectrum(const std::string& filename) {
    std::ifstream input(filename.c_str());
    if (!input) {
        throw std::runtime_error("could not open '" + filename + "'");
    }

    ComponentSpectrum spectrum;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        std::istringstream row(line);
        double wavenumber;
        double imaginary;
        double real;
        const auto location = filename + ':' + std::to_string(line_number) + ": ";
        double magnitude_squared;
        if (!mbspec::read_finite(row, wavenumber) || !mbspec::read_finite(row, imaginary) ||
            !mbspec::read_finite(row, real) ||
            (!(row >> std::ws).eof() &&
             (!mbspec::read_finite(row, magnitude_squared) || !(row >> std::ws).eof())))
            throw std::runtime_error(location + "expected three or four finite spectrum columns");
        if (!spectrum.wavenumber.empty() && wavenumber <= spectrum.wavenumber.back())
            throw std::runtime_error(location + "wavenumbers must be strictly increasing");
        spectrum.wavenumber.push_back(wavenumber);
        spectrum.real.push_back(real);
        spectrum.imaginary.push_back(imaginary);
    }

    if (input.bad())
        throw std::runtime_error("could not read spectrum file '" + filename + "'");

    if (spectrum.wavenumber.empty()) {
        throw std::runtime_error("'" + filename + "' contains no spectrum rows");
    }
    return spectrum;
}

RefractiveIndexTable ReadRefractiveIndexTable(const std::string& filename) {
    std::ifstream input(filename.c_str());
    if (!input) {
        throw std::runtime_error("could not open '" + filename + "'");
    }

    RefractiveIndexTable table;
    double wavelength;
    double value;
    while (input >> wavelength >> value) {
        table.wavelength_um.push_back(wavelength);
        table.value.push_back(value);
    }

    if (table.wavelength_um.empty()) {
        throw std::runtime_error("'" + filename + "' contains no index rows");
    }
    if (!input.eof()) {
        throw std::runtime_error("could not parse an index row in '" + filename + "'");
    }
    for (size_t i = 1; i < table.wavelength_um.size(); ++i) {
        if (table.wavelength_um[i] <= table.wavelength_um[i - 1]) {
            throw std::runtime_error("wavelengths in '" + filename +
                                     "' must be strictly increasing");
        }
    }

    return table;
}

double Interpolate(const RefractiveIndexTable& table, double wavelength_um) {
    if (wavelength_um <= table.wavelength_um.front()) {
        return table.value.front();
    }
    if (wavelength_um >= table.wavelength_um.back()) {
        return table.value.back();
    }

    const std::vector<double>::const_iterator upper = std::upper_bound(
        table.wavelength_um.begin(), table.wavelength_um.end(), wavelength_um);
    const size_t right = static_cast<size_t>(upper - table.wavelength_um.begin());
    const size_t left = right - 1;
    const double fraction = (wavelength_um - table.wavelength_um[left]) /
        (table.wavelength_um[right] - table.wavelength_um[left]);
    return table.value[left] + fraction*(table.value[right] - table.value[left]);
}

Complex FresnelFactor(const std::string& polarization, const Complex& n1,
                      const Complex& n2, const Complex& ni, double theta) {
    const Complex gamma = std::asin(n1*std::sin(theta)/n2);
    if (polarization == "xx") {
        return 2.0*n1*std::cos(gamma) /
            (n1*std::cos(gamma) + n2*std::cos(theta));
    }
    if (polarization == "yy") {
        return 2.0*n1*std::cos(theta) /
            (n1*std::cos(theta) + n2*std::cos(gamma));
    }
    return 2.0*n2*std::cos(theta)*std::pow(n1/ni, 2.0) /
        (n1*std::cos(gamma) + n2*std::cos(theta));
}

Complex InterfacialIndex(const Complex& n) {
    return std::sqrt((std::pow(n, 4.0) + 5.0*std::pow(n, 2.0)) /
                     (4.0*std::pow(n, 2.0) + 2.0));
}

void WriteIntensity(std::ofstream& output, double wavenumber, const Complex& amplitude) {
    output << std::setprecision(17) << std::setw(8) << wavenumber << "     "
           << std::setw(16) << std::norm(amplitude) << '\n';
}

void CheckMatchingGrid(const ComponentSpectrum& reference,
                       const ComponentSpectrum& spectrum,
                       const std::string& filename) {
    if (spectrum.wavenumber.size() != reference.wavenumber.size()) {
        throw std::runtime_error("'" + filename +
                                 "' does not have the same number of rows as the zzz file");
    }
    for (size_t i = 0; i < reference.wavenumber.size(); ++i) {
        const double expected = reference.wavenumber[i];
        const double actual = spectrum.wavenumber[i];
        // mbspec-sfg-spectrum prints six significant digits.
        if (std::abs(actual - expected) > mbspec::printed_roundoff(expected, 6) +
                                         mbspec::printed_roundoff(actual, 6))
            throw std::runtime_error(filename + ':' + std::to_string(i + 1) +
                                     ": wavenumber grid differs from the zzz file");
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    if (const auto status = mbspec::check_arguments(argc, argv, 2, "<fresnel.json>"))
        return *status;

    try {
        std::ifstream config_input(argv[1]);
        if (!config_input) {
            throw std::runtime_error("could not open '" + std::string(argv[1]) + "'");
        }
        nlohmann::json config;
        config_input >> config;
        if (!config.is_object())
            throw std::runtime_error("Fresnel configuration must be a JSON object");

        const auto xxz_file = ReadOption<std::string>(config, "xxz_file");
        const auto zzz_file = ReadOption<std::string>(config, "zzz_file");
        const auto xzx_file = ReadOption<std::string>(config, "xzx_file");
        const auto ssp_output_file = ReadOption<std::string>(config, "ssp_output_file");
        const auto pss_output_file = ReadOption<std::string>(config, "pss_output_file");
        const auto sps_output_file = ReadOption<std::string>(config, "sps_output_file");
        const auto ppp_output_file = ReadOption<std::string>(config, "ppp_output_file");
        const auto water_n_file = ReadOption<std::string>(config, "water_refractive_index_n_file");
        const auto water_k_file = ReadOption<std::string>(config, "water_refractive_index_k_file");
        const double visible_angle_deg = ReadOption<double>(config, "visible_incidence_angle_deg");
        const double infrared_angle_deg = ReadOption<double>(config, "infrared_incidence_angle_deg");
        const double visible_wavelength_um = ReadOption<double>(config, "visible_wavelength_um");
        const double nonresonant_background = ReadOption<double>(config, "nonresonant_background");
        const double redshift_cm1 = ReadOption<double>(config, "redshift_cm1");
        if (!config.empty())
            throw std::runtime_error("unknown Fresnel option " + config.begin().key());

        const ComponentSpectrum xxz = ReadComponentSpectrum(xxz_file);
        const ComponentSpectrum zzz = ReadComponentSpectrum(zzz_file);
        const ComponentSpectrum xzx = ReadComponentSpectrum(xzx_file);
        CheckMatchingGrid(zzz, xxz, xxz_file);
        CheckMatchingGrid(zzz, xzx, xzx_file);

        const auto water_n = ReadRefractiveIndexTable(water_n_file);
        const auto water_k = ReadRefractiveIndexTable(water_k_file);

        std::ofstream ssp_output(ssp_output_file.c_str());
        std::ofstream pss_output(pss_output_file.c_str());
        std::ofstream sps_output(sps_output_file.c_str());
        std::ofstream ppp_output(ppp_output_file.c_str());
        if (!ssp_output || !pss_output || !sps_output || !ppp_output) {
            throw std::runtime_error("could not open one or more SFG output files");
        }

        const double pi = std::acos(-1.0);
        const double beta1 = visible_angle_deg*pi/180.0;
        const double beta2 = infrared_angle_deg*pi/180.0;
        const double speed_of_light = constants::c0; // m/s
        const double visible_frequency_ghz = speed_of_light/(visible_wavelength_um*1000.0);
        const double cm1_to_ghz = constants::cm_per_m*speed_of_light/1.0e9;
        const Complex air_index(1.0, 0.0);

        for (size_t i = 0; i < zzz.wavenumber.size(); ++i) {
            const double wavenumber_cm1 = zzz.wavenumber[i] - redshift_cm1;
            if (wavenumber_cm1 == 0.0) {
                continue;
            }
            const double infrared_wavelength_um = 1.0e4/wavenumber_cm1;
            const double infrared_frequency_ghz = wavenumber_cm1*cm1_to_ghz;
            const double sfg_frequency_ghz = visible_frequency_ghz + infrared_frequency_ghz;
            if (sfg_frequency_ghz == 0.0) {
                throw std::runtime_error("visible and infrared frequencies sum to zero");
            }
            const double sfg_wavelength_um = speed_of_light/(sfg_frequency_ghz*1000.0);
            const double beta = std::asin((visible_frequency_ghz*std::sin(beta1) + infrared_frequency_ghz*std::sin(beta2)) / sfg_frequency_ghz);

            const Complex n_visible(Interpolate(water_n, visible_wavelength_um),
                                    Interpolate(water_k, visible_wavelength_um));
            const Complex n_infrared(Interpolate(water_n, infrared_wavelength_um),
                                     Interpolate(water_k, infrared_wavelength_um));
            const Complex n_sfg(Interpolate(water_n, sfg_wavelength_um),
                                Interpolate(water_k, sfg_wavelength_um));
            const Complex ni_visible = InterfacialIndex(n_visible);
            const Complex ni_infrared = InterfacialIndex(n_infrared);
            const Complex ni_sfg = InterfacialIndex(n_sfg);

            const Complex lxx_sfg = FresnelFactor("xx", air_index, n_sfg, ni_sfg, beta);
            const Complex lxx_visible = FresnelFactor("xx", air_index, n_visible, ni_visible, beta1);
            const Complex lxx_infrared = FresnelFactor("xx", air_index, n_infrared, ni_infrared, beta2);
            const Complex lyy_sfg = FresnelFactor("yy", air_index, n_sfg, ni_sfg, beta);
            const Complex lyy_visible = FresnelFactor("yy", air_index, n_visible, ni_visible, beta1);
            const Complex lyy_infrared = FresnelFactor("yy", air_index, n_infrared, ni_infrared, beta2);
            const Complex lzz_sfg = FresnelFactor("zz", air_index, n_sfg, ni_sfg, beta);
            const Complex lzz_visible = FresnelFactor("zz", air_index, n_visible, ni_visible, beta1);
            const Complex lzz_infrared = FresnelFactor("zz", air_index, n_infrared, ni_infrared, beta2);

            const Complex x_xxz(xxz.real[i] + nonresonant_background, xxz.imaginary[i]);
            const Complex x_zzz(zzz.real[i] + nonresonant_background, zzz.imaginary[i]);
            const Complex x_xzx(xzx.real[i] + nonresonant_background, xzx.imaginary[i]);
            const Complex x_zxx = x_xzx;

            const Complex ssp = lzz_infrared*lyy_visible*lyy_sfg*
                std::sin(beta2)*x_xxz;
            const Complex pss = lyy_infrared*lyy_visible*lzz_sfg*
                std::sin(beta)*x_zxx;
            const Complex sps = lyy_infrared*lzz_visible*lyy_sfg*
                std::sin(beta1)*x_xzx;
            const Complex ppp =
                -lzz_infrared*lxx_visible*lxx_sfg*
                    std::cos(beta)*std::cos(beta1)*std::sin(beta2)*x_xxz
                -lxx_infrared*lzz_visible*lxx_sfg*
                    std::cos(beta)*std::sin(beta1)*std::cos(beta2)*x_xzx
                +lxx_infrared*lxx_visible*lzz_sfg*
                    std::sin(beta)*std::cos(beta1)*std::cos(beta2)*x_zxx
                +lzz_infrared*lzz_visible*lzz_sfg*
                    std::sin(beta)*std::sin(beta1)*std::sin(beta2)*x_zzz;

            WriteIntensity(ssp_output, wavenumber_cm1, ssp);
            WriteIntensity(pss_output, wavenumber_cm1, pss);
            WriteIntensity(sps_output, wavenumber_cm1, sps);
            WriteIntensity(ppp_output, wavenumber_cm1, ppp);
        }
    } catch (const std::exception& error) {
        std::cerr << "fresnel: " << error.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
