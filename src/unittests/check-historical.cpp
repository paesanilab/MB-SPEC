#include <set>

#include "test_support.h"

namespace {

using namespace test;

const fs::path historical = root / "historical";

struct PrintedNumber {
    double value;
    double unit;
};

std::vector<std::vector<PrintedNumber>> read_numbers(const fs::path& file, size_t columns)
{
    std::istringstream input(contents(file));
    std::vector<std::vector<PrintedNumber>> rows;
    for (std::string line; std::getline(input, line);) {
        const auto location = file.string() + ":" + std::to_string(rows.size() + 1);
        const auto fields = words(line);
        require(fields.size() == columns, location + ": wrong column count");
        auto& row = rows.emplace_back();
        for (const auto& field : fields) {
            std::istringstream number(field);
            double value;
            require(bool(number >> value) && number.eof() && std::isfinite(value),
                    location + ": invalid number '" + field + "'");

            double unit = 0.0;
            const auto exponent = field.find_first_of("eE");
            const auto decimal = field.find('.');
            if (decimal != std::string::npos) {
                const auto mantissa_end = exponent == std::string::npos ? field.size() : exponent;
                const auto digits = static_cast<int>(mantissa_end - decimal - 1);
                const auto power = exponent == std::string::npos ? 0 : std::stoi(field.substr(exponent + 1));
                unit = std::pow(10.0, double(power) - digits);
            }
            row.push_back({value, unit});
        }
    }
    require(!rows.empty(), file.string() + ": empty table");
    return rows;
}

void compare(const fs::path& actual_file, const fs::path& reference_file, size_t columns,
             double absolute_scale)
{
    const auto actual = read_numbers(actual_file, columns);
    const auto reference = read_numbers(reference_file, columns);
    require(actual.size() == reference.size(), reference_file.string() + ": different row counts");

    std::vector<double> scales(columns, 0.0);
    for (const auto& row : reference)
        for (size_t column = 0; column < columns; ++column)
            scales[column] = std::max(scales[column], std::abs(row[column].value));
    const double file_scale = *std::max_element(scales.begin(), scales.end());

    for (size_t row = 0; row < reference.size(); ++row) {
        for (size_t column = 0; column < columns; ++column) {
            const auto expected = reference[row][column];
            const double value = actual[row][column].value;
            const double scale = scales[column] > 0.0 ? scales[column] : file_scale;
            // Use the saved precision so coarser program output cannot widen the tolerance.
            const double tolerance = 1e-8 * std::max(std::abs(expected.value), std::abs(value))
                                   + absolute_scale * scale + expected.unit;
            if (std::abs(expected.value - value) > tolerance) {
                std::ostringstream message;
                message << reference_file << ':' << row + 1 << ':' << column + 1
                        << std::setprecision(17) << ": expected " << expected.value
                        << ", got " << value << " (tolerance " << tolerance << ") in " << actual_file;
                throw std::runtime_error(message.str());
            }
        }
    }
}

struct HistoricalSuite {
    Suite suite;
    fs::path references;
    Json config;
    std::set<fs::path>& checked_files;

    HistoricalSuite(const std::string& workflow, const fs::path& legacy_root,
                    std::set<fs::path>& checked)
        : suite("historical-" + workflow, "", legacy_root),
          references(historical / "references" / workflow),
          config(Json::parse(contents(historical / "inputs" / (workflow + ".json")))),
          checked_files(checked) {}

    void check_output(const fs::path& actual, const fs::path& name, size_t columns,
                      double absolute_scale = 1e-10)
    {
        const auto reference = references / name;
        compare(actual, reference, columns, absolute_scale);
        require(checked_files.insert(reference).second, "duplicate reference check: " + reference.string());
    }

    fs::path stage(const std::string& executable, const fs::path& input, const fs::path& name,
                   size_t columns, double absolute_scale = 1e-10)
    {
        const auto actual = suite.run(executable, input, config);
        check_output(actual, name, columns, absolute_scale);
        return actual;
    }
};

void check_ir(HistoricalSuite& test)
{
    const auto dipole = test.stage("infrared/dipole", historical / "inputs/water4.xyz", "dipole.dat", 4);
    // Mean subtraction and autocorrelation propagate rounding from the printed dipoles.
    const auto response = test.stage("infrared/response", dipole, "response.dat", 2, 1e-8);
    test.stage("infrared/spectrum", response, "spectrum.dat", 2);
}

void check_raman(HistoricalSuite& test)
{
    const auto polarizability = test.stage("raman/polarizability", historical / "inputs/water4.xyz",
                                           "polarizability.dat", 7);
    const auto depolarized = test.stage("raman/response", polarizability, "response.dat", 2);
    test.stage("raman/spectrum", depolarized, "spectrum.dat", 2);
    test.config["RAMAN"]["raman_response"] = "isotropic";
    const auto isotropic = test.stage("raman/response", polarizability, "response-iso.dat", 2);
    test.stage("raman/spectrum", isotropic, "spectrum-iso.dat", 2);
}

void check_sfg(HistoricalSuite& test)
{
    const auto directory = test.suite.work / "fresnel";
    fs::create_directories(directory);
    const auto optical = Json::parse(contents(historical / "inputs/fresnel.json"));
    for (const std::string component : {"zzz", "xxz", "xzx"}) {
        test.config["SFG"]["ijk"] = component;
        const auto response = test.stage("sfg/response", historical / "inputs/water4-sfg.lammpstrj",
                                         "response-" + component + ".dat", 4);
        const std::string name = optical.at(component + "_file");
        const auto spectrum = test.stage("sfg/spectrum", response, name, 4);
        fs::copy_file(spectrum, directory / name);
    }
    for (const std::string key : {"water_refractive_index_n_file", "water_refractive_index_k_file"}) {
        const fs::path table = optical.at(key).get<std::string>();
        fs::create_directories((directory / table).parent_path());
        fs::copy_file(root / "../sfg" / table, directory / table);
    }
    test.suite.run("sfg/fresnel", {}, optical, "", 1, "", directory);
    for (const std::string polarization : {"ssp", "pss", "sps", "ppp"}) {
        const fs::path name = optical.at(polarization + "_output_file").get<std::string>();
        test.check_output(directory / name, name, 2);
    }
}

void check_power(HistoricalSuite& test)
{
    const auto response = test.stage("power/response", historical / "inputs/water4-vel.lammpstrj",
                                     "response.dat", 2);
    test.stage("power/spectrum", response, "spectrum.dat", 2);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc > 2) {
        std::cerr << "Usage: " << argv[0] << " [legacy-src-directory]\n";
        return EXIT_FAILURE;
    }

    try {
        const auto legacy_root = argc == 2 ? fs::absolute(argv[1]) : fs::path{};
        const std::pair<const char*, void (*)(HistoricalSuite&)> cases[] = {
            {"ir", check_ir}, {"raman", check_raman}, {"sfg", check_sfg}, {"power", check_power},
        };
        int status = 0;
        std::set<fs::path> checked;
        for (const auto& workflow : cases) {
            HistoricalSuite test(workflow.first, legacy_root, checked);
            test.suite.check("outputs match the pinned master reference", [&] { workflow.second(test); });
            status |= test.suite.finish();
        }
        if (status == 0) {
            for (const auto& file : fs::recursive_directory_iterator(historical / "references")) {
                if (file.is_regular_file() && file.path().extension() == ".dat")
                    require(checked.count(file.path()), "unchecked reference: " + file.path().string());
            }
            std::cout << "Matched " << checked.size() << " historical reference files\n";
        }
        return status;
    } catch (const std::exception& error) {
        std::cerr << "Historical reference check: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
