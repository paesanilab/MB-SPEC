#ifndef MB_SPEC_TEST_SUPPORT_H
#define MB_SPEC_TEST_SUPPORT_H

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
#include <cerrno>
#include <csignal>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>
#include "json/json.h"

namespace test {
namespace fs = std::filesystem;
using Json = nlohmann::json;
using Table = std::vector<std::vector<double>>;
inline const fs::path root = UNITTEST_DIR;
inline const fs::path build_root = UNITTEST_BUILDDIR;
inline const fs::path program_root = UNITTEST_BINDIR;
inline constexpr double dt = 0.004;

inline void require(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}

inline void near(double actual, double expected, double relative = 2e-5, double absolute = 1e-12) {
    if (!std::isfinite(actual) || !std::isfinite(expected) ||
        std::abs(actual - expected) > absolute + relative * std::abs(expected)) {
        std::ostringstream error;
        error << std::setprecision(15) << "expected " << expected << ", got " << actual;
        throw std::runtime_error(error.str());
    }
}

inline std::string contents(const fs::path& path) {
    std::ifstream in(path);
    require(bool(in), "cannot open " + path.string());
    return std::string(std::istreambuf_iterator<char>(in), {});
}

inline void write(const fs::path& path, const std::string& text) {
    std::ofstream out(path);
    require(bool(out), "cannot write " + path.string());
    out << text;
    require(bool(out), "write failed: " + path.string());
}

inline Table table(const fs::path& path, size_t columns) {
    std::istringstream in(contents(path));
    Table result;
    for (std::string line; std::getline(in, line);) {
        if (line.empty() || line[0] == '#')
            continue;
        std::istringstream row(line);
        std::vector<double> values;
        double value;
        while (row >> value) {
            require(std::isfinite(value), "non-finite output in " + path.string());
            values.push_back(value);
        }
        require(row.eof() && values.size() == columns, "wrong output columns: " + path.string());
        result.push_back(values);
    }
    require(!result.empty(), "empty output: " + path.string());
    return result;
}

inline Table table(const Json& saved, size_t columns) {
    require(saved.at("columns").size() == columns, "wrong saved-reference column count");
    Table result;
    for (const auto& row : saved.at("rows")) {
        require(row.size() == columns, "wrong saved-reference row width");
        std::vector<double> values;
        for (const auto& value : row) {
            require(value.is_number(), "non-numeric saved-reference value");
            double number = value.get<double>();
            require(std::isfinite(number), "non-finite saved-reference value");
            values.push_back(number);
        }
        result.push_back(std::move(values));
    }
    require(!result.empty(), "empty saved reference");
    return result;
}

inline void same(const Table& actual, const Table& expected, double relative = 2e-5,
                 double absolute = 1e-10) {
    require(actual.size() == expected.size(), "different row counts");
    for (size_t i = 0; i < actual.size(); ++i) {
        require(actual[i].size() == expected[i].size(), "different column counts");
        for (size_t j = 0; j < actual[i].size(); ++j)
            near(actual[i][j], expected[i][j], relative, absolute);
    }
}

inline std::string serialize(const Table& rows) {
    std::ostringstream out;
    out << std::setprecision(17);
    for (const auto& row : rows) {
        for (double value : row)
            out << value << ' ';
        out << '\n';
    }
    return out.str();
}

struct Frame {
    std::vector<std::string> header, atoms;
    std::vector<std::string> columns;
    size_t column(const std::string& name) const {
        auto found = std::find(columns.begin(), columns.end(), name);
        require(found != columns.end(), "missing trajectory column " + name);
        return found - columns.begin();
    }
};

inline std::vector<std::string> words(const std::string& line) {
    std::istringstream in(line);
    std::vector<std::string> result;
    for (std::string word; in >> word;)
        result.push_back(word);
    return result;
}

inline std::vector<Frame> trajectory(const fs::path& path) {
    std::ifstream in(path);
    require(bool(in), "missing input " + path.string());
    const bool xyz = path.extension() == ".xyz";
    std::vector<Frame> frames;
    for (std::string line; std::getline(in, line);) {
        Frame frame;
        frame.header.push_back(line);
        for (int i = 1; i < (xyz ? 2 : 9); ++i) {
            require(bool(std::getline(in, line)), "truncated trajectory header");
            frame.header.push_back(line);
        }
        const size_t count = std::stoul(frame.header[xyz ? 0 : 3]);
        if (!xyz) {
            require(frame.header[0] == "ITEM: TIMESTEP", "invalid LAMMPS frame");
            auto header = words(frame.header[8]);
            frame.columns.assign(header.begin() + 2, header.end());
        }
        for (size_t i = 0; i < count; ++i) {
            require(bool(std::getline(in, line)), "truncated trajectory atoms");
            require(words(line).size() == (xyz ? 4 : frame.columns.size()), "wrong atom columns");
            frame.atoms.push_back(line);
        }
        frames.push_back(frame);
    }
    require(frames.size() == 10, path.filename().string() + " must have exactly ten frames");
    return frames;
}

inline std::string encode(std::vector<Frame> frames) {
    std::ostringstream out;
    for (auto& frame : frames) {
        const bool xyz = frame.columns.empty();
        frame.header[xyz ? 0 : 3] = std::to_string(frame.atoms.size());
        if (!xyz) {
            frame.header[8] = "ITEM: ATOMS";
            for (const auto& column : frame.columns)
                frame.header[8] += " " + column;
        }
        for (const auto& line : frame.header)
            out << line << '\n';
        for (const auto& line : frame.atoms)
            out << line << '\n';
    }
    return out.str();
}

inline std::vector<Frame> subset(std::vector<Frame> frames, size_t skip, size_t count) {
    for (auto& frame : frames) {
        require(skip + count <= frame.atoms.size(), "invalid fixture subset");
        frame.atoms = std::vector<std::string>(frame.atoms.begin() + skip,
                                               frame.atoms.begin() + skip + count);
    }
    return frames;
}

inline std::vector<Frame> nearby_waters(std::vector<Frame> frames) {
    frames = subset(std::move(frames), 0, 6);
    // Tilt the water plane so the internal PS dipole has a nonzero z component.
    const double positions[18] = {5, 5, 12, 5.9572, 5, 12, 4.760013, 5, 12.927297,
                                  8, 5, 12, 8.9572, 5, 12, 7.760013, 5, 12.927297};
    for (size_t i = 0; i < frames.size(); ++i) {
        auto& frame = frames[i];
        for (size_t atom = 0; atom < 6; ++atom) {
            auto row = words(frame.atoms[atom]);
            for (size_t c = 0; c < 3; ++c) {
                std::ostringstream value;
                value << std::setprecision(17)
                      << positions[3 * atom + c] + ((atom % 3 == 1 && c == 0) ? .001 * i : 0);
                row[frame.column(std::array<std::string, 3>{"x", "y", "z"}[c])] = value.str();
            }
            frame.atoms[atom].clear();
            for (const auto& value : row)
                frame.atoms[atom] += value + " ";
        }
    }
    return frames;
}

inline Table vectors(const std::vector<Frame>& frames, const std::array<std::string, 3>& fields,
                     size_t skip = 0) {
    Table result;
    for (size_t i = 0; i < frames.size(); ++i) {
        std::vector<double> row = {i * dt, 0, 0, 0};
        const auto& frame = frames[i];
        for (size_t atom = skip; atom < frame.atoms.size(); ++atom) {
            auto values = words(frame.atoms[atom]);
            for (size_t c = 0; c < 3; ++c)
                row[c + 1] +=
                    std::stod(values[frame.columns.empty() ? c + 1 : frame.column(fields[c])]);
        }
        result.push_back(row);
    }
    return result;
}

// Deliberately calculate expectations from parsed input, not from saved program output.
inline Table correlation(Table properties, bool subtract_mean, double unit = 1.0) {
    const size_t dimensions = properties[0].size() - 1;
    for (size_t c = 1; c <= dimensions; ++c) {
        double mean = 0;
        if (subtract_mean) {
            for (const auto& row : properties)
                mean += row[c] / properties.size();
        }
        for (auto& row : properties)
            row[c] = (row[c] - mean) * unit;
    }
    Table result;
    for (size_t lag = 0; lag < 4; ++lag) {
        double value = 0;
        for (size_t i = 0; i + lag < properties.size(); ++i)
            for (size_t c = 1; c <= dimensions; ++c)
                value += properties[i][c] * properties[i + lag][c];
        result.push_back({lag * dt, value / (properties.size() - lag)});
    }
    return result;
}

struct Suite {
    std::string name, section;
    fs::path work;
    fs::path legacy_program_root;
    size_t passed = 0, failed = 0, sequence = 0;
    explicit Suite(std::string workflow, std::string config_section, fs::path legacy_root = {})
        : name(std::move(workflow)), section(std::move(config_section)),
          legacy_program_root(std::move(legacy_root)) {
        const auto run_id = std::chrono::system_clock::now().time_since_epoch().count();
        work = build_root / ".build" / "runs" /
               (name + "-" + std::to_string(getpid()) + "-" + std::to_string(run_id));
        require(fs::create_directories(work), "run directory already exists: " + work.string());
    }
    fs::path input(const std::string& filename) const { return root / "inputs" / filename; }
    Json config() const {
        auto config = Json::parse(contents(input(name + ".json")));
        auto& options = config[section];
        options["int_max"] = 4 * dt;
        options["spectrum_transform_time_ps"] = 4 * dt;
        options["trajectory_duration_ps"] = 0.0;
        options["dump_frequency"] = dt;
        options["skip_initial_atoms"] = 0;
        options["tau"] = 0.0;
        options["smoothing_fwhm_cm1"] = 0.0;
        if (section == "SFG")
            options["spectrum_num_wavenumbers"] = 33;
        return config;
    }
    template <typename F> void check(const std::string& label, F action) {
        try {
            action();
            ++passed;
            std::cout << "  PASS " << label << '\n';
        } catch (const std::exception& error) {
            ++failed;
            std::cerr << "  FAIL " << label << ": " << error.what() << '\n';
        }
    }
    fs::path fixture(const std::string& filename, const std::string& text) const {
        fs::path path = work / filename;
        write(path, text);
        return path;
    }
    fs::path run(const std::string& executable, const fs::path& input_path, const Json& config,
                 const std::string& expected_error = "", int threads = 1,
                 const std::string& raw_config = "", const fs::path& working_directory = {}) {
        const auto directory = working_directory.empty() ? work / std::to_string(++sequence)
                                                         : working_directory;
        fs::create_directories(directory);
        auto config_path = directory / "config.json";
        write(config_path, raw_config.empty() ? config.dump(2) : raw_config);
        auto output = directory / "stdout.dat", log = directory / "stderr.log";
        const fs::path stage = executable;
        const auto workflow = stage.parent_path().string();
        const auto program = legacy_program_root.empty()
            ? program_root / workflow /
                ("mbspec-" + (workflow == "infrared" ? "ir" : workflow) + "-" + stage.filename().string())
            : legacy_program_root / stage;
        require(fs::exists(program), "missing executable " + program.string());
        pid_t child = fork();
        require(child >= 0, "fork failed");
        if (child == 0) {
            if (chdir(directory.c_str()) != 0)
                _exit(125);
            int out = open(output.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
            int err = open(log.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0600);
            if (out < 0 || err < 0 || dup2(out, STDOUT_FILENO) < 0 || dup2(err, STDERR_FILENO) < 0)
                _exit(125);
            close(out);
            close(err);
            setenv("OMP_NUM_THREADS", std::to_string(threads).c_str(), 1);
            setenv("OMP_DYNAMIC", "FALSE", 1);
            setenv("MKL_NUM_THREADS", "1", 1);
            if (executable == "sfg/fresnel")
                execl(program.c_str(), program.c_str(), config_path.c_str(),
                      static_cast<char*>(nullptr));
            else
                execl(program.c_str(), program.c_str(), input_path.c_str(), config_path.c_str(),
                      static_cast<char*>(nullptr));
            _exit(125);
        }
        int status = 0;
        const int timeout_seconds = std::getenv("MB_SPEC_EXTENDED") ? 180 : 20;
        auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(timeout_seconds);
        while (true) {
            auto result = waitpid(child, &status, WNOHANG);
            if (result == child)
                break;
            if (result < 0 && errno != EINTR)
                throw std::runtime_error("waitpid failed");
            if (std::chrono::steady_clock::now() > deadline) {
                kill(child, SIGKILL);
                waitpid(child, &status, 0);
                throw std::runtime_error(std::to_string(timeout_seconds) + "-second timeout; " +
                                         log.string());
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        const std::string diagnostic = contents(log);
        require(WIFEXITED(status), "program crashed; " + log.string());
        if (expected_error.empty())
            require(WEXITSTATUS(status) == 0,
                    "program failed; " + log.string() + "\n" + diagnostic);
        else {
            require(WEXITSTATUS(status) != 0 && WEXITSTATUS(status) != 125,
                    "expected rejection; " + log.string());
            require(diagnostic.find(expected_error) != std::string::npos,
                    "wrong diagnostic; " + log.string() + "\n" + diagnostic);
        }
        return output;
    }
    int finish() const {
        std::cout << name << ": " << passed << " passed, " << failed << " failed\n"
                  << "  Logs and generated fixtures: " << work << '\n';
        return failed ? 1 : 0;
    }
};

inline void check_times(const Table& rows, size_t expected) {
    require(rows.size() == expected, "unexpected number of time samples");
    for (size_t i = 0; i < rows.size(); ++i)
        near(rows[i][0], i * dt);
}

inline Table synthetic_properties(size_t dimensions) {
    Table rows;
    for (size_t i = 0; i < 10; ++i) {
        std::vector<double> row = {i * dt};
        for (size_t c = 0; c < dimensions; ++c)
            row.push_back(std::sin(0.7 * i + c) + 0.15 * i + c);
        rows.push_back(row);
    }
    return rows;
}

// Independent Gaussian reference, including the documented exclusion of DC.
inline Table smooth_reference(Table rows, double fwhm, size_t column) {
    const double sigma = fwhm / (2 * std::sqrt(2 * std::log(2.0)));
    const double spacing = rows[1][0] - rows[0][0];
    const int radius = static_cast<int>(std::ceil(3 * sigma / spacing));
    auto result = rows;
    result[0][column] = 0;
    for (size_t i = 1; i < rows.size(); ++i) {
        double sum = 0, weights = 0;
        for (size_t j = 1; j < rows.size(); ++j) {
            if (std::abs(static_cast<int>(i) - static_cast<int>(j)) > radius)
                continue;
            const double d = (rows[i][0] - rows[j][0]) / sigma;
            const double weight = std::exp(-0.5 * d * d);
            sum += weight * rows[j][column];
            weights += weight;
        }
        result[i][column] = sum / weights;
    }
    return result;
}

inline std::string first_stage(const Suite& suite) {
    if (suite.name == "ir")
        return "infrared/dipole";
    if (suite.name == "raman")
        return "raman/polarizability";
    return suite.name + "/response";
}

inline std::string join_words(const std::vector<std::string>& values) {
    std::string line;
    for (const auto& value : values)
        line += value + " ";
    return line;
}

inline std::string decimal(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

// Metamorphic checks: change input representation while preserving water geometry.
inline void geometry_checks(Suite& suite, const std::vector<Frame>& source) {
    auto frames = nearby_waters(source);
    frames.resize(2);
    auto config = suite.config();
    config[suite.section]["int_max"] = 2 * dt;
    if (suite.name == "ir") {
        config[suite.section]["DMS"] = true;
        config[suite.section]["dipole_mb"] = "1b+2b";
    } else if (suite.name == "raman") {
        config[suite.section]["polarizability_mb"] = "1b+2b";
    } else {
        config[suite.section]["dipole_mb"] = "1b+2b";
        config[suite.section]["polarizability_mb"] = "1b+2b";
        config[suite.section]["center_slab"] = false;
        config[suite.section]["zc1"] = 0.;
        config[suite.section]["zc2"] = .1;
        config[suite.section]["include_cross_correlation"] = true;
        config[suite.section]["cross_correlation_sampling_stride"] = 1;
    }
    const size_t columns = suite.name == "raman" ? 7 : 4;
    auto original = suite.fixture("geometry-original.lammpstrj", encode(frames));
    auto evaluate = [&](const fs::path& input, const Json& options, int threads = 1) {
        return table(suite.run(first_stage(suite), input, options, "", threads), columns);
    };
    Table baseline;
    auto get_baseline = [&]() -> const Table& {
        if (baseline.empty()) {
            auto rows = evaluate(original, config);
            require(rows.size() == 2, "expected two frames/lags");
            require(std::abs(rows[0][1]) > 1e-10, "geometry baseline must be nonzero");
            baseline = std::move(rows);
        }
        return baseline;
    };
    auto compare = [&](std::vector<Frame> changed, const std::string& filename) {
        get_baseline();
        same(evaluate(suite.fixture(filename, encode(changed)), config), baseline, 3e-5, 1e-8);
    };
    suite.check("periodic imaging of individually shifted atoms and whole waters", [&] {
        auto shifted = frames;
        for (auto& frame : shifted)
            for (size_t atom = 0; atom < frame.atoms.size(); ++atom) {
                auto values = words(frame.atoms[atom]);
                for (size_t c = 0; c < 3; ++c) {
                    auto bounds = words(frame.header[5 + c]);
                    double length = std::stod(bounds[1]) - std::stod(bounds[0]);
                    // Oxygen z positions fix SFG's interface; move hydrogens in z only.
                    int images = c == 2 && atom % 3 == 0 ? 0 : (atom % 2 ? 2 : -1);
                    size_t index = frame.column(std::array<std::string, 3>{"x", "y", "z"}[c]);
                    values[index] = decimal(std::stod(values[index]) + images * length);
                }
                frame.atoms[atom] = join_words(values);
            }
        compare(shifted, "geometry-images.lammpstrj");
    });
    auto compare_coordinate_format = [&](const std::string& format) {
        auto converted = frames;
        for (auto& frame : converted) {
            for (size_t c = 0; c < 3; ++c) {
                const std::string axis = std::array<std::string, 3>{"x", "y", "z"}[c];
                const size_t index = frame.column(axis);
                auto bounds = words(frame.header[5 + c]);
                double origin = std::stod(bounds[0]), length = std::stod(bounds[1]) - origin;
                for (auto& atom : frame.atoms) {
                    auto values = words(atom);
                    double value = std::stod(values[index]);
                    if (format != "xu")
                        value = (value - origin) / length;
                    values[index] = decimal(value);
                    atom = join_words(values);
                }
                frame.columns[index] = axis + (format == "xu"   ? "u"
                                               : format == "xs" ? "s"
                                                                : "su");
            }
        }
        compare(converted, "geometry-" + format + ".lammpstrj");
    };
    suite.check("equivalent coordinate formats xu, xs, and xsu", [&] {
        for (const std::string format : {"xu", "xs", "xsu"})
            compare_coordinate_format(format);
    });
    suite.check("atom-column order does not change coordinate/dipole interpretation", [&] {
        auto reversed = frames;
        for (auto& frame : reversed) {
            std::reverse(frame.columns.begin(), frame.columns.end());
            for (auto& atom : frame.atoms) {
                auto values = words(atom);
                std::reverse(values.begin(), values.end());
                atom = join_words(values);
            }
        }
        compare(reversed, "geometry-columns-reversed.lammpstrj");
    });
    suite.check("water geometry atom skipping equals physical removal", [&] {
        auto prefixed = frames;
        for (auto& frame : prefixed)
            frame.atoms.insert(frame.atoms.begin(), frame.atoms.front());
        auto options = config;
        options[suite.section]["skip_initial_atoms"] = 1;
        get_baseline();
        same(evaluate(suite.fixture("geometry-prefix.lammpstrj", encode(prefixed)), options),
             baseline);
    });
    suite.check("two-body geometry agrees with two and four threads", [&] {
        get_baseline();
        for (int threads : {2, 4})
            same(evaluate(original, config, threads), baseline);
    });
}

inline void regression_checks(Suite& suite) {
    const auto manifest = Json::parse(contents(root / "references" / (suite.name + ".json")));
    const bool extended = std::getenv("MB_SPEC_EXTENDED") != nullptr;
    std::map<std::string, fs::path> actual_outputs;
    for (const auto& entry : manifest.at("cases")) {
        if (entry.at("extended").get<bool>() && !extended)
            continue;
        suite.check("trajectory reference: " + entry.at("name").get<std::string>(), [&] {
            require(!entry.value("pending", false),
                    "reference values are pending; generate this suite's references before testing");
            auto input = suite.input(entry.at("input").get<std::string>());
            auto frames = trajectory(input);
            frames.resize(entry.at("frames").get<size_t>());
            const Json config = entry.at("config");
            if (suite.name == "sfg" && config.at("SFG").at("use_ps_dipole").get<bool>() &&
                config.at("SFG").at("dipole_mb") == "1b") {
                // PS replaces permanent dipoles; induced terms are disabled in 1b mode.
                for (auto& frame : frames)
                    for (const std::string column :
                         {"v_dipx_perm", "v_dipy_perm", "v_dipz_perm",
                          "v_dipx_ind", "v_dipy_ind", "v_dipz_ind"})
                        if (std::find(frame.columns.begin(), frame.columns.end(), column) ==
                            frame.columns.end()) {
                            frame.columns.push_back(column);
                            for (auto& atom : frame.atoms)
                                atom += " 0";
                        }
            }
            if (entry.count("cell") && entry.at("cell") == "matrix" &&
                input.extension() == ".lammpstrj") {
                for (auto& frame : frames) {
                    frame.header[4] = "ITEM: BOX BOUNDS xy xz yz pp pp pp";
                    for (size_t axis = 0; axis < 3; ++axis)
                        frame.header[5 + axis] += " 0";
                }
            }
            auto current = suite.fixture("reference-" + input.filename().string(), encode(frames));
            for (const auto& output : entry.at("outputs")) {
                const std::string filename = output.at("file").get<std::string>();
                const std::string stage = fs::path(filename).stem().string();
                const std::string program =
                    (suite.name == "ir" ? "infrared" : suite.name) + "/" + stage;
                const size_t columns = stage == "dipole" || suite.name == "sfg" ? 4
                                       : stage == "polarizability"              ? 7
                                                                                : 2;
                if (entry.count("reuse") && entry.at("reuse").count(filename)) {
                    const std::string previous = entry.at("reuse").at(filename);
                    require(actual_outputs.count(previous),
                            "upstream reference did not pass: " + previous);
                    current = actual_outputs.at(previous);
                } else if ((stage == "dipole" || stage == "polarizability") &&
                           entry.at("cell") == "per-frame-xyz") {
                    // XYZ cells are supplied through JSON; retain each source frame's cell.
                    Table properties;
                    for (size_t i = 0; i < frames.size(); ++i) {
                        const auto& comment = frames[i].header[1];
                        auto start = comment.find("Lattice=\"");
                        require(start != std::string::npos, "missing XYZ lattice metadata");
                        start += 9;
                        auto end = comment.find('"', start);
                        require(end != std::string::npos, "unterminated XYZ lattice metadata");
                        std::istringstream lattice(comment.substr(start, end - start));
                        std::vector<double> cell;
                        for (double value; lattice >> value;)
                            cell.push_back(value);
                        require(cell.size() == 9, "expected nine XYZ lattice components");
                        auto frame_config = config;
                        frame_config[suite.section]["box"] = cell;
                        auto frame_input = suite.fixture("reference-frame.xyz", encode({frames[i]}));
                        auto rows = table(suite.run(program, frame_input, frame_config), columns);
                        require(rows.size() == 1, "expected one XYZ property sample per frame");
                        rows[0][0] = i * config.at(suite.section).at("dump_frequency").get<double>();
                        properties.push_back(rows[0]);
                    }
                    current = suite.fixture(entry.at("name").get<std::string>() + "-" + filename,
                                            serialize(properties));
                } else {
                    current = suite.run(program, current, config);
                }
                const double absolute = program == "infrared/response" ? 1e-70
                                        : program == "raman/spectrum"  ? 1e-35
                                                                       : 1e-12;
                same(table(current, columns), table(output, columns), 2e-5, absolute);
                if ((stage == "dipole" || stage == "polarizability") &&
                    entry.count("molecular")) {
                    const auto& molecular = entry.at("molecular");
                    const std::string molecular_name = molecular.at("filename");
                    same(table(current.parent_path() / molecular_name,
                               stage == "dipole" ? 5 : 11),
                         table(molecular.at("table"), stage == "dipole" ? 5 : 11));
                }
                actual_outputs[entry.at("name").get<std::string>() + "/" + filename] = current;
            }
        });
    }
    if (manifest.count("fresnel")) {
        const auto& entry = manifest.at("fresnel");
        suite.check("Fresnel reference: " + entry.at("name").get<std::string>(), [&] {
            auto config = entry.at("config");
            for (const std::string component : {"xxz", "zzz", "xzx"}) {
                const std::string source = entry.at("components").at(component);
                require(actual_outputs.count(source), "component reference did not pass: " + source);
                // Pass the full grid to verify that Fresnel skips zero wavenumber.
                config[component + "_file"] =
                    suite.fixture("fresnel-" + component + ".dat",
                                  contents(actual_outputs.at(source))).string();
            }
            for (const std::string key :
                 {"water_refractive_index_n_file", "water_refractive_index_k_file"}) {
                const std::string filename = config.at(key);
                const auto path = root / ".." / "sfg" / "fresnel_factors" / filename;
                config[key] = path.string();
            }
            const auto result = suite.run("sfg/fresnel", {}, config);
            for (const auto& output : entry.at("outputs")) {
                auto rows = table(result.parent_path() / output.at("file").get<std::string>(), 2);
                require(std::any_of(rows.begin(), rows.end(),
                                    [](const auto& row) { return row[1] > 0; }),
                        "Fresnel intensity must not be entirely zero");
                same(rows, table(output, 2), 2e-5, 1e-12);
            }
        });
    }
}
} // namespace test
#endif
