#include "test_support.h"
#include "avila.h"
#include "json_options.h"

using namespace test;

static Table weighted_self_response(const std::vector<Frame>& frames, const Json& options) {
    struct Molecule {
        double dipole[3], alpha[9], weight;
    };
    std::vector<std::vector<Molecule>> properties;
    for (const auto& frame : frames) {
        std::vector<Molecule> molecules;
        for (size_t m = 0; m < 2; ++m) {
            Molecule molecule{};
            double xyz[9];
            for (size_t atom = 0; atom < 3; ++atom) {
                auto row = words(frame.atoms[3 * m + atom]);
                for (size_t c = 0; c < 3; ++c)
                    xyz[3 * atom + c] =
                        std::stod(row[frame.column(std::array<std::string, 3>{"x", "y", "z"}[c])]);
            }
            x2o::avila::alpha(xyz, molecule.alpha);
            auto oxygen = words(frame.atoms[3 * m]);
            for (size_t c = 0; c < 3; ++c)
                molecule.dipole[c] = std::stod(oxygen[frame.column(
                    "v_dip" + std::array<std::string, 3>{"x", "y", "z"}[c] + "_perm")]);
            const std::string surface = options.at("surface");
            const double lower = options.at("zc1"), upper = options.at("zc2"), z = xyz[2];
            double weight =
                std::abs(z) <= lower ? 0
                : std::abs(z) >= upper
                    ? 1
                    : std::pow(std::sin(M_PI * (std::abs(z) - lower) / (2 * (upper - lower))), 2);
            if ((surface == "positive" && z < 0) || (surface == "negative" && z > 0))
                weight = 0;
            molecule.weight = std::copysign(weight, z);
            molecules.push_back(molecule);
        }
        properties.push_back(molecules);
    }
    const std::string component = options.at("ijk");
    const std::array<int, 4> indices = component == "xxz"   ? std::array<int, 4>{2, 0, 2, 4}
                                       : component == "zzz" ? std::array<int, 4>{2, 8, 2, 8}
                                       : component == "xzx" ? std::array<int, 4>{0, 2, 1, 5}
                                                            : std::array<int, 4>{0, 6, 1, 7};
    Table result;
    for (size_t lag = 0; lag < 2; ++lag) {
        double sum = 0;
        for (size_t i = 0; i + lag < 2; ++i)
            for (size_t m = 0; m < 2; ++m) {
                const auto& a = properties[i][m];
                const auto& b = properties[i + lag][m];
                double weight =
                    a.weight *
                    std::pow(options.at("time_dependent_switch").get<bool>() ? b.weight : a.weight,
                             2);
                sum += weight * .5 *
                       (a.dipole[indices[0]] * b.alpha[indices[1]] +
                        a.dipole[indices[2]] * b.alpha[indices[3]]);
            }
        result.push_back({lag * dt, sum / (2 - lag), sum, double(2 - lag)});
    }
    return result;
}

static void sfg_option_checks(Suite& suite, const std::vector<Frame>& air) {
    auto frames = nearby_waters(air);
    frames.resize(2);
    auto close = suite.fixture("sfg-options-neighbors.lammpstrj", encode(frames));
    for (size_t i = 0; i < 2; ++i)
        for (size_t atom = 0; atom < 6; ++atom) {
            auto values = words(frames[i].atoms[atom]);
            double target = (atom < 3 ? 1 : -1) * (i ? 13.5 : 12.);
            size_t c = frames[i].column("z");
            values[c] = decimal(std::stod(values[c]) + target - 12.);
            frames[i].atoms[atom] = join_words(values);
        }
    auto input = suite.fixture("sfg-options-interfaces.lammpstrj", encode(frames));
    auto config = suite.config();
    auto& options = config["SFG"];
    options["int_max"] = 2 * dt;
    options["center_slab"] = false;
    options["zc1"] = 11.;
    options["zc2"] = 13.;
    options["dipole_mb"] = "1b";
    options["polarizability_mb"] = "1b";
    options["include_cross_correlation"] = false;
    options["include_autocorrelation"] = true;
    options["use_ps_dipole"] = false;
    for (const std::string ijk : {"xxz", "zzz", "xzx", "zxx"}) {
        suite.check("two-frame SFG " + ijk + ": all surfaces and fixed/moving weights", [&] {
            for (const std::string surface : {"positive", "negative", "both"})
                for (bool moving : {false, true}) {
                    auto settings = config;
                    settings["SFG"]["ijk"] = ijk;
                    settings["SFG"]["surface"] = surface;
                    settings["SFG"]["time_dependent_switch"] = moving;
                    same(table(suite.run("sfg/response", input, settings), 4),
                         weighted_self_response(frames, settings["SFG"]));
                }
        });
    }
    suite.check("two-frame interface weights: partial, full, and excluded", [&] {
        for (const auto& bounds :
             std::vector<std::array<double, 2>>{{10, 13}, {11, 14}, {0, 1}, {20, 21}}) {
            auto settings = config;
            settings["SFG"]["zc1"] = bounds[0];
            settings["SFG"]["zc2"] = bounds[1];
            same(table(suite.run("sfg/response", input, settings), 4),
                 weighted_self_response(frames, settings["SFG"]));
        }
    });
    suite.check("neighbors outside or exactly at cutoff give zero cross response", [&] {
        auto settings = config;
        settings["SFG"]["zc1"] = 0.;
        settings["SFG"]["zc2"] = .1;
        settings["SFG"]["include_autocorrelation"] = false;
        settings["SFG"]["include_cross_correlation"] = true;
        settings["SFG"]["cross_correlation_sampling_stride"] = 1;
        for (double cutoff : {2., 3.}) {
            settings["SFG"]["cross_correlation_cutoff"] = cutoff;
            for (const auto& row : table(suite.run("sfg/response", close, settings), 4)) {
                near(row[1], 0.);
                near(row[2], 0.);
            }
        }
    });
    suite.check(
        "two-frame slab centering removes a uniform z translation and half-cell shift", [&] {
            auto settings = config;
            settings["SFG"]["zc1"] = 0.;
            settings["SFG"]["zc2"] = 1.;
            settings["SFG"]["center_slab"] = true;
            auto expected = table(suite.run("sfg/response", close, settings), 4);
            auto translated = nearby_waters(air);
            translated.resize(2);
            for (auto& frame : translated)
                for (auto& atom : frame.atoms) {
                    auto values = words(atom);
                    auto c = frame.column("z");
                    values[c] = decimal(std::stod(values[c]) + 3);
                    atom = join_words(values);
                }
            auto shifted = suite.fixture("sfg-options-translated.lammpstrj", encode(translated));
            same(table(suite.run("sfg/response", shifted, settings), 4), expected);
            settings["SFG"]["shift_slab_by_half_cell"] = true;
            same(table(suite.run("sfg/response", shifted, settings), 4), expected);
        });
}

int main() {
    Suite suite("sfg", "SFG");
    const auto air = trajectory(suite.input("air-water.lammpstrj"));
    const auto graphene = trajectory(suite.input("graphene-water.lammpstrj"));
    suite.check("one-body response agrees with one and two threads for all components", [&] {
        for (const std::string component : {"xxz", "zzz", "xzx", "zxx"}) {
            auto config = suite.config();
            config["SFG"]["ijk"] = component;
            auto baseline = table(
                suite.run("sfg/response", suite.input("air-water.lammpstrj"), config), 4);
            auto parallel = table(
                suite.run("sfg/response", suite.input("air-water.lammpstrj"), config, "", 2), 4);
            same(parallel, baseline);
        }
    });
    suite.check("skipping 180 leading carbons matches a water-only graphene fixture", [&] {
        auto water_only =
            suite.fixture("graphene-water-only.lammpstrj", encode(subset(graphene, 180, 1536)));
        auto config = suite.config();
        config["SFG"]["skip_initial_atoms"] = 180;
        auto skipped =
            table(suite.run("sfg/response", suite.input("graphene-water.lammpstrj"), config), 4);
        config["SFG"]["skip_initial_atoms"] = 0;
        same(table(suite.run("sfg/response", water_only, config), 4), skipped);
    });
    geometry_checks(suite, air);
    json_option_checks(suite, air);
    sfg_option_checks(suite, air);
    regression_checks(suite);
    return suite.finish();
}
