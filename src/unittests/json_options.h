#ifndef MB_SPEC_JSON_OPTIONS_H
#define MB_SPEC_JSON_OPTIONS_H

#include "test_support.h"

namespace test {

inline void json_option_checks(Suite& suite, const std::vector<Frame>& source) {
    auto frames = suite.name == "power" ? source : nearby_waters(source);
    frames.resize(2);
    auto input = suite.fixture("options-two-frames.lammpstrj", encode(frames));
    auto config = suite.config();
    auto& options = config[suite.section];
    options["int_max"] = 2 * dt;
    options["spectrum_transform_time_ps"] = 2 * dt;
    if (suite.name == "ir")
        options["DMS"] = true;
    if (suite.name == "sfg") {
        options["center_slab"] = false;
        options["zc1"] = 0.;
        options["zc2"] = .1;
    }
    if (suite.name == "power")
        options["skip_initial_atoms"] = 1;
    const size_t columns = suite.name == "raman" ? 7 : suite.name == "power" ? 2 : 4;
    auto evaluate = [&](const Json& settings, const fs::path& path) {
        return table(suite.run(first_stage(suite), path, settings), columns);
    };
    auto baseline = evaluate(config, input);
    if (suite.name == "sfg") {
        suite.check("two-body switch weights 0, 0.5, and 1", [&] {
            auto settings = config;
            settings[suite.section]["dipole_mb"] = "1b+2b";
            auto full = evaluate(settings, input);
            for (const auto& radii :
                 std::vector<std::array<double, 3>>{{1, 3, 0}, {2, 4, .5}, {3, 4, 1}}) {
                settings[suite.section]["cutoff_1"] = radii[0];
                settings[suite.section]["cutoff_2"] = radii[1];
                auto expected = baseline;
                for (size_t i = 0; i < expected.size(); ++i)
                    for (size_t c = 1; c < 3; ++c)
                        expected[i][c] += radii[2] * (full[i][c] - baseline[i][c]);
                same(evaluate(settings, input), expected);
            }
        });
    }

    suite.check("dump_frequency changes the time axis", [&] {
        auto changed = config;
        changed[suite.section]["dump_frequency"] = 2 * dt;
        changed[suite.section]["int_max"] = 4 * dt;
        auto expected = baseline;
        for (auto& row : expected)
            row[0] *= 2;
        same(evaluate(changed, input), expected);
    });
    if (suite.name == "power") {
        suite.check("velocity correlation is independent of box and cell convention", [&] {
            auto settings = config;
            settings[suite.section]["box"] = {0, 50, 0, 50, 0, 50};
            same(evaluate(settings, input), baseline);
            auto velocity_xyz = trajectory(suite.input("f-water.xyz"));
            velocity_xyz.resize(2);
            auto xyz = suite.fixture("options-two-frames.xyz", encode(velocity_xyz));
            auto expected = evaluate(config, xyz);
            same(expected, baseline);
            for (int imcon : {1, 2, 3}) {
                settings[suite.section]["imcon"] = imcon;
                settings[suite.section]["box"] = {20, 0, 0, 2, 20, 0, 1, .5, 20};
                same(evaluate(settings, xyz), expected);
            }
        });
    }

    if (suite.name == "sfg")
        return;

    // Two response samples suffice except for Raman's explicit four-point reader limit.
    Table response = {{0, 4}, {dt, 2}};
    if (suite.name == "ir")
        for (auto& row : response)
            row[1] *= 1e-59;
    if (suite.name == "raman") {
        response.push_back({2 * dt, -1});
        response.push_back({3 * dt, .5});
    }
    auto spectrum_input = suite.fixture("options-response.dat", serialize(response));
    auto spectral_config = config;
    spectral_config[suite.section]["int_max"] = response.size() * dt;
    spectral_config[suite.section]["spectrum_transform_time_ps"] = response.size() * dt;
    const std::string spectrum_program =
        (suite.name == "ir" ? "infrared" : suite.name) + "/spectrum";
    auto spectrum = [&](const Json& settings) {
        return table(suite.run(spectrum_program, spectrum_input, settings), 2);
    };
    auto original_spectrum = spectrum(spectral_config);
    if (suite.name == "ir" || suite.name == "raman") {
        suite.check("temperature dependence of the spectrum", [&] {
            auto changed = spectral_config;
            changed[suite.section]["temperature"] = 600.;
            auto actual = spectrum(changed), expected = original_spectrum;
            for (size_t i = 0; i < expected.size(); ++i) {
                if (suite.name == "ir")
                    expected[i][1] *= .5;
                if (suite.name == "raman" && i) {
                    double x = 1.054571726e-34 * 2 * M_PI * expected[i][0] * 29979245800. /
                               (2 * 1.3806488e-23 * 300);
                    expected[i][1] *= std::tanh(x / 2) / std::tanh(x);
                }
            }
            same(actual, expected, 3e-5, suite.name == "raman" ? 1e-35 : 1e-12);
        });
    }
    if (suite.name == "ir" || suite.name == "raman") {
        const std::string label = suite.name == "raman"
                                      ? "Raman spectrum is independent of volume and box representation"
                                      : "box volume dependence and matrix representation";
        suite.check(label, [&] {
            auto changed = spectral_config;
            auto box = changed[suite.section]["box"].get<std::vector<double>>();
            changed[suite.section]["box"] = {box[1] - box[0], 0, 0, 0, box[3] - box[2], 0, 0, 0,
                                             box[5] - box[4]};
            same(spectrum(changed), original_spectrum);
            box[1] = box[0] + 2 * (box[1] - box[0]);
            changed[suite.section]["box"] = box;
            auto expected = original_spectrum;
            if (suite.name == "ir")
                for (auto& row : expected)
                    row[1] *= .5;
            same(spectrum(changed), expected);
        });
    }
}
} // namespace test
#endif
