#include "test_support.h"
#include "json_options.h"

using namespace test;

int main() {
    Suite suite("power", "POWER");
    const auto xyz = trajectory(suite.input("f-water.xyz"));
    const auto lammps = trajectory(suite.input("f-water.lammpstrj"));
    suite.check("fluoride skipping equals removal in XYZ and LAMMPS; threads agree", [&] {
        for (const auto& filename : {"f-water.xyz", "f-water.lammpstrj"}) {
            auto config = suite.config();
            config["POWER"]["skip_initial_atoms"] = 1;
            auto baseline = table(suite.run("power/response", suite.input(filename), config), 2);
            const auto& frames = std::string(filename) == "f-water.xyz" ? xyz : lammps;
            auto water = suite.fixture("water-only-" + std::string(filename),
                                       encode(subset(frames, 1, 831)));
            config["POWER"]["skip_initial_atoms"] = 0;
            for (int threads : {1, 2, 4})
                same(table(suite.run("power/response", water, config, "", threads), 2), baseline);
        }
    });
    const Table response = {{0, 4}, {dt, 2}, {2 * dt, -1}, {3 * dt, .5}};
    const auto response_path = suite.fixture("normalization-response.dat", serialize(response));
    suite.check("normalization is unchanged by response amplitude", [&] {
        auto scaled = response;
        for (auto& row : scaled)
            row[1] *= 7;
        auto input = suite.fixture("scaled-response.dat", serialize(scaled));
        same(table(suite.run("power/spectrum", input, suite.config()), 2),
             table(suite.run("power/spectrum", response_path, suite.config()), 2));
    });
    suite.check("velocity column order is arbitrary", [&] {
        auto reversed = lammps;
        for (auto& frame : reversed) {
            std::reverse(frame.columns.begin(), frame.columns.end());
            for (auto& atom : frame.atoms) {
                auto values = words(atom);
                std::reverse(values.begin(), values.end());
                atom = join_words(values);
            }
        }
        auto input = suite.fixture("velocities-columns-reversed.lammpstrj", encode(reversed));
        same(table(suite.run("power/response", input, suite.config()), 2),
             table(suite.run("power/response", suite.input("f-water.lammpstrj"), suite.config()),
                   2));
    });
    json_option_checks(suite, lammps);
    regression_checks(suite);
    return suite.finish();
}
