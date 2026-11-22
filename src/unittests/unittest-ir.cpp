#include "test_support.h"
#include "ltp2011.h"
#include "json_options.h"

using namespace test;

int main() {
    Suite suite("ir", "INFRARED");
    const auto bulk = trajectory(suite.input("bulk-water.lammpstrj"));
    suite.check("bulk internal water dipoles agree between XYZ and LAMMPS", [&] {
        auto config = suite.config();
        config["INFRARED"]["DMS"] = true;
        auto xyz = table(suite.run("infrared/dipole", suite.input("bulk-water.xyz"), config), 4);
        auto lammps =
            table(suite.run("infrared/dipole", suite.input("bulk-water.lammpstrj"), config), 4);
        check_times(xyz, 10);
        same(xyz, lammps, 1e-7, 1e-7);
    });
    suite.check("one-body dipole rotates as a vector and is translation invariant", [&] {
        double xyz[9] = {0, 0, 0, .9572, 0, 0, -.239987, .927297, 0};
        double rotated[9], original[3], result[3];
        for (size_t i = 0; i < 3; ++i) {
            rotated[3 * i] = -xyz[3 * i + 1] + 4;
            rotated[3 * i + 1] = xyz[3 * i] - 7;
            rotated[3 * i + 2] = xyz[3 * i + 2] + 3;
        }
        h2o::ltp2011 model;
        model.cartesian(1, xyz, original);
        model.cartesian(1, rotated, result);
        near(result[0], -original[1]);
        near(result[1], original[0]);
        near(result[2], original[2]);
    });
    suite.check("molecular dipoles sum to total; one and two threads agree", [&] {
        auto config = suite.config();
        config["INFRARED"]["DMS"] = true;
        config["INFRARED"]["write_molecular_properties"] = true;
        auto output = suite.run("infrared/dipole", suite.input("bulk-water.xyz"), config);
        auto molecules = table(output.parent_path() / "molecular_dipoles.dat", 5);
        require(molecules.size() == 2560, "expected 256 waters times ten frames");
        Table sum(10, std::vector<double>(4, 0));
        for (size_t i = 0; i < 10; ++i)
            sum[i][0] = i * dt;
        for (const auto& row : molecules) {
            const size_t i = std::lround(row[1] / dt);
            require(i < 10, "invalid molecular time");
            for (size_t c = 0; c < 3; ++c)
                sum[i][c + 1] += row[c + 2];
        }
        same(table(output, 4), sum);
        same(table(output, 4),
             table(suite.run("infrared/dipole", suite.input("bulk-water.xyz"), config, "", 2), 4));
    });

    auto small = suite.fixture("two-waters.lammpstrj", encode(nearby_waters(bulk)));
    std::map<std::string, Table> modes;
    suite.check("two-body switch weights 0, 0.5, and 1", [&] {
        for (const auto& radii : std::vector<std::array<double, 3>>{{1, 3, 0}, {2, 4, .5}, {3, 4, 1}}) {
            auto config = suite.config();
            for (const std::string mode : {"1b", "1b+2b"}) {
                if (!modes.count(mode)) {
                    config["INFRARED"]["dipole_mb"] = mode;
                    auto rows = table(suite.run("infrared/dipole", small, config), 4);
                    check_times(rows, 10);
                    modes[mode] = rows;
                }
            }
            config["INFRARED"]["dipole_mb"] = "1b+2b";
            config["INFRARED"]["cutoff_1"] = radii[0];
            config["INFRARED"]["cutoff_2"] = radii[1];
            auto expected = modes.at("1b");
            double correction_norm = 0;
            for (size_t i = 0; i < 10; ++i)
                for (size_t c = 1; c < 4; ++c) {
                    double delta = modes.at("1b+2b")[i][c] - modes.at("1b")[i][c];
                    correction_norm += std::abs(delta);
                    expected[i][c] += radii[2] * delta;
                }
            require(correction_norm > 1e-6, "fixture must exercise a nonzero two-body correction");
            same(table(suite.run("infrared/dipole", small, config), 4), expected);
        }
    });
    suite.check("leading nonwater record can be skipped without changing water dipoles", [&] {
        auto frames = nearby_waters(bulk);
        for (auto& frame : frames) {
            auto row = words(frame.atoms.front());
            for (auto& value : row)
                value = "0";
            row[frame.column("v_dipx_perm")] = "100";
            std::string atom;
            for (const auto& value : row)
                atom += value + " ";
            frame.atoms.insert(frame.atoms.begin(), atom);
        }
        auto with_ion = suite.fixture("synthetic-leading-ion.lammpstrj", encode(frames));
        auto config = suite.config();
        auto baseline = table(suite.run("infrared/dipole", small, config), 4);
        config["INFRARED"]["skip_initial_atoms"] = 1;
        same(table(suite.run("infrared/dipole", with_ion, config), 4), baseline);
    });
    geometry_checks(suite, bulk);
    json_option_checks(suite, bulk);
    regression_checks(suite);
    return suite.finish();
}
