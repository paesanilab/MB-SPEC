#include "test_support.h"
#include "avila.h"
#include "json_options.h"

using namespace test;

int main() {
    Suite suite("raman", "RAMAN");
    const auto bulk = trajectory(suite.input("bulk-water.lammpstrj"));
    suite.check("bulk XYZ and LAMMPS one-body tensors agree", [&] {
        auto config = suite.config();
        auto xyz =
            table(suite.run("raman/polarizability", suite.input("bulk-water.xyz"), config), 7);
        auto lammps = table(
            suite.run("raman/polarizability", suite.input("bulk-water.lammpstrj"), config), 7);
        check_times(xyz, 10);
        same(xyz, lammps, 1e-8, 1e-8);
    });
    suite.check("molecular tensors sum to total and remain symmetric", [&] {
        auto config = suite.config();
        auto disabled = suite.run("raman/polarizability", suite.input("bulk-water.xyz"), config);
        require(!fs::exists(disabled.parent_path() / "molecular_polarizabilities.dat"),
                "disabled molecular output was written");
        config["RAMAN"]["write_molecular_properties"] = true;
        auto output = suite.run("raman/polarizability", suite.input("bulk-water.xyz"), config);
        same(table(output, 7), table(disabled, 7));
        auto molecules = table(output.parent_path() / "molecular_polarizabilities.dat", 11);
        require(molecules.size() == 2560, "expected 256 waters times ten frames");
        Table sum(10, std::vector<double>(7, 0));
        const std::array<size_t, 6> components = {0, 1, 2, 4, 5, 8};
        for (size_t i = 0; i < 10; ++i)
            sum[i][0] = i * dt;
        for (const auto& row : molecules) {
            size_t i = std::lround(row[1] / dt);
            require(i < 10, "invalid tensor time");
            near(row[3], row[5]);
            near(row[4], row[8]);
            near(row[7], row[9]);
            for (size_t c = 0; c < 6; ++c)
                sum[i][c + 1] += row[components[c] + 2];
        }
        same(table(output, 7), sum);
        same(table(output, 7),
             table(suite.run("raman/polarizability", suite.input("bulk-water.xyz"), config, "", 2),
                   7));
    });
    suite.check("Avila tensor transforms under rotation and translation", [&] {
        double xyz[9] = {0, 0, 0, .9572, 0, 0, -.239987, .927297, 0};
        double rotated[9], a[9], b[9];
        for (size_t i = 0; i < 3; ++i) {
            rotated[3 * i] = -xyz[3 * i + 1] + 3;
            rotated[3 * i + 1] = xyz[3 * i] - 5;
            rotated[3 * i + 2] = xyz[3 * i + 2] + 8;
        }
        x2o::avila::alpha(xyz, a);
        x2o::avila::alpha(rotated, b);
        const double r[9] = {0, -1, 0, 1, 0, 0, 0, 0, 1};
        for (size_t i = 0; i < 3; ++i)
            for (size_t j = 0; j < 3; ++j) {
                double expected = 0;
                for (size_t k = 0; k < 3; ++k)
                    for (size_t l = 0; l < 3; ++l)
                        expected += r[3 * i + k] * a[3 * k + l] * r[3 * j + l];
                near(b[3 * i + j], expected);
            }
    });
    auto small = suite.fixture("two-waters.lammpstrj", encode(nearby_waters(bulk)));
    std::map<std::string, Table> modes;
    suite.check("two-body polarizability switch weights 0, 0.5, and 1", [&] {
        for (const auto& radii : std::vector<std::array<double, 3>>{{1, 3, 0}, {2, 4, .5}, {3, 4, 1}}) {
            auto config = suite.config();
            for (const std::string mode : {"1b", "1b+2b"}) {
                if (!modes.count(mode)) {
                    config["RAMAN"]["polarizability_mb"] = mode;
                    auto rows = table(suite.run("raman/polarizability", small, config), 7);
                    check_times(rows, 10);
                    modes[mode] = rows;
                }
            }
            config["RAMAN"]["polarizability_mb"] = "1b+2b";
            config["RAMAN"]["cutoff_1"] = radii[0];
            config["RAMAN"]["cutoff_2"] = radii[1];
            auto expected = modes.at("1b");
            double correction_norm = 0;
            for (size_t i = 0; i < 10; ++i)
                for (size_t c = 1; c < 7; ++c) {
                    double delta = modes.at("1b+2b")[i][c] - modes.at("1b")[i][c];
                    correction_norm += std::abs(delta);
                    expected[i][c] += radii[2] * delta;
                }
            require(correction_norm > 1e-6, "fixture must exercise a nonzero two-body correction");
            same(table(suite.run("raman/polarizability", small, config), 7), expected);
        }
    });
    geometry_checks(suite, bulk);
    json_option_checks(suite, bulk);
    regression_checks(suite);
    return suite.finish();
}
