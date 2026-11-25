#include "test_support.h"

namespace {

using namespace test;

Json historical_config(const std::string& workflow)
{
    return Json::parse(contents(root / "historical/inputs" / (workflow + ".json")));
}

void check_usage(Suite& suite)
{
    for (const auto* executable : {"infrared/dipole", "infrared/response", "infrared/spectrum",
                                   "raman/polarizability", "raman/response", "raman/spectrum",
                                   "power/response", "power/spectrum", "sfg/response", "sfg/spectrum",
                                   "sfg/fresnel"}) {
        suite.check(std::string(executable) + " accepts help and rejects incorrect argument counts", [&] {
            for (const auto* flag : {"-h", "--help"}) {
                const auto output = suite.run_arguments(executable, {flag});
                require(contents(output).find("Usage: mbspec-") == 0, "missing help on stdout");
                require(contents(output.parent_path() / "stderr.log").empty(), "unexpected help diagnostic");
            }
            for (const auto& arguments : {std::vector<std::string>{}, std::vector<std::string>{"a", "b", "c"}}) {
                const auto output = suite.run_arguments(executable, arguments, "Usage: mbspec-");
                require(contents(output).empty(), "unexpected output with incorrect arguments");
            }
        });
    }
}

void check_configuration(Suite& suite)
{
    suite.check("workflow JSON rejects wrong sections and unsupported options", [&] {
        const auto input = root / "historical/references/ir/dipole.dat";
        const auto program = "infrared/response";
        auto config = historical_config("ir");
        const auto expected = contents(suite.run(program, input, config));
        config["POWER"] = historical_config("power").at("POWER");
        require(contents(suite.run(program, input, config)) == expected, "combined config changed output");
        config["INFRARED"]["subtract_mean_diplole"] = true;
        suite.run(program, input, config, "unknown option INFRARED.subtract_mean_diplole");
        config["INFRARED"].erase("subtract_mean_diplole");
        config["INFRARED"]["subtract_mean_velocity"] = true;
        suite.run(program, input, config, "unknown option INFRARED.subtract_mean_velocity");
        config.erase("INFRARED");
        suite.run(program, input, config, "missing JSON workflow section INFRARED");
    });
}

void check_response_tables(Suite& suite, const std::string& workflow)
{
    const auto program = (workflow == "ir" ? "infrared" : workflow) + "/spectrum";
    const auto section = workflow == "ir" ? "INFRARED" : workflow == "raman" ? "RAMAN" : "POWER";
    auto config = historical_config(workflow);
    config[section]["spectrum_transform_time_ps"] = 0.5;

    std::ostringstream rounded;
    if (workflow == "power")
        rounded << std::setprecision(6);
    else
        rounded << std::scientific << std::setprecision(8);
    for (int i = 0; i < 9; ++i)
        rounded << i * 0.123456789123 << ' ' << 9 - i << '\n';

    suite.check(workflow + " accepts timestamps rounded by the response writer", [&] {
        const auto input = suite.fixture(workflow + "-rounded.dat", rounded.str());
        require(!table(suite.run(program, input, config), 2).empty(), "empty spectrum");
    });

    suite.check(workflow + " rejects malformed, nonfinite and uneven response rows", [&] {
        const auto bad_row = workflow + "-invalid.dat";
        for (const auto* tail : {"malformed\n", "2 nan\n", "2 1e999\n", "2 1 2\n", "2\n", "2 1garbage\n"}) {
            const auto input = suite.fixture(bad_row, rounded.str() + tail);
            suite.run(program, input, config, bad_row + ":10: expected exactly two finite numbers");
        }
        const auto input = suite.fixture(bad_row, "0 4\n0.125 3\n0.3 2\n0.375 1\n");
        suite.run(program, input, config, bad_row + ":3: response times must be evenly spaced");
        suite.fixture(bad_row, "0 4\n0.125 3\n0.125 2\n0.375 1\n");
        suite.run(program, input, config, bad_row + ":3: response times must be finite and strictly increasing");
    });
}

void check_fresnel(Suite& suite)
{
    auto config = historical_config("fresnel");
    const auto aligned = "500 1 2\n1000 1 2\n";
    config["xxz_file"] = suite.fixture("xxz.dat", aligned).string();
    config["xzx_file"] = suite.fixture("xzx.dat", aligned).string();
    config["zzz_file"] = suite.fixture("zzz.dat", "500 1 2 5\n1000 1 2 5\n").string();
    for (const auto* key : {"water_refractive_index_n_file", "water_refractive_index_k_file"})
        config[key] = (root.parent_path() / "sfg" / config.at(key).get<std::string>()).string();

    suite.check("Fresnel accepts matching grids with three or four columns", [&] {
        const auto output = suite.run("sfg/fresnel", {}, config);
        require(table(output.parent_path() / config.at("ppp_output_file").get<std::string>(), 2).size() == 2,
                "expected two Fresnel samples");
    });
    suite.check("Fresnel rejects shifted or malformed grids before writing output", [&] {
        suite.fixture("xxz.dat", "519 1 2\n1019 1 2\n");
        const auto output = suite.run("sfg/fresnel", {}, config, "xxz.dat:1: wavenumber grid differs");
        require(!fs::exists(output.parent_path() / config.at("ppp_output_file").get<std::string>()),
                "invalid Fresnel grid produced an output file");
        suite.fixture("xxz.dat", "500 1 2 garbage\n1000 1 2\n");
        suite.run("sfg/fresnel", {}, config, "xxz.dat:1: expected three or four finite spectrum columns");
        config["redhsift_cm1"] = 0.0;
        suite.run("sfg/fresnel", {}, config, "unknown Fresnel option redhsift_cm1");
    });
}

} // namespace

int main()
{
    Suite suite("inputs", "");
    check_usage(suite);
    check_configuration(suite);
    for (const auto* workflow : {"ir", "raman", "power"})
        check_response_tables(suite, workflow);
    check_fresnel(suite);
    return suite.finish();
}
