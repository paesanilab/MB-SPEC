#include "test_support.h"

namespace {

using namespace test;

Json historical_config(const std::string& workflow)
{
    return Json::parse(contents(root / "historical/inputs" / (workflow + ".json")));
}

void check_usage(Suite& suite)
{
    struct Command {
        const char* executable;
        const char* synopsis;
        bool config_only = false;
    };
    const Command commands[] = {
        {"infrared/dipole", "mbspec-ir-dipole <trajectory> <mb-spec.json>"},
        {"infrared/response", "mbspec-ir-response <dipole_file> <mb-spec.json>"},
        {"infrared/spectrum", "mbspec-ir-spectrum <response.dat> <mb-spec.json>"},
        {"raman/polarizability", "mbspec-raman-polarizability <trajectory> <mb-spec.json>"},
        {"raman/response", "mbspec-raman-response <polarizability_file> <mb-spec.json>"},
        {"raman/spectrum", "mbspec-raman-spectrum <response.dat> <mb-spec.json>"},
        {"power/response", "mbspec-power-response <velocity_trajectory> <mb-spec.json>"},
        {"power/spectrum", "mbspec-power-spectrum <response.dat> <mb-spec.json>"},
        {"sfg/response", "mbspec-sfg-response <trajectory.lammpstrj> <mb-spec.json>"},
        {"sfg/spectrum", "mbspec-sfg-spectrum <response.dat> <mb-spec.json>"},
        {"sfg/fresnel", "mbspec-sfg-fresnel <fresnel.json>", true},
    };
    for (const auto& command : commands) {
        const auto executable = command.executable;
        suite.check(std::string(executable) + " accepts help and rejects incorrect argument counts", [&] {
            const auto usage = std::string("Usage: ") + command.synopsis + '\n';
            for (const auto* flag : {"-h", "--help"}) {
                const auto output = suite.run_arguments(executable, {flag});
                require(contents(output) == usage, "wrong usage on stdout");
                require(contents(output.parent_path() / "stderr.log").empty(), "unexpected help diagnostic");
            }
            const auto check_bad_count = [&](const std::vector<std::string>& arguments) {
                const std::string expected_count = command.config_only ? "1 positional argument" : "2 positional arguments";
                const auto diagnostic = "Expected " + expected_count + ", got " +
                                        std::to_string(arguments.size()) + ".\n" + usage;
                const auto output = suite.run_arguments(executable, arguments, diagnostic);
                require(contents(output).empty(), "unexpected output with incorrect arguments");
                require(contents(output.parent_path() / "stderr.log") == diagnostic, "wrong argument diagnostic");
            };
            check_bad_count({});
            if (command.config_only) {
                check_bad_count({"a", "b"});
            } else {
                check_bad_count({"a"});
                check_bad_count({"a", "b", "c"});
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

    suite.check(workflow + " stops reading at the exclusive cutoff", [&] {
        auto cutoff_config = config;
        cutoff_config[section]["int_max"] = 1.5;
        cutoff_config[section]["spectrum_transform_time_ps"] = 1.5;
        const std::string prefix = "0 6\n0.25 5\n0.5 4\n0.75 3\n1 2\n1.25 1\n";
        const auto input = suite.fixture(workflow + "-cutoff.dat", prefix);
        const auto expected = contents(suite.run(program, input, cutoff_config));
        for (const auto* tail : {"1.5 0\n1.75 0\n2.25 0\nmalformed\n", "1.5 unused\n"}) {
            suite.fixture(workflow + "-cutoff.dat", prefix + tail);
            require(contents(suite.run(program, input, cutoff_config)) == expected,
                    "response rows at or above the cutoff changed the spectrum");
        }
    });

    suite.check(workflow + " rejects malformed, nonfinite and uneven response rows", [&] {
        const auto bad_row = workflow + "-invalid.dat";
        for (const auto* tail : {"malformed\n", "nan 1\n", "1.125 nan\n", "1.125 1e999\n",
                                 "1.125 1 2\n", "1.125\n", "1.125 1garbage\n"}) {
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
