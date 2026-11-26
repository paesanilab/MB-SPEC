#ifndef MB_SPEC_CLI_H
#define MB_SPEC_CLI_H

#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <optional>
#include <string_view>

namespace mbspec {

[[nodiscard]] inline std::optional<int> check_arguments(
    int argc, char** argv, std::initializer_list<std::string_view> positional_names)
{
    const bool help = argc == 2 && (std::string_view(argv[1]) == "--help" ||
                                   std::string_view(argv[1]) == "-h");
    const int actual_count = argc - 1;
    const auto expected_count = static_cast<int>(positional_names.size());
    if (!help && actual_count == expected_count)
        return std::nullopt;

    std::string_view program(argv[0]);
    if (const auto slash = program.find_last_of('/'); slash != std::string_view::npos)
        program.remove_prefix(slash + 1);
    auto& output = help ? std::cout : std::cerr;
    if (!help) {
        output << "Expected " << expected_count << " positional argument"
               << (expected_count == 1 ? "" : "s") << ", got " << actual_count << ".\n";
    }
    output << "Usage: " << program;
    for (const auto name : positional_names)
        output << " <" << name << '>';
    output << '\n';
    return help ? EXIT_SUCCESS : EXIT_FAILURE;
}

} // namespace mbspec

#endif
