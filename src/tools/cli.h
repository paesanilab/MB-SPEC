#ifndef MB_SPEC_CLI_H
#define MB_SPEC_CLI_H

#include <cstdlib>
#include <iostream>
#include <optional>
#include <string_view>

namespace mbspec {

inline std::optional<int> check_arguments(int argc, char** argv, int expected,
                                          std::string_view arguments)
{
    const bool help = argc == 2 && (std::string_view(argv[1]) == "--help" ||
                                   std::string_view(argv[1]) == "-h");
    if (!help && argc == expected)
        return std::nullopt;

    std::string_view program(argv[0]);
    if (const auto slash = program.find_last_of('/'); slash != std::string_view::npos)
        program.remove_prefix(slash + 1);
    auto& output = help ? std::cout : std::cerr;
    output << "Usage: " << program << ' ' << arguments << '\n';
    return help ? EXIT_SUCCESS : EXIT_FAILURE;
}

} // namespace mbspec

#endif
