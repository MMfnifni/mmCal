#pragma once

#include "mathematics/angle.hpp"

#include <charconv>
#include <cstddef>
#include <optional>
#include <ostream>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::cli {

enum class InputMode {
    Interactive,
    Evaluate,
    Batch
};

enum class ExitCode : int {
    Success = 0,
    Argument = 2,
    Syntax = 3,
    Evaluation = 4,
    Internal = 5
};

struct StartupOptions final {
    std::optional<std::size_t> fixedDigits;
    std::optional<mathematics::AngleUnit> angleUnit;
    InputMode inputMode = InputMode::Interactive;
    std::optional<std::string> expression;
    bool showHelp = false;
};

[[nodiscard]] inline std::size_t parseFixedDigitsOption(std::string_view text) {
    std::size_t digits = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), digits);
    if (text.empty() || result.ec != std::errc{}
        || result.ptr != text.data() + text.size() || digits > 1000)
        throw std::invalid_argument("--fix expects an integer from 0 to 1000");
    return digits;
}

[[nodiscard]] inline StartupOptions parseStartupOptions(
    std::span<const std::string_view> arguments) {
    StartupOptions options;

    const auto selectInputMode = [&options](InputMode mode, std::string_view option) {
        if (options.inputMode != InputMode::Interactive)
            throw std::invalid_argument(
                std::string{option} + " cannot be combined with another input mode");
        options.inputMode = mode;
    };

    for (std::size_t index = 0; index < arguments.size(); ++index) {
        const std::string_view argument = arguments[index];
        if (argument == "--help" || argument == "-h") {
            options.showHelp = true;
            continue;
        }

        if (argument == "--fix") {
            if (++index >= arguments.size())
                throw std::invalid_argument("--fix requires a value");
            options.fixedDigits = parseFixedDigitsOption(arguments[index]);
            continue;
        }

        if (argument == "--angle") {
            if (++index >= arguments.size())
                throw std::invalid_argument("--angle requires deg, rad, or grad");
            const auto unit = mathematics::AngleSemantics::parseUnit(arguments[index]);
            if (!unit)
                throw std::invalid_argument("--angle expects deg, rad, or grad");
            options.angleUnit = *unit;
            continue;
        }

        if (argument == "--eval") {
            selectInputMode(InputMode::Evaluate, argument);
            if (++index >= arguments.size())
                throw std::invalid_argument("--eval requires an expression");
            options.expression = std::string{arguments[index]};
            continue;
        }

        if (argument == "--batch" || argument == "--bach") {
            selectInputMode(InputMode::Batch, argument);
            continue;
        }

        throw std::invalid_argument("Unknown command-line option: " + std::string{argument});
    }

    return options;
}

[[nodiscard]] inline StartupOptions parseStartupOptions(int argc, char* argv[]) {
    std::vector<std::string_view> arguments;
    arguments.reserve(argc > 1 ? static_cast<std::size_t>(argc - 1) : 0);
    for (int index = 1; index < argc; ++index)
        arguments.emplace_back(argv[index]);
    return parseStartupOptions(arguments);
}

inline void printUsage(std::ostream& output) {
    output
        << "Usage: mmCal [--fix <0..1000>] [--angle <deg|rad|grad>]\n"
        << "       mmCal [options] --eval <expression>\n"
        << "       mmCal [options] --batch    (alias: --bach)\n"
        << "Interactive help: :help [function]\n";
}

} // namespace mmcal::cli
