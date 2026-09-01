// CLIの起動・入力ループ・表示設定
#include "approximation/approximation_context.hpp"
#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "cli/repl_help.hpp"
#include "cli/repl_layout.hpp"
#include "cli/startup_options.hpp"
#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "expression/array_utils.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "../../version.h"

#include <algorithm>
#include <charconv>
#include <csignal>
#include <cstddef>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>
#if defined(__linux__) || defined(__APPLE__)
#include <sys/ioctl.h>
#endif

#if defined(_WIN32)
#define NOMINMAX
#include <windows.h>
#elif defined(__linux__) || defined(__APPLE__)
#include <unistd.h>
#endif

namespace {

struct DisplaySettings final {
    std::optional<std::size_t> fixedDigits;
    mmcal::cli::OutputLayout layout = mmcal::cli::OutputLayout::Auto;
};

mmcal::evaluation::EvaluationCancellationToken consoleCancellation;

#if defined(_WIN32)
BOOL WINAPI consoleControlHandler(DWORD controlType) {
    if (controlType != CTRL_C_EVENT && controlType != CTRL_BREAK_EVENT)
        return FALSE;
    consoleCancellation.requestCancellation();
    return TRUE;
}

class ConsoleCancellationScope final {
public:
    ConsoleCancellationScope() {
        consoleCancellation.reset();
        installed_ = ::SetConsoleCtrlHandler(consoleControlHandler, TRUE) != FALSE;
    }

    ConsoleCancellationScope(const ConsoleCancellationScope&) = delete;
    ConsoleCancellationScope& operator=(const ConsoleCancellationScope&) = delete;

    ~ConsoleCancellationScope() {
        if (installed_)
            static_cast<void>(::SetConsoleCtrlHandler(consoleControlHandler, FALSE));
    }

private:
    bool installed_ = false;
};
#else
void consoleInterruptHandler(int) {
    consoleCancellation.requestSignalCancellation();
}

class ConsoleCancellationScope final {
public:
    ConsoleCancellationScope() {
        consoleCancellation.reset();
        previous_ = std::signal(SIGINT, consoleInterruptHandler);
    }

    ConsoleCancellationScope(const ConsoleCancellationScope&) = delete;
    ConsoleCancellationScope& operator=(const ConsoleCancellationScope&) = delete;

    ~ConsoleCancellationScope() {
        static_cast<void>(std::signal(SIGINT, previous_));
    }

private:
    using SignalHandler = void (*)(int);
    SignalHandler previous_ = SIG_DFL;
};
#endif

[[nodiscard]] mmcal::expression::Expr evaluateWithConsoleCancellation(
    mmcal::kernel::KernelSession& session,
    std::string_view source) {
    ConsoleCancellationScope cancellationScope;
    return session.evaluate(source, consoleCancellation);
}

[[nodiscard]] std::string displayModeName(const DisplaySettings& settings) {
    if (!settings.fixedDigits)
        return "Exact";
    return "Fixed(" + std::to_string(*settings.fixedDigits) + ")";
}

[[nodiscard]] std::string angleModeName(const mmcal::kernel::KernelSession& session) {
    return std::string{
        mmcal::mathematics::AngleSemantics::canonicalName(session.defaultAngleUnit())};
}

void setConsoleTitle(std::string_view title) {
#if defined(_WIN32)
    const std::string owned{title};
    static_cast<void>(::SetConsoleTitleA(owned.c_str()));
#elif defined(__linux__) || defined(__APPLE__)
    if (::isatty(STDOUT_FILENO))
        std::cout << "\033]0;" << title << '\a' << std::flush;
#else
    static_cast<void>(title);
#endif
}

void updateConsoleTitle(
    const mmcal::kernel::KernelSession& session,
    const DisplaySettings& settings) {
    setConsoleTitle(
        "mmCal " MMCAL_VERSION_STRING " - " + angleModeName(session) + " - "
        + displayModeName(settings) + " - Layout("
        + std::string{mmcal::cli::outputLayoutName(settings.layout)} + ")");
}

[[nodiscard]] std::size_t nextGuardDigits(std::size_t current) {
    const std::size_t growth = std::max<std::size_t>(8, current / 2);
    if (growth > std::numeric_limits<std::size_t>::max() - current)
        throw std::overflow_error("Approximation precision is too large");
    return current + growth;
}

[[nodiscard]] bool intervalIsExactZero(
    const mmcal::approximation::RealInterval& interval) noexcept {
    return interval.isPoint() && interval.lower().isZero();
}

[[nodiscard]] std::optional<mmcal::expression::Expr> fixedFromCertified(
    const mmcal::approximation::CertifiedValue& value,
    std::size_t digits) {
    using mmcal::expression::Expr;
    using mmcal::numeric::ComplexDecimalApproximation;
    using mmcal::numeric::DecimalApproximation;

    if (value.isReal()) {
        const auto decimal = DecimalApproximation::fromCertifiedInterval(
            value.asReal().lower().toRational(),
            value.asReal().upper().toRational(),
            digits);
        return decimal ? std::optional<Expr>{Expr{*decimal}} : std::nullopt;
    }

    const auto& complex = value.asComplex();
    const auto real = DecimalApproximation::fromCertifiedInterval(
        complex.real().lower().toRational(),
        complex.real().upper().toRational(),
        digits);
    const auto imaginary = DecimalApproximation::fromCertifiedInterval(
        complex.imaginary().lower().toRational(),
        complex.imaginary().upper().toRational(),
        digits);
    if (!real || !imaginary)
        return std::nullopt;

    return Expr{ComplexDecimalApproximation::fromComponents(
        *real,
        *imaginary,
        intervalIsExactZero(complex.real()),
        intervalIsExactZero(complex.imaginary()))};
}

[[nodiscard]] mmcal::expression::Expr fixedApproximation(
    const mmcal::expression::Expr& expression,
    std::size_t digits,
    const mmcal::kernel::KernelSession& session) {
    using mmcal::expression::Expr;
    using mmcal::numeric::ComplexDecimalApproximation;
    using mmcal::numeric::DecimalApproximation;

    if (expression.isArray()) {
        const auto& array = expression.asArray();
        std::vector<Expr> elements;
        elements.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i)
            elements.push_back(fixedApproximation(array.element(i), digits, session));
        return Expr::array(array.shape, std::move(elements));
    }
    if (expression.isList()) {
        const auto& list = expression.asList();
        std::vector<Expr> elements;
        elements.reserve(list.elements.size());
        for (const Expr& element : list.elements)
            elements.push_back(fixedApproximation(element, digits, session));
        return mmcal::expression::braceValue(std::move(elements));
    }

    if (expression.isCall()) {
        const auto& call = expression.asCall();
        const auto* definition = session.builtinRegistry().find(call.head);
        if (definition && definition->id == mmcal::evaluation::BuiltinId::UnitApplied
            && call.arguments.size() == 2 && call.arguments[1].isString()) {
            return Expr::call(call.head, {
                fixedApproximation(call.arguments[0], digits, session),
                call.arguments[1]
            });
        }
    }

    if (expression.isNumber()) {
        const auto& number = expression.asNumber();
        if (number.isReal())
            return Expr{DecimalApproximation::fromRealFixed(number.asReal(), digits)};

        const auto& complex = number.asComplex();
        return Expr{ComplexDecimalApproximation::fromComponents(
            DecimalApproximation::fromRealFixed(complex.real, digits),
            DecimalApproximation::fromRealFixed(complex.imaginary, digits),
            complex.real.isZero(),
            complex.imaginary.isZero())};
    }

    if (expression.isDecimalApproximation())
        return Expr{DecimalApproximation::fromRealFixed(
            mmcal::numeric::RealNumber{expression.asDecimalApproximation().displayedValue()},
            digits)};

    if (expression.isComplexDecimalApproximation()) {
        const auto& complex = expression.asComplexDecimalApproximation();
        return Expr{ComplexDecimalApproximation::fromComponents(
            DecimalApproximation::fromRealFixed(
                mmcal::numeric::RealNumber{complex.real().displayedValue()}, digits),
            DecimalApproximation::fromRealFixed(
                mmcal::numeric::RealNumber{complex.imaginary().displayedValue()}, digits),
            complex.realExactlyZero(),
            complex.imaginaryExactlyZero())};
    }

    const mmcal::mathematics::AngleSemantics angles{session.defaultAngleUnit()};
    const mmcal::approximation::CertifiedEvaluator certified{
        session.builtinRegistry(), session.mathRegistry(), angles};
    mmcal::approximation::ApproximationContext context{std::max<std::size_t>(digits, 1)};

    for (std::size_t attempt = 0; attempt < 16; ++attempt) {
        try {
            const auto enclosed = certified.enclose(expression, context.workingBinaryBits());
            if (!enclosed)
                return expression;
            if (const auto fixed = fixedFromCertified(*enclosed, digits))
                return *fixed;
        }
        catch (const mmcal::approximation::PrecisionInsufficient&) {
            // 表示のための精度不足なので、ガード桁だけ増やして再試行する。
        }
        catch (const mmcal::approximation::CertifiedBackendUnsupported&) {
            // :fix/--fix は表示設定であり、未対応backendを新しいエラーへ変換しない。
            return expression;
        }
        catch (const std::domain_error&) {
            // :fix/--fix は表示設定であり、未評価式へ新しいDomain/InternalErrorを導入しない。
            return expression;
        }
        context.setGuardDigits(nextGuardDigits(context.guardDigits()));
    }

    return expression;
}

[[nodiscard]] mmcal::expression::Expr displayExpression(
    const mmcal::expression::Expr& expression,
    const mmcal::kernel::KernelSession& session,
    const DisplaySettings& settings) {
    if (!settings.fixedDigits)
        return expression;
    return fixedApproximation(expression, *settings.fixedDigits, session);
}

[[nodiscard]] std::string formatForAutomatedDisplay(
    const mmcal::expression::Expr& expression,
    const mmcal::kernel::KernelSession& session,
    const DisplaySettings& settings) {
    const auto displayed = displayExpression(expression, session, settings);
    std::string text = mmcal::formatting::formatExpr(displayed);
    if (settings.fixedDigits)
        text = mmcal::formatting::trimRedundantFractionalZeros(text);
    return text;
}

[[nodiscard]] bool consoleOutputIsTerminal() noexcept {
#if defined(_WIN32)
    DWORD mode = 0;
    return ::GetConsoleMode(::GetStdHandle(STD_OUTPUT_HANDLE), &mode) != FALSE;
#elif defined(__linux__) || defined(__APPLE__)
    return ::isatty(STDOUT_FILENO) != 0;
#else
    return false;
#endif
}

[[nodiscard]] std::size_t consoleWidth() noexcept {
#if defined(_WIN32)
    CONSOLE_SCREEN_BUFFER_INFO info{};
    if (::GetConsoleScreenBufferInfo(::GetStdHandle(STD_OUTPUT_HANDLE), &info))
        return static_cast<std::size_t>(info.srWindow.Right - info.srWindow.Left + 1);
#elif defined(__linux__) || defined(__APPLE__)
    winsize size{};
    if (::isatty(STDOUT_FILENO) && ::ioctl(STDOUT_FILENO, TIOCGWINSZ, &size) == 0
        && size.ws_col != 0)
        return size.ws_col;
#endif
    return 100;
}

[[nodiscard]] std::string formatForInteractiveDisplay(
    const mmcal::expression::Expr& expression,
    const mmcal::kernel::KernelSession& session,
    const DisplaySettings& settings,
    std::size_t prefixWidth) {
    const auto displayed = displayExpression(expression, session, settings);
    const std::size_t width = consoleWidth();
    const std::size_t available = width > prefixWidth ? width - prefixWidth : width;
    const auto layout = settings.layout == mmcal::cli::OutputLayout::Auto
            && !consoleOutputIsTerminal()
        ? mmcal::cli::OutputLayout::Single
        : settings.layout;
    std::string text = mmcal::cli::formatReplExpression(displayed, layout, available);
    if (settings.fixedDigits)
        text = mmcal::formatting::trimRedundantFractionalZeros(text);
    return text;
}

void printInteractiveResult(
    std::ostream& output,
    std::string_view prefix,
    std::string_view formatted) {
    output << prefix;
    const std::string continuation(prefix.size(), ' ');
    for (const char character : formatted) {
        output.put(character);
        if (character == '\n')
            output << continuation;
    }
    output << '\n';
}

[[nodiscard]] std::string_view trim(std::string_view text) noexcept {
    while (!text.empty() && (text.front() == ' ' || text.front() == '\t'))
        text.remove_prefix(1);
    while (!text.empty() && (text.back() == ' ' || text.back() == '\t'))
        text.remove_suffix(1);
    return text;
}

[[nodiscard]] bool handleFixCommand(
    std::string_view line,
    DisplaySettings& settings,
    const mmcal::kernel::KernelSession& session,
    std::ostream& output,
    bool updateTitle) {
    line = trim(line);
    if (!line.starts_with(":fix"))
        return false;
    if (line.size() > 4 && line[4] != ' ' && line[4] != '\t')
        return false;

    std::string_view argument = trim(line.substr(4));
    if (argument.empty()) {
        output << "Display: " << displayModeName(settings) << '\n';
        return true;
    }

    if (argument == "off") {
        settings.fixedDigits.reset();
        output << "Display: Exact\n";
        if (updateTitle)
            updateConsoleTitle(session, settings);
        return true;
    }

    std::size_t digits = 0;
    const auto result = std::from_chars(
        argument.data(), argument.data() + argument.size(), digits);
    if (result.ec != std::errc{} || result.ptr != argument.data() + argument.size()
        || digits > 1000) {
        output << "Usage: :fix <0..1000>|off\n";
        return true;
    }

    settings.fixedDigits = digits;
    output << "Display: " << displayModeName(settings) << '\n';
    if (updateTitle)
        updateConsoleTitle(session, settings);
    return true;
}

[[nodiscard]] bool handleLayoutCommand(
    std::string_view line,
    DisplaySettings& settings,
    const mmcal::kernel::KernelSession& session,
    std::ostream& output) {
    line = trim(line);
    if (!line.starts_with(":layout"))
        return false;
    if (line.size() > 7 && line[7] != ' ' && line[7] != '\t')
        return false;

    const std::string_view argument = trim(line.substr(7));
    if (argument.empty()) {
        output << "Layout: " << mmcal::cli::outputLayoutName(settings.layout) << '\n';
        return true;
    }
    if (argument == "auto")
        settings.layout = mmcal::cli::OutputLayout::Auto;
    else if (argument == "single")
        settings.layout = mmcal::cli::OutputLayout::Single;
    else if (argument == "multi")
        settings.layout = mmcal::cli::OutputLayout::Multi;
    else {
        output << "Usage: :layout <auto|single|multi>\n";
        return true;
    }
    output << "Layout: " << mmcal::cli::outputLayoutName(settings.layout) << '\n';
    updateConsoleTitle(session, settings);
    return true;
}

[[nodiscard]] bool isExitCommand(std::string_view line) noexcept {
    line = trim(line);
    return line == ":quit" || line == ":exit";
}

[[nodiscard]] bool handleStatusCommand(
    std::string_view line,
    const mmcal::kernel::KernelSession& session,
    const DisplaySettings& settings,
    std::ostream& output) {
    if (trim(line) != ":status")
        return false;

    output << "Angle: " << angleModeName(session) << '\n'
           << "Display: " << displayModeName(settings) << '\n'
           << "Layout: " << mmcal::cli::outputLayoutName(settings.layout) << '\n'
           << "Evaluation: Exact-first\n"
           << "Definitions: "
           << session.environment().size() + session.userFunctions().size() << '\n'
           << "History: " << session.historySize() << '\n';
    return true;
}

[[nodiscard]] mmcal::cli::ExitCode exitCodeFor(mmcal::error::CalcErrorType type) noexcept {
    using mmcal::cli::ExitCode;
    using mmcal::error::CalcErrorType;

    switch (type) {
    case CalcErrorType::Syntax:
    case CalcErrorType::ResourceLimit:
        return ExitCode::Syntax;
    case CalcErrorType::Internal:
        return ExitCode::Internal;
    case CalcErrorType::Domain:
    case CalcErrorType::Type:
    case CalcErrorType::Overflow:
    case CalcErrorType::Name:
    case CalcErrorType::Evaluation:
        return ExitCode::Evaluation;
    }
    return ExitCode::Internal;
}

void printDiagnostics(
    const mmcal::kernel::KernelSession& session,
    std::ostream& output) {
    for (const mmcal::evaluation::EvaluationDiagnostic& diagnostic : session.diagnostics()) {
        output << (diagnostic.severity == mmcal::evaluation::DiagnosticSeverity::Warning
            ? "WARN: " : "INFO: ") << diagnostic.message;
        if (diagnostic.previousExpression)
            output << " (was "
                   << mmcal::formatting::formatExpr(*diagnostic.previousExpression) << ')';
        output << '\n';
    }
}

struct AutomatedLineResult final {
    mmcal::cli::ExitCode exitCode = mmcal::cli::ExitCode::Success;
    bool stop = false;
};

[[nodiscard]] AutomatedLineResult runAutomatedLine(
    std::string_view line,
    mmcal::kernel::KernelSession& session,
    DisplaySettings& settings,
    std::ostream& output,
    std::ostream& diagnostics) {
    using namespace mmcal;

    if (mmcal::cli::handleReplHelpCommand(line, session.builtinRegistry(), output)
        || handleFixCommand(line, settings, session, output, false)
        || handleStatusCommand(line, session, settings, output))
        return {};

    if (isExitCommand(line))
        return {cli::ExitCode::Success, true};

    const std::string_view commandLine = trim(line);
    if (!commandLine.empty() && commandLine.front() == ':') {
        diagnostics << "Unknown command\n";
        return {cli::ExitCode::Evaluation, false};
    }

    try {
        const expression::Expr result = evaluateWithConsoleCancellation(session, line);
        if (session.exitRequested())
            return {cli::ExitCode::Success, true};
        if (session.clearRequested()) {
            output << "Cleared\n";
            return {};
        }

        printDiagnostics(session, diagnostics);
        output << formatForAutomatedDisplay(result, session, settings) << '\n';
        return {};
    }
    catch (const error::CalcError& exception) {
        diagnostics << error::errorMessage(exception) << '\n';
        return {exitCodeFor(exception.type()), false};
    }
    catch (const std::exception& exception) {
        diagnostics << error::errorMessage(
            error::CalcError{error::CalcErrorType::Internal, exception.what()}) << '\n';
        return {cli::ExitCode::Internal, false};
    }
}

[[nodiscard]] int runAutomated(
    const mmcal::cli::StartupOptions& startup,
    mmcal::kernel::KernelSession& session,
    DisplaySettings& settings) {
    using mmcal::cli::ExitCode;

    ExitCode overall = ExitCode::Success;
    const auto runLine = [&](std::string_view line) {
        const AutomatedLineResult result = runAutomatedLine(
            line, session, settings, std::cout, std::cerr);
        if (static_cast<int>(result.exitCode) > static_cast<int>(overall))
            overall = result.exitCode;
        return result.stop;
    };

    if (startup.inputMode == mmcal::cli::InputMode::Evaluate) {
        static_cast<void>(runLine(*startup.expression));
        return static_cast<int>(overall);
    }

    std::string line;
    while (std::getline(std::cin, line)) {
        if (trim(line).empty())
            continue;
        if (runLine(line))
            break;
    }
    return static_cast<int>(overall);
}

} // namespace

int main(int argc, char* argv[]) {
    using namespace mmcal;

    cli::StartupOptions startup;
    try {
        startup = cli::parseStartupOptions(argc, argv);
    }
    catch (const std::exception& error) {
        std::cerr << "Argument error: " << error.what() << '\n';
        cli::printUsage(std::cerr);
        return static_cast<int>(cli::ExitCode::Argument);
    }

    if (startup.showHelp) {
        cli::printUsage(std::cout);
        return 0;
    }

    kernel::KernelSession session;
    if (startup.angleUnit)
        session.setDefaultAngleUnit(*startup.angleUnit);

    DisplaySettings displaySettings;
    displaySettings.fixedDigits = startup.fixedDigits;
    if (startup.outputLayout)
        displaySettings.layout = *startup.outputLayout;

    if (startup.inputMode != cli::InputMode::Interactive)
        return runAutomated(startup, session, displaySettings);

    updateConsoleTitle(session, displaySettings);

    std::cout << "================================\n"
                 "  mm Calculator " MMCAL_VERSION_STRING "\n"
                 "        mmKreutzef 2021-2026\n"
                 "================================\n";

    std::string line;
    while (true) {
        const std::size_t inputNumber = session.nextInputNumber();
        std::cout << "\nIn [" << inputNumber << "]> ";
        if (!std::getline(std::cin, line))
            break;
        if (trim(line).empty())
            continue;

        if (cli::handleReplHelpCommand(line, session.builtinRegistry(), std::cout)
            || handleFixCommand(line, displaySettings, session, std::cout, true)
            || handleLayoutCommand(line, displaySettings, session, std::cout)
            || handleStatusCommand(line, session, displaySettings, std::cout))
            continue;
        if (isExitCommand(line))
            break;

        const std::string_view commandLine = trim(line);
        if (!commandLine.empty() && commandLine.front() == ':') {
            std::cout << "Unknown command\n";
            continue;
        }

        try {
            const expression::Expr result = evaluateWithConsoleCancellation(session, line);
            if (session.exitRequested())
                break;
            if (session.clearRequested()) {
                std::cout << "Cleared\n";
                updateConsoleTitle(session, displaySettings);
                continue;
            }

            for (const evaluation::EvaluationDiagnostic& diagnostic : session.diagnostics()) {
                std::cout << (diagnostic.severity == evaluation::DiagnosticSeverity::Warning
                    ? "WARN: " : "INFO: ") << diagnostic.message;
                if (diagnostic.previousExpression)
                    std::cout << " (was "
                              << formatting::formatExpr(*diagnostic.previousExpression) << ')';
                std::cout << '\n';
            }
            const std::string prefix = "Out[" + std::to_string(inputNumber) + "]> ";
            printInteractiveResult(
                std::cout, prefix,
                formatForInteractiveDisplay(
                    result, session, displaySettings, prefix.size()));
            updateConsoleTitle(session, displaySettings);
        }
        catch (const error::CalcError& exception) {
            std::cout << error::errorMessage(exception) << '\n';
        }
        catch (const std::exception& exception) {
            std::cout << error::errorMessage(
                error::CalcError{error::CalcErrorType::Internal, exception.what()}) << '\n';
        }
    }

    std::cout << "\nbye..nara...\n";
    return 0;
}
