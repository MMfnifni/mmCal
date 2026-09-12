// 組込み函数登録と属性の回帰テスト
#include "builtin_registry_tests.hpp"

#include "builtins/names.hpp"
#include "evaluation/builtin_registry.hpp"
#include "test_framework.hpp"
#include "symbols/symbol_table.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::tests {
namespace {

[[nodiscard]] std::filesystem::path projectRoot() {
    std::filesystem::path current = std::filesystem::current_path();
    for (int depth = 0; depth < 10; ++depth) {
        if (std::filesystem::exists(current / "docs/reference.ja.md")
            && std::filesystem::exists(current / "src"))
            return current;
        if (!current.has_parent_path() || current.parent_path() == current)
            break;
        current = current.parent_path();
    }

    const std::filesystem::path sourcePath = std::filesystem::path{__FILE__};
    current = sourcePath.is_absolute()
        ? sourcePath.parent_path()
        : std::filesystem::absolute(sourcePath).parent_path();
    for (int depth = 0; depth < 10; ++depth) {
        if (std::filesystem::exists(current / "docs/reference.ja.md"))
            return current;
        if (!current.has_parent_path() || current.parent_path() == current)
            break;
        current = current.parent_path();
    }
    throw std::runtime_error("Cannot locate mmCal project root for Reference registry test");
}

[[nodiscard]] std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
        return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

[[nodiscard]] std::vector<std::string> referenceSourceFunctions(
    const std::filesystem::path& path,
    std::string_view headingMarker) {
    std::ifstream stream{path};
    if (!stream)
        throw std::runtime_error("Cannot open Reference file: " + path.string());
    std::ostringstream buffer;
    buffer << stream.rdbuf();
    const std::string text = buffer.str();

    const std::size_t heading = text.find(headingMarker);
    if (heading == std::string::npos)
        throw std::runtime_error("Source-callable Reference heading was not found");
    const std::size_t fence = text.find("```text", heading);
    if (fence == std::string::npos)
        throw std::runtime_error("Source-callable Reference code block was not found");
    const std::size_t contentBegin = text.find('\n', fence);
    const std::size_t contentEnd = text.find("```", contentBegin);
    if (contentBegin == std::string::npos || contentEnd == std::string::npos)
        throw std::runtime_error("Source-callable Reference code block is incomplete");

    std::vector<std::string> names;
    std::stringstream list{text.substr(contentBegin + 1, contentEnd - contentBegin - 1)};
    std::string item;
    while (std::getline(list, item, ',')) {
        item = trim(std::move(item));
        if (!item.empty())
            names.push_back(std::move(item));
    }
    return names;
}

[[nodiscard]] std::vector<std::pair<std::string, std::string>> referenceAliases(
    const std::filesystem::path& path,
    std::string_view headingMarker) {
    std::ifstream stream{path};
    if (!stream)
        throw std::runtime_error("Cannot open Reference file: " + path.string());

    std::vector<std::pair<std::string, std::string>> aliases;
    std::string line;
    bool inSection = false;
    while (std::getline(stream, line)) {
        if (!inSection) {
            if (line.find(headingMarker) != std::string::npos)
                inSection = true;
            continue;
        }
        if (line.rfind("## ", 0) == 0)
            break;
        if (line.empty() || line.front() != '|')
            continue;

        const std::size_t firstTick = line.find('`');
        if (firstTick == std::string::npos)
            continue;
        const std::size_t secondTick = line.find('`', firstTick + 1);
        const std::size_t thirdTick = line.find('`', secondTick + 1);
        const std::size_t fourthTick = line.find('`', thirdTick + 1);
        if (secondTick == std::string::npos || thirdTick == std::string::npos
            || fourthTick == std::string::npos)
            continue;
        aliases.emplace_back(
            line.substr(firstTick + 1, secondTick - firstTick - 1),
            line.substr(thirdTick + 1, fourthTick - thirdTick - 1));
    }
    return aliases;
}

[[nodiscard]] std::string setDifferenceReport(
    const std::unordered_set<std::string>& registry,
    const std::vector<std::string>& reference) {
    const std::set<std::string> documented(reference.begin(), reference.end());
    std::vector<std::string> registryOnly;
    std::vector<std::string> referenceOnly;
    for (const std::string& name : registry)
        if (!documented.contains(name))
            registryOnly.push_back(name);
    for (const std::string& name : documented)
        if (!registry.contains(name))
            referenceOnly.push_back(name);
    std::sort(registryOnly.begin(), registryOnly.end());
    std::sort(referenceOnly.begin(), referenceOnly.end());

    std::string report;
    const auto append = [&](std::string_view label, const std::vector<std::string>& names) {
        if (names.empty())
            return;
        if (!report.empty())
            report += "; ";
        report += label;
        report += ": ";
        for (std::size_t i = 0; i < names.size(); ++i) {
            if (i != 0)
                report += ", ";
            report += names[i];
        }
    };
    append("registry-only", registryOnly);
    append("reference-only", referenceOnly);
    return report;
}

} // namespace

void runBuiltinRegistryTests(TestRunner& tests) {
    using evaluation::ArgumentEvaluation;
    using evaluation::BuiltinDefinition;
    using evaluation::BuiltinId;
    using evaluation::BuiltinRegistry;
    using symbols::SymbolTable;

    SymbolTable table;
    BuiltinRegistry registry = BuiltinRegistry::defaults(table);
    tests.expectEqual(registry.size(), std::size_t{285},
        "BuiltinRegistry: registers all current builtins");
    tests.expect(registry.contains(builtins::names::sqrt),
        "BuiltinRegistry: contains sqrt");
    tests.expect(registry.contains(builtins::names::abs)
        && registry.contains(builtins::names::sign)
        && registry.contains(builtins::names::re)
        && registry.contains(builtins::names::im)
        && registry.contains(builtins::names::conj),
        "BuiltinRegistry: contains real/complex elementary functions");
    tests.expect(!registry.contains("missing"),
        "BuiltinRegistry: rejects missing name");

    const BuiltinDefinition* sqrt = registry.find(builtins::names::sqrt);
    tests.expect(sqrt && sqrt->sourceCallable,
        "BuiltinRegistry: marks sqrt as source-callable");
    tests.expect(sqrt && sqrt->acceptsArity(1) && !sqrt->acceptsArity(2),
        "BuiltinRegistry: validates fixed arity");

    const BuiltinDefinition* numerical = registry.find(builtins::names::numericalApproximation);
    tests.expect(numerical && numerical->sourceCallable,
        "BuiltinRegistry: marks N as source-callable");
    tests.expect(numerical && numerical->acceptsArity(1)
        && numerical->acceptsArity(2) && !numerical->acceptsArity(3),
        "BuiltinRegistry: validates N arity range");

    const BuiltinDefinition* derivative = registry.find(builtins::names::derivative);
    tests.expect(derivative && derivative->sourceCallable
        && derivative->argumentEvaluation == ArgumentEvaluation::HoldAll
        && derivative->acceptsArity(2) && derivative->acceptsArity(3),
        "BuiltinRegistry: D holds the expression and accepts sequential derivative specifications");

    const BuiltinDefinition* numericIntegral = registry.find(builtins::names::numericIntegral);
    tests.expect(numericIntegral && numericIntegral->sourceCallable
        && numericIntegral->argumentEvaluation == ArgumentEvaluation::HoldFirstAndIteratorSpec
        && numericIntegral->acceptsArity(2) && numericIntegral->acceptsArity(3)
        && !numericIntegral->acceptsArity(4),
        "BuiltinRegistry: nintegrate uses iterator-spec binding semantics");
    const BuiltinDefinition* symbolicIntegral = registry.find(builtins::names::symbolicIntegral);
    tests.expect(symbolicIntegral && symbolicIntegral->sourceCallable
        && symbolicIntegral->argumentEvaluation == ArgumentEvaluation::HoldFirstAndIteratorSpec
        && symbolicIntegral->acceptsArity(2) && symbolicIntegral->acceptsArity(3),
        "BuiltinRegistry: integrate uses symbolic binder semantics with optional assumptions");
    const BuiltinDefinition* limit = registry.find(builtins::names::limit);
    tests.expect(limit && limit->sourceCallable
        && limit->argumentEvaluation == ArgumentEvaluation::HoldFirstAndIteratorSpec
        && limit->acceptsArity(2) && limit->acceptsArity(3) && limit->acceptsArity(4),
        "BuiltinRegistry: limit holds the expression and variable and accepts one-sided direction");

    const BuiltinDefinition* simplify = registry.find(builtins::names::simplify);
    tests.expect(simplify && simplify->acceptsArity(1) && simplify->acceptsArity(2)
        && !simplify->acceptsArity(3),
        "BuiltinRegistry: simplify accepts optional assumptions");

    const BuiltinDefinition* numericalApproximation = registry.find(
        builtins::names::numericalApproximation);
    tests.expect(numericalApproximation
        && numericalApproximation->argumentEvaluation == ArgumentEvaluation::HoldFirst,
        "BuiltinRegistry: N holds its expression so precision-aware builtins can see the request");

    const BuiltinDefinition* solve = registry.find(builtins::names::solve);
    tests.expect(solve && solve->argumentEvaluation == ArgumentEvaluation::HoldAll
        && solve->acceptsArity(2) && solve->acceptsArity(3)
        && !solve->acceptsArity(4),
        "BuiltinRegistry: solve accepts an optional domain/constraint specification");

    const BuiltinDefinition* set = registry.find(builtins::names::set);
    tests.expect(set && set->argumentEvaluation == ArgumentEvaluation::HoldFirst,
        "BuiltinRegistry: Set holds its first argument");
    const BuiltinDefinition* rule = registry.find(builtins::names::rule);
    tests.expect(rule && rule->sourceCallable && rule->acceptsArity(2)
        && rule->argumentEvaluation == ArgumentEvaluation::HoldFirst,
        "BuiltinRegistry: Rule is source-callable and holds only its left argument");

    const BuiltinDefinition* ifDefinition = registry.find(builtins::names::ifThenElse);
    tests.expect(ifDefinition && ifDefinition->sourceCallable
        && ifDefinition->argumentEvaluation == ArgumentEvaluation::HoldAll,
        "BuiltinRegistry: If is a held source-callable special form");

    const BuiltinDefinition* map = registry.find(builtins::names::map);
    tests.expect(map && map->sourceCallable
        && map->argumentEvaluation == ArgumentEvaluation::HoldFirst
        && map->acceptsArity(2),
        "BuiltinRegistry: map holds its function symbol while evaluating the container");
    const BuiltinDefinition* range = registry.find(builtins::names::range);
    tests.expect(range && range->sourceCallable
        && range->argumentEvaluation == ArgumentEvaluation::All
        && range->acceptsArity(1) && range->acceptsArity(2) && range->acceptsArity(3),
        "BuiltinRegistry: range accepts one to three evaluated exact bounds");
    const BuiltinDefinition* tableBuiltin = registry.find(builtins::names::table);
    tests.expect(tableBuiltin && tableBuiltin->sourceCallable
        && tableBuiltin->argumentEvaluation == ArgumentEvaluation::HoldFirstAndTableIteratorSpec
        && tableBuiltin->acceptsArity(2),
        "BuiltinRegistry: table holds the body and binds its iterator specification");

    const auto sourceFunctions = registry.sourceFunctionNames();
    tests.expectEqual(sourceFunctions.size(), std::size_t{265},
        "BuiltinRegistry: exports the documented source-callable name count");
    tests.expect(sourceFunctions.contains("explain")
        && sourceFunctions.contains("sqrt") && sourceFunctions.contains("sin")
        && sourceFunctions.contains("cos") && sourceFunctions.contains("tan")
        && sourceFunctions.contains("cot") && sourceFunctions.contains("sec")
        && sourceFunctions.contains("csc") && sourceFunctions.contains("asin")
        && sourceFunctions.contains("acos") && sourceFunctions.contains("atan")
        && sourceFunctions.contains("atan2") && sourceFunctions.contains("sinh")
        && sourceFunctions.contains("cosh") && sourceFunctions.contains("tanh")
        && sourceFunctions.contains("asinh") && sourceFunctions.contains("acosh")
        && sourceFunctions.contains("atanh") && sourceFunctions.contains("csch")
        && sourceFunctions.contains("sech") && sourceFunctions.contains("coth")
        && sourceFunctions.contains("arg")
        && sourceFunctions.contains("abs") && sourceFunctions.contains("sign")
        && sourceFunctions.contains("re") && sourceFunctions.contains("im")
        && sourceFunctions.contains("conj") && sourceFunctions.contains("element")
        && sourceFunctions.contains("log") && sourceFunctions.contains("exp")
        && sourceFunctions.contains("N") && sourceFunctions.contains("Rule")
        && sourceFunctions.contains("if")
        && sourceFunctions.contains("In") && sourceFunctions.contains("Out")
        && sourceFunctions.contains("simplify") && sourceFunctions.contains("fullSimplify")
        && sourceFunctions.contains("expand") && sourceFunctions.contains("factor")
        && sourceFunctions.contains("collect") && sourceFunctions.contains("solve")
        && sourceFunctions.contains("D") && sourceFunctions.contains("floor")
        && sourceFunctions.contains("ceil") && sourceFunctions.contains("trunc")
        && sourceFunctions.contains("round") && sourceFunctions.contains("frac")
        && sourceFunctions.contains("gcd") && sourceFunctions.contains("lcm")
        && sourceFunctions.contains("mod") && sourceFunctions.contains("rem")
        && sourceFunctions.contains("quotient")
        && sourceFunctions.contains("perm") && sourceFunctions.contains("comb")
        && sourceFunctions.contains("fib") && sourceFunctions.contains("dft")
        && sourceFunctions.contains("fft") && sourceFunctions.contains("ifft")
        && sourceFunctions.contains("convolve")
        && sourceFunctions.contains("log2") && sourceFunctions.contains("log10")
        && sourceFunctions.contains("gamma") && sourceFunctions.contains("lgamma")
        && sourceFunctions.contains("lambertw")
        && sourceFunctions.contains("erf") && sourceFunctions.contains("erfc")
        && sourceFunctions.contains("fresnelc") && sourceFunctions.contains("fresnels")
        && sourceFunctions.contains("beta") && sourceFunctions.contains("betaln")
        && sourceFunctions.contains("binom") && sourceFunctions.contains("fallingfact")
        && sourceFunctions.contains("risingfact")
        && sourceFunctions.contains("randSeed") && sourceFunctions.contains("rand")
        && sourceFunctions.contains("randint") && sourceFunctions.contains("choice")
        && sourceFunctions.contains("randn")
        && sourceFunctions.contains("transpose") && sourceFunctions.contains("madd")
        && sourceFunctions.contains("matmul") && sourceFunctions.contains("det")
        && sourceFunctions.contains("inverse") && sourceFunctions.contains("rref")
        && sourceFunctions.contains("rank") && sourceFunctions.contains("solveLinear")
        && sourceFunctions.contains("diff")
        && sourceFunctions.contains("nintegrate")
        && sourceFunctions.contains("integrate") && sourceFunctions.contains("limit")
        && sourceFunctions.contains("series") && sourceFunctions.contains("normal")
        && sourceFunctions.contains("map") && sourceFunctions.contains("range")
        && sourceFunctions.contains("table")
        && sourceFunctions.contains("root")
        && sourceFunctions.contains("cbrt") && sourceFunctions.contains("hypot")
        && sourceFunctions.contains("cis") && sourceFunctions.contains("polar")
        && sourceFunctions.contains("nextpow2")
        && sourceFunctions.contains("DtoR") && sourceFunctions.contains("DtoG")
        && sourceFunctions.contains("RtoD") && sourceFunctions.contains("RtoG")
        && sourceFunctions.contains("GtoD") && sourceFunctions.contains("GtoR")
        && sourceFunctions.contains("pow") && sourceFunctions.contains("fact")
        && sourceFunctions.contains("fract") && sourceFunctions.contains("ln")
        && sourceFunctions.contains("real") && sourceFunctions.contains("imag")
        && sourceFunctions.contains("mag") && sourceFunctions.contains("unit")
        && sourceFunctions.contains("csgn") && sourceFunctions.contains("rect")
        && sourceFunctions.contains("mmul") && sourceFunctions.contains("mtranspose")
        && sourceFunctions.contains("mdet") && sourceFunctions.contains("minverse")
        && sourceFunctions.contains("mrank")
        && sourceFunctions.contains("angleMode")
        && sourceFunctions.contains("Plot") && sourceFunctions.contains("plot")
        && sourceFunctions.contains("ParametricPlot")
        && sourceFunctions.contains("Show")
        && sourceFunctions.contains("Export")
        && sourceFunctions.contains("hypergeometric2F1")
        && sourceFunctions.contains("ellipticF")
        && sourceFunctions.contains("ellipticE")
        && sourceFunctions.contains("ellipticPi")
        && sourceFunctions.contains("Ei") && sourceFunctions.contains("Si")
        && sourceFunctions.contains("Ci") && sourceFunctions.contains("li")
        && sourceFunctions.contains("polylog")
        && sourceFunctions.contains("rejection")
        && sourceFunctions.contains("orthogonalQ")
        && sourceFunctions.contains("orthonormalQ")
        && sourceFunctions.contains("linearIndependentQ")
        && sourceFunctions.contains("gramSchmidt")
        && !sourceFunctions.contains("vdot") && !sourceFunctions.contains("vnorm")
        && !sourceFunctions.contains("vlength") && !sourceFunctions.contains("vnormalize")
        && !sourceFunctions.contains("vunit") && !sourceFunctions.contains("vdistance")
        && !sourceFunctions.contains("gradient")
        && !sourceFunctions.contains("Sin") && !sourceFunctions.contains("ArcTan")
        && !sourceFunctions.contains("Integrate") && !sourceFunctions.contains("Limit")
        && !sourceFunctions.contains("Solve") && !sourceFunctions.contains("Rationalize"),
        "BuiltinRegistry: exports canonical source-callable names with the documented Plot compatibility alias");

    const BuiltinDefinition* ln = registry.find("ln");
    const BuiltinDefinition* real = registry.find("real");
    const BuiltinDefinition* mmul = registry.find("mmul");
    tests.expect(ln && ln->id == BuiltinId::Log
        && real && real->id == BuiltinId::Re
        && mmul && mmul->id == BuiltinId::MatrixMultiply,
        "BuiltinRegistry: aliases resolve to the canonical builtin IDs");
    tests.expectEqual(registry.symbol(BuiltinId::Log).name(), std::string{"log"},
        "BuiltinRegistry: alias registration preserves the canonical symbol");
    tests.expect(!sourceFunctions.contains("Add"),
        "BuiltinRegistry: hides internal heads from source calls");


    const std::filesystem::path root = projectRoot();
    const std::vector<std::string> japaneseReference = referenceSourceFunctions(
        root / "docs/reference.ja.md", "現在のsource-callable函数一覧");
    const std::vector<std::string> englishReference = referenceSourceFunctions(
        root / "docs/reference.md", "Current source-callable function list");
    tests.expectEqual(
        std::set<std::string>(japaneseReference.begin(), japaneseReference.end()).size(),
        japaneseReference.size(),
        "BuiltinRegistry: Japanese Reference source-callable list has no duplicates");
    tests.expectEqual(
        std::set<std::string>(englishReference.begin(), englishReference.end()).size(),
        englishReference.size(),
        "BuiltinRegistry: English Reference source-callable list has no duplicates");
    tests.expectEqual(setDifferenceReport(sourceFunctions, japaneseReference), std::string{},
        "BuiltinRegistry: Japanese Reference source-callable list matches registry");
    tests.expectEqual(setDifferenceReport(sourceFunctions, englishReference), std::string{},
        "BuiltinRegistry: English Reference source-callable list matches registry");
    tests.expect(
        std::set<std::string>(japaneseReference.begin(), japaneseReference.end())
            == std::set<std::string>(englishReference.begin(), englishReference.end()),
        "BuiltinRegistry: Japanese and English Reference source-callable lists match");
    tests.expect(japaneseReference == englishReference,
        "BuiltinRegistry: Japanese and English Reference source-callable lists use the same order");
    tests.expect(std::is_sorted(japaneseReference.begin(), japaneseReference.end()),
        "BuiltinRegistry: Reference source-callable list is kept in lexical order");
    const auto japaneseAliases = referenceAliases(root / "docs/reference.ja.md", "主要alias");
    const auto englishAliases = referenceAliases(root / "docs/reference.md", "Major aliases");
    std::vector<std::pair<std::string, std::string>> registryAliases;
    for (const std::string& name : sourceFunctions) {
        const BuiltinDefinition* definition = registry.find(name);
        if (!definition)
            continue;
        const std::string canonical = registry.symbol(definition->id).name();
        if (canonical != name)
            registryAliases.emplace_back(name, canonical);
    }
    std::sort(registryAliases.begin(), registryAliases.end());
    tests.expect(japaneseAliases == englishAliases,
        "BuiltinRegistry: Japanese and English Reference alias tables match");
    tests.expect(std::is_sorted(japaneseAliases.begin(), japaneseAliases.end()),
        "BuiltinRegistry: Reference alias table is kept in lexical order");
    tests.expect(japaneseAliases == registryAliases,
        "BuiltinRegistry: Reference alias table matches source-callable registry aliases");

    tests.expectThrows<std::invalid_argument>([&] {
        registry.addAlias("ln", BuiltinId::Log);
    }, "BuiltinRegistry: rejects duplicate alias registrations");

    tests.expectThrows<std::invalid_argument>([&] {
        BuiltinRegistry isolated{table};
        isolated.addAlias("badAlias", BuiltinId::Log);
    }, "BuiltinRegistry: aliases require an existing canonical builtin ID");

    tests.expectThrows<std::invalid_argument>([&] {
        registry.add(
            builtins::names::sqrt,
            BuiltinId::Sqrt,
            1,
            1,
            ArgumentEvaluation::All,
            true);
    }, "BuiltinRegistry: rejects duplicate registrations");
}

} // namespace mmcal::tests
