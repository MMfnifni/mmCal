// 予約定数・domain Symbolの登録の回帰テスト
#include "symbol_registry_tests.hpp"

#include "symbols/symbol_registry.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>

namespace mmcal::tests {

void runSymbolRegistryTests(TestRunner& tests) {
    using symbols::PredefinedSymbolDefinition;
    using symbols::PredefinedSymbolId;
    using symbols::PredefinedSymbolKind;
    using symbols::SymbolRegistry;
    using symbols::SymbolTable;

    SymbolTable table;
    SymbolRegistry registry = SymbolRegistry::defaults(table);
    tests.expectEqual(registry.size(), std::size_t{14},
        "SymbolRegistry: registers all predefined symbols");

    const PredefinedSymbolDefinition* pi = registry.find("Pi");
    tests.expect(pi && pi->kind == PredefinedSymbolKind::SymbolicConstant,
        "SymbolRegistry: Pi is an exact symbolic constant");
    tests.expect(pi && pi->protectedName,
        "SymbolRegistry: Pi is protected");
    tests.expect(registry.isSymbolicConstant("E"),
        "SymbolRegistry: E is an exact symbolic constant");

    const PredefinedSymbolDefinition* imaginary = registry.find("I");
    tests.expect(imaginary && imaginary->kind == PredefinedSymbolKind::ImaginaryUnit,
        "SymbolRegistry: I is the imaginary-unit literal");

    const PredefinedSymbolDefinition* realDomain = registry.find("Real");
    tests.expect(realDomain && realDomain->kind == PredefinedSymbolKind::MathematicalDomain
        && realDomain->protectedName,
        "SymbolRegistry: Real is a protected mathematical domain symbol");
    tests.expect(registry.contains("Complex") && registry.contains("Rational")
        && registry.contains("Integer"),
        "SymbolRegistry: all numeric domain symbols are predefined");

    const PredefinedSymbolDefinition* booleanTrue = registry.find("True");
    tests.expect(booleanTrue && booleanTrue->kind == PredefinedSymbolKind::BooleanTrue,
        "SymbolRegistry: True is a Boolean literal");
    tests.expect(!registry.isSymbolicConstant("True"),
        "SymbolRegistry: Boolean literals are not symbolic constants");

    const PredefinedSymbolDefinition* radians = registry.find("Rad");
    tests.expect(radians && radians->kind == PredefinedSymbolKind::EnumeratedValue
        && radians->protectedName && registry.contains("Deg") && registry.contains("Grad"),
        "SymbolRegistry: angle-mode values are protected enumerated symbols");

    const auto constants = registry.sourcePredefinedNames();
    tests.expect(constants.contains("Pi") && constants.contains("I")
        && constants.contains("True") && constants.contains("False"),
        "SymbolRegistry: exports all protected literal and constant names to the parser");
    tests.expect(!registry.contains("Tau") && !registry.contains("NA") && !registry.contains("ESP"),
        "SymbolRegistry: removed legacy constants are no longer predefined");
    tests.expect(!registry.contains("x") && !registry.isProtected("x"),
        "SymbolRegistry: ordinary variables are not predefined");

    tests.expectThrows<std::invalid_argument>([&] {
        registry.add(
            "Pi",
            PredefinedSymbolId::Pi,
            PredefinedSymbolKind::SymbolicConstant,
            true);
    }, "SymbolRegistry: rejects duplicate registrations");
}

} // namespace mmcal::tests
