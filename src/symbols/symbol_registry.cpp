// 予約定数・domain Symbolの登録
#include "symbol_registry.hpp"

#include <stdexcept>
#include <utility>

namespace mmcal::symbols {

SymbolRegistry::SymbolRegistry(SymbolTable& symbolTable)
    : symbolTable_(symbolTable) {}

SymbolRegistry SymbolRegistry::defaults(SymbolTable& symbolTable) {
    SymbolRegistry registry{symbolTable};

    registry.add("Pi", PredefinedSymbolId::Pi, PredefinedSymbolKind::SymbolicConstant);
    registry.add("E", PredefinedSymbolId::E, PredefinedSymbolKind::SymbolicConstant);
    registry.add("Phi", PredefinedSymbolId::Phi, PredefinedSymbolKind::SymbolicConstant);
    registry.add("I", PredefinedSymbolId::I, PredefinedSymbolKind::ImaginaryUnit);
    registry.add("True", PredefinedSymbolId::True, PredefinedSymbolKind::BooleanTrue);
    registry.add("False", PredefinedSymbolId::False, PredefinedSymbolKind::BooleanFalse);

    // solve / assuming / refine が共有する数体系domain。数学定数ではなく、
    // 「探索する集合」を表す保護Symbolとして登録する。
    registry.add("Integer", PredefinedSymbolId::IntegerDomain, PredefinedSymbolKind::MathematicalDomain);
    registry.add("Rational", PredefinedSymbolId::RationalDomain, PredefinedSymbolKind::MathematicalDomain);
    registry.add("Real", PredefinedSymbolId::RealDomain, PredefinedSymbolKind::MathematicalDomain);
    registry.add("Complex", PredefinedSymbolId::ComplexDomain, PredefinedSymbolKind::MathematicalDomain);
    // ではprecision/accuracyのexact値を表すsentinelとして導入する。
    // 拡張実数算術そのものはまだ自動簡約しない。
    registry.add("Infinity", PredefinedSymbolId::Infinity, PredefinedSymbolKind::SymbolicConstant);

    // session設定と角度単位指定で共有する列挙値。数学定数としては扱わない。
    registry.add("Deg", PredefinedSymbolId::DegreeUnit, PredefinedSymbolKind::EnumeratedValue);
    registry.add("Rad", PredefinedSymbolId::RadianUnit, PredefinedSymbolKind::EnumeratedValue);
    registry.add("Grad", PredefinedSymbolId::GradianUnit, PredefinedSymbolKind::EnumeratedValue);

    return registry;
}

void SymbolRegistry::add(
    std::string_view name,
    PredefinedSymbolId id,
    PredefinedSymbolKind kind,
    bool protectedName) {
    const expression::Symbol symbol = symbolTable_.intern(name);
    const auto [iterator, inserted] = definitions_.emplace(
        symbol.id(),
        PredefinedSymbolDefinition{symbol, id, kind, protectedName});
    static_cast<void>(iterator);
    if (!inserted)
        throw std::invalid_argument("Predefined symbol is already registered: " + symbol.name());
}

const PredefinedSymbolDefinition* SymbolRegistry::find(
    const expression::Symbol& symbol) const noexcept {
    const auto iterator = definitions_.find(symbol.id());
    return iterator == definitions_.end() ? nullptr : &iterator->second;
}

const PredefinedSymbolDefinition* SymbolRegistry::find(std::string_view name) const noexcept {
    const expression::Symbol symbol = symbolTable_.find(name);
    return symbol.valid() ? find(symbol) : nullptr;
}

bool SymbolRegistry::contains(const expression::Symbol& symbol) const noexcept {
    return find(symbol) != nullptr;
}

bool SymbolRegistry::contains(std::string_view name) const noexcept {
    return find(name) != nullptr;
}

bool SymbolRegistry::isProtected(const expression::Symbol& symbol) const noexcept {
    const PredefinedSymbolDefinition* definition = find(symbol);
    return definition && definition->protectedName;
}

bool SymbolRegistry::isProtected(std::string_view name) const noexcept {
    const PredefinedSymbolDefinition* definition = find(name);
    return definition && definition->protectedName;
}

bool SymbolRegistry::isSymbolicConstant(const expression::Symbol& symbol) const noexcept {
    const PredefinedSymbolDefinition* definition = find(symbol);
    return definition && definition->kind == PredefinedSymbolKind::SymbolicConstant;
}

bool SymbolRegistry::isSymbolicConstant(std::string_view name) const noexcept {
    const PredefinedSymbolDefinition* definition = find(name);
    return definition && definition->kind == PredefinedSymbolKind::SymbolicConstant;
}

std::unordered_set<std::string> SymbolRegistry::sourcePredefinedNames() const {
    std::unordered_set<std::string> result;
    result.reserve(definitions_.size());
    for (const auto& [id, definition] : definitions_) {
        static_cast<void>(id);
        result.insert(definition.symbol.name());
    }
    return result;
}

std::size_t SymbolRegistry::size() const noexcept {
    return definitions_.size();
}

const SymbolRegistry& defaultSymbolRegistry() {
    static const SymbolRegistry registry = SymbolRegistry::defaults(defaultSymbolTable());
    return registry;
}

} // namespace mmcal::symbols
