#pragma once

#include "expression/symbol.hpp"
#include "symbol_id.hpp"
#include "symbol_table.hpp"

#include <cstddef>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>

namespace mmcal::symbols {

// ソース上で予約される組み込みシンボルを識別する。
enum class PredefinedSymbolId {
    Pi,
    E,
    Phi,
    I,
    True,
    False,
    IntegerDomain,
    RationalDomain,
    RealDomain,
    ComplexDomain,
    Infinity,
    ComplexInfinity,
    Indeterminate,
    DegreeUnit,
    RadianUnit,
    GradianUnit
};

// 名前をLoweringしたときの意味を定義する。
enum class PredefinedSymbolKind {
    SymbolicConstant,
    ImaginaryUnit,
    BooleanTrue,
    BooleanFalse,
    MathematicalDomain,
    ExceptionalValue,
    EnumeratedValue
};

struct PredefinedSymbolDefinition final {
    expression::Symbol symbol;
    PredefinedSymbolId id = PredefinedSymbolId::Pi;
    PredefinedSymbolKind kind = PredefinedSymbolKind::SymbolicConstant;
    bool protectedName = true;
};

// Symbolのidentityと、そのシンボルが持つ組み込み意味・保護属性を分離して管理する。
class SymbolRegistry final {
public:
    explicit SymbolRegistry(SymbolTable& symbolTable);

    [[nodiscard]] static SymbolRegistry defaults(SymbolTable& symbolTable);

    void add(
        std::string_view name,
        PredefinedSymbolId id,
        PredefinedSymbolKind kind,
        bool protectedName = true);
    [[nodiscard]] const PredefinedSymbolDefinition* find(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] const PredefinedSymbolDefinition* find(std::string_view name) const noexcept;
    [[nodiscard]] const PredefinedSymbolDefinition* find(PredefinedSymbolId id) const noexcept;
    [[nodiscard]] bool contains(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] bool contains(std::string_view name) const noexcept;
    [[nodiscard]] bool isProtected(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] bool isProtected(std::string_view name) const noexcept;
    [[nodiscard]] bool isSymbolicConstant(const expression::Symbol& symbol) const noexcept;
    [[nodiscard]] bool isSymbolicConstant(std::string_view name) const noexcept;
    [[nodiscard]] std::unordered_set<std::string> sourcePredefinedNames() const;
    [[nodiscard]] std::size_t size() const noexcept;

private:
    SymbolTable& symbolTable_;
    std::unordered_map<SymbolId, PredefinedSymbolDefinition, SymbolIdHash> definitions_;
};

// 単体利用向けの既定表。通常のKernelSessionは専用表を所有する。
[[nodiscard]] const SymbolRegistry& defaultSymbolRegistry();

} // namespace mmcal::symbols
