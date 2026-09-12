#pragma once

#include "expression/expr.hpp"
#include "expression/symbol.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <vector>

namespace mmcal::symbolic::risch {

// Risch tower の超越拡大種別。代数拡大は residue 表現側で独立に扱う。
enum class DifferentialExtensionKind {
    Primitive,
    Exponential
};

// Expr は認識・結果構築の境界表現であり、係数体演算そのものには使わない。
// Primitive では differentialCoefficient = D(generator)、Exponential では
// differentialCoefficient = D(generator)/generator を表す。
class DifferentialExtension final {
public:
    DifferentialExtension(
        DifferentialExtensionKind kind,
        expression::Symbol generator,
        expression::Expr source,
        expression::Expr differentialCoefficient,
        std::size_t level);

    [[nodiscard]] DifferentialExtensionKind kind() const noexcept;
    [[nodiscard]] const expression::Symbol& generator() const noexcept;
    [[nodiscard]] const expression::Expr& source() const noexcept;
    [[nodiscard]] const expression::Expr& differentialCoefficient() const noexcept;
    [[nodiscard]] std::size_t level() const noexcept;

private:
    DifferentialExtensionKind kind_;
    expression::Symbol generator_;
    expression::Expr source_;
    expression::Expr differentialCoefficient_;
    std::size_t level_ = 0;
};

enum class DifferentialTowerError {
    None,
    InvalidBaseVariable,
    InvalidGenerator,
    GeneratorIsBaseVariable,
    DuplicateGenerator,
    SelfDependentDefinition,
    DepthLimit
};

struct DifferentialTowerAppendResult final {
    DifferentialTowerError error = DifferentialTowerError::None;
    std::optional<std::size_t> level;

    [[nodiscard]] explicit operator bool() const noexcept {
        return error == DifferentialTowerError::None && level.has_value();
    }
};

class DifferentialTower final {
public:
    explicit DifferentialTower(
        expression::Symbol baseVariable,
        std::size_t maximumDepth = 16);

    [[nodiscard]] const expression::Symbol& baseVariable() const noexcept;
    [[nodiscard]] std::size_t depth() const noexcept;
    [[nodiscard]] std::size_t maximumDepth() const noexcept;
    [[nodiscard]] bool valid() const noexcept;
    [[nodiscard]] std::span<const DifferentialExtension> extensions() const noexcept;
    [[nodiscard]] std::optional<std::size_t> levelOf(
        const expression::Symbol& generator) const noexcept;

    [[nodiscard]] DifferentialTowerAppendResult appendPrimitive(
        expression::Symbol generator,
        expression::Expr source,
        expression::Expr derivative);
    [[nodiscard]] DifferentialTowerAppendResult appendExponential(
        expression::Symbol generator,
        expression::Expr source,
        expression::Expr logarithmicDerivative);

private:
    [[nodiscard]] DifferentialTowerAppendResult append(
        DifferentialExtensionKind kind,
        expression::Symbol generator,
        expression::Expr source,
        expression::Expr differentialCoefficient);

    expression::Symbol baseVariable_;
    std::size_t maximumDepth_ = 16;
    std::vector<DifferentialExtension> extensions_;
};

} // namespace mmcal::symbolic::risch
