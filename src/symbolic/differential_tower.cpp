#include "differential_tower.hpp"

#include "symbolic/polynomial.hpp"

#include <utility>

namespace mmcal::symbolic::risch {

DifferentialExtension::DifferentialExtension(
    DifferentialExtensionKind kind,
    expression::Symbol generator,
    expression::Expr source,
    expression::Expr differentialCoefficient,
    std::size_t level)
    : kind_(kind),
      generator_(std::move(generator)),
      source_(std::move(source)),
      differentialCoefficient_(std::move(differentialCoefficient)),
      level_(level) {}

DifferentialExtensionKind DifferentialExtension::kind() const noexcept {
    return kind_;
}

const expression::Symbol& DifferentialExtension::generator() const noexcept {
    return generator_;
}

const expression::Expr& DifferentialExtension::source() const noexcept {
    return source_;
}

const expression::Expr& DifferentialExtension::differentialCoefficient() const noexcept {
    return differentialCoefficient_;
}

std::size_t DifferentialExtension::level() const noexcept {
    return level_;
}

DifferentialTower::DifferentialTower(
    expression::Symbol baseVariable,
    std::size_t maximumDepth)
    : baseVariable_(std::move(baseVariable)), maximumDepth_(maximumDepth) {}

const expression::Symbol& DifferentialTower::baseVariable() const noexcept {
    return baseVariable_;
}

std::size_t DifferentialTower::depth() const noexcept {
    return extensions_.size();
}

std::size_t DifferentialTower::maximumDepth() const noexcept {
    return maximumDepth_;
}

bool DifferentialTower::valid() const noexcept {
    return baseVariable_.valid();
}

std::span<const DifferentialExtension> DifferentialTower::extensions() const noexcept {
    return extensions_;
}

std::optional<std::size_t> DifferentialTower::levelOf(
    const expression::Symbol& generator) const noexcept {
    if (generator == baseVariable_)
        return 0;
    for (const DifferentialExtension& extension : extensions_)
        if (extension.generator() == generator)
            return extension.level();
    return std::nullopt;
}

DifferentialTowerAppendResult DifferentialTower::appendPrimitive(
    expression::Symbol generator,
    expression::Expr source,
    expression::Expr derivative) {
    return append(
        DifferentialExtensionKind::Primitive,
        std::move(generator), std::move(source), std::move(derivative));
}

DifferentialTowerAppendResult DifferentialTower::appendExponential(
    expression::Symbol generator,
    expression::Expr source,
    expression::Expr logarithmicDerivative) {
    return append(
        DifferentialExtensionKind::Exponential,
        std::move(generator), std::move(source), std::move(logarithmicDerivative));
}

DifferentialTowerAppendResult DifferentialTower::append(
    DifferentialExtensionKind kind,
    expression::Symbol generator,
    expression::Expr source,
    expression::Expr differentialCoefficient) {
    if (!baseVariable_.valid())
        return {DifferentialTowerError::InvalidBaseVariable, std::nullopt};
    if (!generator.valid())
        return {DifferentialTowerError::InvalidGenerator, std::nullopt};
    if (generator == baseVariable_)
        return {DifferentialTowerError::GeneratorIsBaseVariable, std::nullopt};
    if (levelOf(generator))
        return {DifferentialTowerError::DuplicateGenerator, std::nullopt};
    if (containsSymbol(source, generator)
        || containsSymbol(differentialCoefficient, generator))
        return {DifferentialTowerError::SelfDependentDefinition, std::nullopt};
    if (extensions_.size() >= maximumDepth_)
        return {DifferentialTowerError::DepthLimit, std::nullopt};

    const std::size_t level = extensions_.size() + 1;
    extensions_.emplace_back(
        kind, std::move(generator), std::move(source),
        std::move(differentialCoefficient), level);
    return {DifferentialTowerError::None, level};
}

} // namespace mmcal::symbolic::risch
