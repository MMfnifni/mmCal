#include "risch_tower_recognizer.hpp"

#include "expression/exact_value.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "symbolic/differentiation.hpp"

#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mmcal::symbolic::risch {
namespace {

using evaluation::BuiltinId;
using expression::Expr;

[[nodiscard]] bool isArithmeticCall(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    return builtins.isCallTo(expression, BuiltinId::Add)
        || builtins.isCallTo(expression, BuiltinId::Subtract)
        || builtins.isCallTo(expression, BuiltinId::Multiply)
        || builtins.isCallTo(expression, BuiltinId::Divide)
        || builtins.isCallTo(expression, BuiltinId::Negate);
}

[[nodiscard]] bool containsAnyDependency(
    const Expr& expression,
    const DifferentialTower& tower) {
    if (expression.isSymbol()) {
        if (expression.asSymbol() == tower.baseVariable())
            return true;
        return tower.levelOf(expression.asSymbol()).has_value();
    }
    if (expression.isCall()) {
        for (const Expr& argument : expression.asCall().arguments)
            if (containsAnyDependency(argument, tower))
                return true;
    }
    else if (expression.isList()) {
        for (const Expr& element : expression.asList().elements)
            if (containsAnyDependency(element, tower))
                return true;
    }
    else if (expression.isArray()
        && expression.asArray().storageKind()
            == expression::ArrayStorageKind::Generic) {
        for (const Expr& element : expression.asArray().storedExpressions())
            if (containsAnyDependency(element, tower))
                return true;
    }
    return false;
}

[[nodiscard]] bool rationalMembership(
    const Expr& expression,
    const DifferentialTower& tower,
    const evaluation::BuiltinRegistry& builtins) {
    if (expression.isNumber())
        return true;
    if (expression.isSymbol())
        return expression.asSymbol() == tower.baseVariable()
            || tower.levelOf(expression.asSymbol()).has_value()
            || !containsAnyDependency(expression, tower);
    if (!expression.isCall())
        return false;

    // base/tower generatorを含まない評価済みscalarはconstant fieldの元として扱う。
    if (!containsAnyDependency(expression, tower))
        return true;
    const auto& arguments = expression.asCall().arguments;
    if (isArithmeticCall(expression, builtins)) {
        for (const Expr& argument : arguments)
            if (!rationalMembership(argument, tower, builtins))
                return false;
        return true;
    }
    if (builtins.isCallTo(expression, BuiltinId::Power)
        && arguments.size() == 2) {
        const auto exponent = expression::exact::realRational(arguments[1]);
        return exponent && exponent->isInteger()
            && rationalMembership(arguments[0], tower, builtins);
    }
    return false;
}

struct RecognizedSource final {
    Expr original;
    expression::Symbol generator;
};

class TowerRecognizer final {
public:
    TowerRecognizer(
        const expression::Symbol& baseVariable,
        symbols::SymbolTable& symbols,
        const evaluation::BuiltinRegistry& builtins,
        const mathematics::MathRegistry& mathematics,
        const mathematics::AngleSemantics& angles,
        std::size_t maximumDepth,
        const Expr& input)
        : symbols_(symbols),
          builtins_(builtins),
          mathematics_(mathematics),
          angles_(angles),
          tower_(baseVariable, maximumDepth) {
        collectSymbolNames(input);
    }

    [[nodiscard]] Expr rewrite(const Expr& expression) {
        for (const RecognizedSource& known : recognized_)
            if (known.original == expression)
                return Expr{known.generator};
        if (!expression.isCall())
            return expression;

        std::vector<Expr> arguments;
        arguments.reserve(expression.asCall().arguments.size());
        for (const Expr& argument : expression.asCall().arguments)
            arguments.push_back(rewrite(argument));
        Expr rebuilt = Expr::rebuildCall(expression.asCall(), std::move(arguments));

        const bool primitive = builtins_.isCallTo(expression, BuiltinId::Log);
        const bool exponential = builtins_.isCallTo(expression, BuiltinId::Exp);
        if ((!primitive && !exponential)
            || expression.asCall().arguments.size() != 1
            || !containsAnyDependency(rebuilt, tower_))
            return rebuilt;

        Expr differential = primitive
            ? differentiateExpression(
                expression, tower_.baseVariable(), builtins_, mathematics_, angles_)
            : differentiateExpression(
                expression.asCall().arguments.front(), tower_.baseVariable(),
                builtins_, mathematics_, angles_);
        differential = rewriteKnown(differential);
        differential = simplification::Simplifier{}.simplify(
            differential,
            simplification::SimplificationContext{
                builtins_, mathematics_, angles_});
        if (!rationalMembership(differential, tower_, builtins_))
            return rebuilt;

        expression::Symbol generator = nextGenerator();
        const DifferentialTowerAppendResult appended = primitive
            ? tower_.appendPrimitive(generator, rebuilt, differential)
            : tower_.appendExponential(generator, rebuilt, differential);
        if (!appended) {
            if (appended.error == DifferentialTowerError::DepthLimit)
                depthLimitReached_ = true;
            return rebuilt;
        }
        recognized_.push_back({expression, generator});
        return Expr{generator};
    }

    [[nodiscard]] DifferentialTower takeTower() {
        return std::move(tower_);
    }

    [[nodiscard]] bool depthLimitReached() const noexcept {
        return depthLimitReached_;
    }

private:
    [[nodiscard]] Expr rewriteKnown(const Expr& expression) const {
        for (const RecognizedSource& known : recognized_)
            if (known.original == expression)
                return Expr{known.generator};
        if (!expression.isCall())
            return expression;
        std::vector<Expr> arguments;
        arguments.reserve(expression.asCall().arguments.size());
        for (const Expr& argument : expression.asCall().arguments)
            arguments.push_back(rewriteKnown(argument));
        return Expr::rebuildCall(expression.asCall(), std::move(arguments));
    }

    void collectSymbolNames(const Expr& expression) {
        if (expression.isSymbol()) {
            usedNames_.insert(std::string{expression.asSymbol().view()});
            return;
        }
        if (expression.isCall()) {
            for (const Expr& argument : expression.asCall().arguments)
                collectSymbolNames(argument);
        }
        else if (expression.isList()) {
            for (const Expr& element : expression.asList().elements)
                collectSymbolNames(element);
        }
        else if (expression.isArray()
            && expression.asArray().storageKind()
                == expression::ArrayStorageKind::Generic) {
            for (const Expr& element : expression.asArray().storedExpressions())
                collectSymbolNames(element);
        }
    }

    [[nodiscard]] expression::Symbol nextGenerator() {
        for (;;) {
            const std::string name = "__risch_t" + std::to_string(nextGeneratorIndex_++);
            if (usedNames_.contains(name) || symbols_.contains(name))
                continue;
            usedNames_.insert(name);
            return symbols_.intern(name);
        }
    }

    symbols::SymbolTable& symbols_;
    const evaluation::BuiltinRegistry& builtins_;
    const mathematics::MathRegistry& mathematics_;
    const mathematics::AngleSemantics& angles_;
    DifferentialTower tower_;
    std::vector<RecognizedSource> recognized_;
    std::unordered_set<std::string> usedNames_;
    std::size_t nextGeneratorIndex_ = 1;
    bool depthLimitReached_ = false;
};

} // namespace

DifferentialTowerRecognitionResult::DifferentialTowerRecognitionResult(
    DifferentialTower recognizedTower,
    expression::Expr rewritten,
    DifferentialTowerRecognitionStatus recognitionStatus)
    : tower(std::move(recognizedTower)),
      rewrittenExpression(std::move(rewritten)),
      status(recognitionStatus) {}

bool DifferentialTowerRecognitionResult::complete() const noexcept {
    return status == DifferentialTowerRecognitionStatus::Complete;
}

DifferentialTowerRecognitionResult recognizeDifferentialTower(
    const expression::Expr& expression,
    const expression::Symbol& baseVariable,
    symbols::SymbolTable& symbols,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    DifferentialTowerRecognitionOptions options) {
    TowerRecognizer recognizer{
        baseVariable, symbols, builtins, mathematics, angles,
        options.maximumTowerDepth, expression};
    Expr rewritten = recognizer.rewrite(expression);
    const bool depthLimit = recognizer.depthLimitReached();
    DifferentialTower tower = recognizer.takeTower();
    DifferentialTowerRecognitionStatus status;
    if (!tower.valid())
        status = DifferentialTowerRecognitionStatus::InvalidBaseVariable;
    else if (depthLimit)
        status = DifferentialTowerRecognitionStatus::ResourceLimit;
    else if (rationalMembership(rewritten, tower, builtins))
        status = DifferentialTowerRecognitionStatus::Complete;
    else
        status = DifferentialTowerRecognitionStatus::Partial;
    return {std::move(tower), std::move(rewritten), status};
}

bool isRationalInDifferentialTower(
    const expression::Expr& expression,
    const DifferentialTower& tower,
    const evaluation::BuiltinRegistry& builtins) {
    return tower.valid() && rationalMembership(expression, tower, builtins);
}

} // namespace mmcal::symbolic::risch
