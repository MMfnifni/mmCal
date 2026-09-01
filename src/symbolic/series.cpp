#include "series.hpp"

#include "builtins/arithmetic.hpp"
#include "builtins/special_functions.hpp"
#include "expression/array_utils.hpp"
#include "mathematics/value_facts.hpp"
#include "mathematics/knowledge_context.hpp"
#include "numeric/integer_algorithms.hpp"
#include "simplification/simplifier.hpp"
#include "solver/solution_set.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace mmcal::symbolic {
namespace {

using expression::Expr;
using evaluation::BuiltinId;

constexpr std::int64_t kMaximumInternalExponent = 8192;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{numeric::Number{numeric::BigInt{value}}};
}

[[nodiscard]] Expr pi(const mathematics::MathRegistry& mathematics) {
    if (const auto* definition = mathematics.findConstant(mathematics::ConstantId::Pi))
        return Expr{definition->symbol};
    return Expr{expression::Symbol{"Pi"}};
}

[[nodiscard]] Expr e(const mathematics::MathRegistry& mathematics) {
    if (const auto* definition = mathematics.findConstant(mathematics::ConstantId::E))
        return Expr{definition->symbol};
    return Expr{expression::Symbol{"E"}};
}

[[nodiscard]] Expr call(
    BuiltinId id,
    std::vector<Expr> arguments,
    const evaluation::BuiltinRegistry& builtins) {
    return Expr::call(builtins.symbol(id), std::move(arguments));
}

[[nodiscard]] std::optional<mathematics::AngleUnit> explicitAngleUnit(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall()) return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::UnitApplied)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() != 2 || !arguments[1].isString())
        return std::nullopt;
    return mathematics::AngleSemantics::parseUnit(arguments[1].asString());
}

[[nodiscard]] bool isZero(const Expr& expression) noexcept {
    return expression.isNumber() && expression.asNumber().isZero();
}

[[nodiscard]] bool isOne(const Expr& expression) noexcept {
    return expression.isNumber()
        && expression.asNumber().isReal()
        && expression.asNumber().asReal().isInteger()
        && expression.asNumber().asReal().asInteger() == numeric::BigInt{1};
}

[[nodiscard]] bool isPositiveInfinityCenter(const Expr& expression) noexcept {
    return expression.isSymbol() && expression.asSymbol().view() == "Infinity";
}

[[nodiscard]] Expr canonicalizeInfinityReciprocals(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall()) return expression;

    const auto& source = expression.asCall();
    std::vector<Expr> arguments;
    arguments.reserve(source.arguments.size());
    for (const Expr& argument : source.arguments)
        arguments.push_back(canonicalizeInfinityReciprocals(argument, builtins));

    Expr rebuilt = Expr::rebuildCall(source, std::move(arguments));
    const auto* definition = builtins.find(rebuilt.asCall().head);
    if (!definition || definition->id != BuiltinId::Divide) return rebuilt;

    const auto& outer = rebuilt.asCall().arguments;
    if (outer.size() != 2 || !isOne(outer[0]) || !outer[1].isCall()) return rebuilt;
    const auto* innerDefinition = builtins.find(outer[1].asCall().head);
    if (!innerDefinition || innerDefinition->id != BuiltinId::Divide) return rebuilt;

    const auto& inner = outer[1].asCall().arguments;
    if (inner.size() != 2 || !isOne(inner[0])) return rebuilt;
    return inner[1];
}

[[nodiscard]] Expr simplify(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    return simplification::Simplifier{}.simplify(
        std::move(expression),
        simplification::SimplificationContext{builtins, mathematics, angles, assumptions});
}

[[nodiscard]] Expr specialFunctionValue(
    BuiltinId id,
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const std::array<Expr, 1> arguments{argument};
    return simplify(
        builtins::evaluateSpecialFunction(
            id, arguments, builtins, mathematics, angles),
        builtins, mathematics, angles, assumptions);
}


[[nodiscard]] Expr add(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (isZero(lhs)) return rhs;
    if (isZero(rhs)) return lhs;
    const std::array<Expr, 2> arguments{lhs, rhs};
    return simplify(
        builtins::evaluateAdd(arguments, builtins),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr subtract(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (isZero(rhs)) return lhs;
    const std::array<Expr, 2> arguments{lhs, rhs};
    return simplify(
        builtins::evaluateSubtract(arguments, builtins),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr multiply(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (isZero(lhs) || isZero(rhs)) return integer(0);
    if (isOne(lhs)) return rhs;
    if (isOne(rhs)) return lhs;
    const std::array<Expr, 2> arguments{lhs, rhs};
    return simplify(
        builtins::evaluateMultiply(arguments, builtins),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr divide(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (isZero(lhs)) return integer(0);
    if (isOne(rhs)) return lhs;
    const std::array<Expr, 2> arguments{lhs, rhs};
    return simplify(
        builtins::evaluateDivide(arguments, builtins),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr directTrigScale(
    mathematics::AngleUnit unit,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    switch (unit) {
    case mathematics::AngleUnit::Radian:
        return integer(1);
    case mathematics::AngleUnit::Degree:
        return divide(pi(mathematics), integer(180),
            builtins, mathematics, angles, assumptions);
    case mathematics::AngleUnit::Gradian:
        return divide(pi(mathematics), integer(200),
            builtins, mathematics, angles, assumptions);
    }
    return integer(1);
}

[[nodiscard]] Expr inverseTrigScale(
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Radian:
        return integer(1);
    case mathematics::AngleUnit::Degree:
        return divide(integer(180), pi(mathematics),
            builtins, mathematics, angles, assumptions);
    case mathematics::AngleUnit::Gradian:
        return divide(integer(200), pi(mathematics),
            builtins, mathematics, angles, assumptions);
    }
    return integer(1);
}

[[nodiscard]] Expr radianTrigArgument(
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    mathematics::AngleUnit unit = angles.defaultUnit();
    Expr value = argument;
    if (const auto explicitUnit = explicitAngleUnit(argument, builtins)) {
        unit = *explicitUnit;
        value = argument.asCall().arguments[0];
    }
    return multiply(
        directTrigScale(unit, builtins, mathematics, angles, assumptions), value,
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr radianTrigValue(
    BuiltinId id,
    const Expr& radianArgument,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr unit = call(BuiltinId::UnitApplied,
        {radianArgument, Expr{std::string{"Rad"}}}, builtins);
    return simplify(
        call(id, {std::move(unit)}, builtins),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<Expr> lowCostSeriesRewrite(
    BuiltinId id,
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (id == BuiltinId::Log && arguments.size() == 2) {
        return call(BuiltinId::Divide, {
            call(BuiltinId::Log, {arguments[1]}, builtins),
            call(BuiltinId::Log, {arguments[0]}, builtins)}, builtins);
    }
    if (arguments.size() != 1) return std::nullopt;
    const Expr& x = arguments.front();
    const auto unary = [&](BuiltinId head, Expr argument) {
        return call(head, {std::move(argument)}, builtins);
    };

    switch (id) {
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc:
    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc: {
        const Expr angle = radianTrigArgument(
            x, builtins, mathematics, angles, assumptions);
        const Expr sine = radianTrigValue(
            BuiltinId::Sin, angle, builtins, mathematics, angles, assumptions);
        const Expr cosine = radianTrigValue(
            BuiltinId::Cos, angle, builtins, mathematics, angles, assumptions);
        if (id == BuiltinId::Tan)
            return call(BuiltinId::Divide, {sine, cosine}, builtins);
        if (id == BuiltinId::Cot)
            return call(BuiltinId::Divide, {cosine, sine}, builtins);
        if (id == BuiltinId::Sec)
            return call(BuiltinId::Divide, {integer(1), cosine}, builtins);
        if (id == BuiltinId::Csc)
            return call(BuiltinId::Divide, {integer(1), sine}, builtins);
        if (id == BuiltinId::Sinc)
            return call(BuiltinId::Divide, {sine, angle}, builtins);
        if (id == BuiltinId::Cosc) {
            const Expr halfAngle = call(
                BuiltinId::Divide, {angle, integer(2)}, builtins);
            const Expr halfSine = radianTrigValue(
                BuiltinId::Sin, halfAngle, builtins, mathematics, angles, assumptions);
            return call(BuiltinId::Divide, {
                call(BuiltinId::Multiply, {
                    integer(2), call(BuiltinId::Power, {halfSine, integer(2)}, builtins)}, builtins),
                angle}, builtins);
        }
        return call(BuiltinId::Divide, {
            call(BuiltinId::Divide, {sine, cosine}, builtins), angle}, builtins);
    }
    case BuiltinId::Tanh:
        return call(BuiltinId::Divide, {
            unary(BuiltinId::Sinh, x), unary(BuiltinId::Cosh, x)}, builtins);
    case BuiltinId::Coth:
        return call(BuiltinId::Divide, {
            unary(BuiltinId::Cosh, x), unary(BuiltinId::Sinh, x)}, builtins);
    case BuiltinId::Sech:
        return call(BuiltinId::Divide, {integer(1), unary(BuiltinId::Cosh, x)}, builtins);
    case BuiltinId::Csch:
        return call(BuiltinId::Divide, {integer(1), unary(BuiltinId::Sinh, x)}, builtins);
    case BuiltinId::Expm1:
        return call(BuiltinId::Subtract, {unary(BuiltinId::Exp, x), integer(1)}, builtins);
    case BuiltinId::Log1p:
        return unary(BuiltinId::Log, call(BuiltinId::Add, {integer(1), x}, builtins));
    case BuiltinId::Sinhc:
        return call(BuiltinId::Divide, {unary(BuiltinId::Sinh, x), x}, builtins);
    case BuiltinId::Tanhc:
        return call(BuiltinId::Divide, {
            call(BuiltinId::Divide, {unary(BuiltinId::Sinh, x), unary(BuiltinId::Cosh, x)}, builtins), x}, builtins);
    case BuiltinId::Expc:
        return call(BuiltinId::Divide, {
            unary(BuiltinId::Expm1, x), x}, builtins);
    case BuiltinId::Log2:
    case BuiltinId::Log10: {
        const Expr base = integer(id == BuiltinId::Log2 ? 2 : 10);
        return call(BuiltinId::Divide, {
            unary(BuiltinId::Log, x), unary(BuiltinId::Log, base)}, builtins);
    }
    default:
        return std::nullopt;
    }
}

[[nodiscard]] Expr negate(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (isZero(value)) return value;
    const std::array<Expr, 1> arguments{value};
    return simplify(
        builtins::evaluateNegate(arguments, builtins),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<std::int64_t> exactInt64(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const std::string text = expression.asNumber().asReal().asInteger().toString();
    std::int64_t value = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size())
        return std::nullopt;
    return value;
}

[[nodiscard]] std::optional<numeric::Rational> exactRational(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    return expression.asNumber().asReal().toRational();
}

struct SmallRational final {
    std::int64_t numerator = 0;
    std::uint32_t denominator = 1;
};

[[nodiscard]] std::optional<SmallRational> smallRational(const Expr& expression) {
    const auto value = exactRational(expression);
    if (!value) return std::nullopt;
    const std::string numeratorText = value->numerator().toString();
    const std::string denominatorText = value->denominator().toString();
    std::int64_t numerator = 0;
    std::uint64_t denominator = 0;
    const auto nr = std::from_chars(
        numeratorText.data(), numeratorText.data() + numeratorText.size(), numerator);
    const auto dr = std::from_chars(
        denominatorText.data(), denominatorText.data() + denominatorText.size(), denominator);
    if (nr.ec != std::errc{} || nr.ptr != numeratorText.data() + numeratorText.size()
        || dr.ec != std::errc{} || dr.ptr != denominatorText.data() + denominatorText.size()
        || denominator == 0 || denominator > 64)
        return std::nullopt;
    return SmallRational{numerator, static_cast<std::uint32_t>(denominator)};
}

[[nodiscard]] Expr rational(std::int64_t numerator, std::int64_t denominator) {
    return Expr{numeric::Number{numeric::Rational{
        numeric::BigInt{numerator}, numeric::BigInt{denominator}}}};
}

[[nodiscard]] std::optional<std::uint32_t> exactUint32(const Expr& expression) {
    const auto value = exactInt64(expression);
    if (!value || *value <= 0
        || static_cast<std::uint64_t>(*value) > std::numeric_limits<std::uint32_t>::max())
        return std::nullopt;
    return static_cast<std::uint32_t>(*value);
}

[[nodiscard]] std::optional<std::vector<Expr>> braceElements(const Expr& expression) {
    if (expression.isArray())
        return expression.asArray().materialize();
    if (expression.isList())
        return expression.asList().elements;
    return std::nullopt;
}

[[nodiscard]] std::optional<BuiltinId> builtinId(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) noexcept {
    if (!expression.isCall()) return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition) return std::nullopt;
    return definition->id;
}

[[nodiscard]] bool provablyNonZero(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    if (expression.isNumber())
        return !expression.asNumber().isZero();
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::relation(
            mathematics::RelationKind::NotEqual, expression, integer(0)))
        == mathematics::TruthValue::True)
        return true;
    const auto facts = knowledge.facts(expression);
    return facts.provablyNonReal
        || facts.sign == mathematics::RealSign::Positive
        || facts.sign == mathematics::RealSign::Negative
        || facts.sign == mathematics::RealSign::NonZero;
}

[[nodiscard]] bool principalPowerAnalyticAt(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(value);
    return facts.provablyNonReal || facts.sign == mathematics::RealSign::Positive;
}

[[nodiscard]] bool logarithmicIntegralAnalyticAt(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    // principal li(z)=Ei(Log(z)) として，Logの切断とEiの切断を同時に避ける。
    // 現段階ではDLMF 6.2.8の実領域x>1，または非実中心だけを安全な局所枝として採用する。
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(value);
    if (facts.provablyNonReal)
        return true;
    return knowledge.prove(mathematics::relation(
               mathematics::RelationKind::Greater, value, integer(1)))
        == mathematics::TruthValue::True;
}

[[nodiscard]] bool polylogRegularCenterAnalyticAt(
    const Expr& order,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const auto integerOrder = exactInt64(order);
    if (!integerOrder || *integerOrder <= 0)
        return false;

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(value);
    if (facts.provablyNonReal)
        return true;
    return knowledge.prove(mathematics::relation(
               mathematics::RelationKind::Less, value, integer(1)))
        == mathematics::TruthValue::True;
}

[[nodiscard]] bool polylogPositiveRealCenter(
    const Expr& order,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    const auto integerOrder = exactInt64(order);
    if (!integerOrder || *integerOrder <= 0)
        return false;

    if (const auto exactCenter = exactRational(value)) {
        const numeric::Rational zero{numeric::BigInt{0}};
        const numeric::Rational one{numeric::BigInt{1}};
        return *exactCenter > zero && *exactCenter < one;
    }

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(value);
    const bool positive = facts.sign == mathematics::RealSign::Positive
        || knowledge.prove(mathematics::relation(
               mathematics::RelationKind::Greater, value, integer(0)))
            == mathematics::TruthValue::True;
    return positive
        && knowledge.prove(mathematics::relation(
               mathematics::RelationKind::Less, value, integer(1)))
            == mathematics::TruthValue::True;
}

[[nodiscard]] bool seriesCenterValueProvablyNonZero(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr centerValue = simplify(
        substituteSymbol(expression, variable, center),
        builtins, mathematics, angles, assumptions);
    if (provablyNonZero(centerValue, builtins, mathematics, assumptions))
        return true;

    const auto id = builtinId(centerValue, builtins);
    if (!id) return false;
    const auto& arguments = centerValue.asCall().arguments;
    if (*id == BuiltinId::Log && arguments.size() == 1) {
        if (const auto value = exactRational(arguments[0])) {
            const numeric::Rational zero{numeric::BigInt{0}};
            const numeric::Rational one{numeric::BigInt{1}};
            return *value > zero && *value != one;
        }
        return false;
    }
    return *id == BuiltinId::Polylog && arguments.size() == 2
        && polylogPositiveRealCenter(
            arguments[0], arguments[1], builtins, mathematics, assumptions);
}

[[nodiscard]] std::optional<std::int64_t> lambertWBranchIndex(
    std::span<const Expr> arguments) {
    if (arguments.size() == 1)
        return 0;
    if (arguments.size() != 2)
        return std::nullopt;
    return exactInt64(arguments[0]);
}

[[nodiscard]] std::optional<bool> lambertWPrincipalBranch(
    std::span<const Expr> arguments) {
    const auto branch = lambertWBranchIndex(arguments);
    if (!branch) return std::nullopt;
    return *branch == 0;
}

[[nodiscard]] bool isLambertWBranchPoint(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const Expr branchPoint = divide(
        integer(-1), e(mathematics), builtins, mathematics, angles, assumptions);
    return isZero(subtract(
        value, branchPoint, builtins, mathematics, angles, assumptions));
}

[[nodiscard]] bool lambertWAnalyticAt(
    bool principal,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    // DLMF 4.13に従い，W_0は(-∞,-1/E]，その他の枝は(-∞,0]を局所Taylor展開から除外する。
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(value);
    if (facts.provablyNonReal)
        return true;
    if (principal && (facts.sign == mathematics::RealSign::Zero
            || facts.sign == mathematics::RealSign::Positive
            || facts.sign == mathematics::RealSign::NonNegative))
        return true;
    if (principal) {
        const Expr branchPoint = divide(
            integer(-1), e(mathematics), builtins, mathematics, angles, assumptions);
        return knowledge.prove(mathematics::relation(
                   mathematics::RelationKind::Greater, value, branchPoint))
            == mathematics::TruthValue::True;
    }
    return knowledge.prove(mathematics::relation(
               mathematics::RelationKind::Greater, value, integer(0)))
        == mathematics::TruthValue::True;
}

enum class GammaSeriesBase { One, Half };

struct GammaSeriesPlan final {
    GammaSeriesBase base = GammaSeriesBase::One;
    std::int64_t shift = 0;
    bool positiveAtCenter = true;
};

[[nodiscard]] std::optional<GammaSeriesPlan> gammaSeriesPlan(const Expr& center) {
    constexpr std::int64_t kMaximumShift = 256;
    const auto value = smallRational(center);
    if (!value) return std::nullopt;

    if (value->denominator == 1) {
        if (value->numerator <= 0 || value->numerator - 1 > kMaximumShift)
            return std::nullopt;
        return GammaSeriesPlan{
            GammaSeriesBase::One, value->numerator - 1, true};
    }

    if (value->denominator != 2 || (value->numerator & 1) == 0)
        return std::nullopt;
    const std::int64_t shift = (value->numerator - 1) / 2;
    if (std::abs(shift) > kMaximumShift)
        return std::nullopt;
    const bool positive = shift >= 0 || ((-shift) % 2 == 0);
    return GammaSeriesPlan{GammaSeriesBase::Half, shift, positive};
}

[[nodiscard]] bool logGammaZeroAtSupportedCenter(const Expr& center) {
    const auto value = exactRational(center);
    if (!value || !value->isInteger()) return false;
    return value->numerator() == numeric::BigInt{1}
        || value->numerator() == numeric::BigInt{2};
}

[[nodiscard]] bool inverseTrigAnalyticAt(
    BuiltinId id,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    // asin/acosのprincipal branchは実軸の(-∞,-1]と[1,∞)を切断する。
    // atanは虚軸上の±Iから外側を切断するため，実部が非零なら正則である。
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(value);
    if (id == BuiltinId::Asin || id == BuiltinId::Acos) {
        if (facts.provablyNonReal)
            return true;
        return knowledge.prove(mathematics::relation(
                   mathematics::RelationKind::Greater, value, integer(-1)))
                == mathematics::TruthValue::True
            && knowledge.prove(mathematics::relation(
                   mathematics::RelationKind::Less, value, integer(1)))
                == mathematics::TruthValue::True;
    }
    if (id != BuiltinId::Atan)
        return false;
    if (facts.isProvablyReal())
        return true;

    const Expr realPart = simplify(
        call(BuiltinId::Re, {value}, builtins),
        builtins, mathematics, angles, assumptions);
    if (provablyNonZero(realPart, builtins, mathematics, assumptions))
        return true;
    if (!isZero(realPart))
        return false;

    const Expr imaginaryPart = simplify(
        call(BuiltinId::Im, {value}, builtins),
        builtins, mathematics, angles, assumptions);
    return knowledge.prove(mathematics::relation(
               mathematics::RelationKind::Greater, imaginaryPart, integer(-1)))
            == mathematics::TruthValue::True
        && knowledge.prove(mathematics::relation(
               mathematics::RelationKind::Less, imaginaryPart, integer(1)))
            == mathematics::TruthValue::True;
}

struct Valuation final {
    enum class State { Finite, Zero, Unknown } state = State::Unknown;
    std::int64_t exponent = 0;

    [[nodiscard]] static Valuation finite(std::int64_t exponent) noexcept {
        return Valuation{State::Finite, exponent};
    }
    [[nodiscard]] static Valuation zero() noexcept {
        return Valuation{State::Zero, 0};
    }
    [[nodiscard]] static Valuation unknown() noexcept {
        return Valuation{State::Unknown, 0};
    }
};

[[nodiscard]] expression::Symbol temporarySeriesVariable(const Expr& expression) {
    for (std::size_t suffix = 0; suffix < 1024; ++suffix) {
        const std::string name = suffix == 0
            ? "__mmcal_series_t"
            : "__mmcal_series_t" + std::to_string(suffix);
        expression::Symbol candidate{name};
        if (!containsSymbol(expression, candidate))
            return candidate;
    }
    return expression::Symbol{"__mmcal_series_t_fallback"};
}

[[nodiscard]] std::optional<std::int64_t> polynomialValuation(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const expression::Symbol temporary = temporarySeriesVariable(expression);
    Expr shiftedVariable = Expr{temporary};
    if (!isZero(center)) {
        const std::array<Expr, 2> arguments{Expr{temporary}, center};
        shiftedVariable = builtins::evaluateAdd(arguments, builtins);
    }
    const Expr shifted = substituteSymbol(expression, variable, shiftedVariable);
    const auto polynomial = toExpressionPolynomial(
        shifted, temporary, builtins, mathematics, angles,
        PolynomialConversionOptions{4096, 4096});
    if (!polynomial)
        return std::nullopt;
    if (polynomial->isZero())
        return std::numeric_limits<std::int64_t>::max();
    for (std::size_t exponent = 0; exponent <= polynomial->degree(); ++exponent)
        if (!isZero(polynomial->coefficient(exponent)))
            return static_cast<std::int64_t>(exponent);
    return std::numeric_limits<std::int64_t>::max();
}

[[nodiscard]] Valuation structuralValuation(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!containsSymbol(expression, variable))
        return isZero(expression) ? Valuation::zero() : Valuation::finite(0);

    // 多項式部分はx=center+tへexactに平行移動して最低次数を読む。
    // これによりx-aや加法相殺を個別規則なしで扱える。
    if (const auto polynomial = polynomialValuation(
            expression, variable, center, builtins, mathematics, angles)) {
        if (*polynomial == std::numeric_limits<std::int64_t>::max())
            return Valuation::zero();
        return Valuation::finite(*polynomial);
    }

    if (expression.isSymbol() && expression.asSymbol() == variable)
        return isZero(center) ? Valuation::finite(1) : Valuation::finite(0);

    if (!expression.isCall())
        return Valuation::unknown();
    const auto id = builtinId(expression, builtins);
    if (!id) return Valuation::unknown();
    const auto& arguments = expression.asCall().arguments;

    const auto childValuation = [&](const Expr& child) {
        return structuralValuation(
            child, variable, center, builtins, mathematics, angles, assumptions);
    };

    switch (*id) {
    case BuiltinId::Negate:
        return arguments.size() == 1 ? childValuation(arguments[0]) : Valuation::unknown();
    case BuiltinId::Add:
    case BuiltinId::Subtract: {
        if (arguments.empty()) return Valuation::zero();
        std::optional<std::int64_t> minimum;
        std::size_t minimumCount = 0;
        for (const Expr& argument : arguments) {
            const auto value = childValuation(argument);
            if (value.state == Valuation::State::Unknown)
                return Valuation::unknown();
            if (value.state == Valuation::State::Zero)
                continue;
            if (!minimum || value.exponent < *minimum) {
                minimum = value.exponent;
                minimumCount = 1;
            }
            else if (value.exponent == *minimum) {
                ++minimumCount;
            }
        }
        if (!minimum) return Valuation::zero();
        return minimumCount == 1 ? Valuation::finite(*minimum) : Valuation::unknown();
    }
    case BuiltinId::Multiply: {
        std::int64_t total = 0;
        for (const Expr& argument : arguments) {
            const auto value = childValuation(argument);
            if (value.state == Valuation::State::Unknown)
                return Valuation::unknown();
            if (value.state == Valuation::State::Zero)
                return Valuation::zero();
            if ((value.exponent > 0 && total > kMaximumInternalExponent - value.exponent)
                || (value.exponent < 0 && total < -kMaximumInternalExponent - value.exponent))
                return Valuation::unknown();
            total += value.exponent;
        }
        return Valuation::finite(total);
    }
    case BuiltinId::Divide: {
        if (arguments.size() != 2) return Valuation::unknown();
        const auto numerator = childValuation(arguments[0]);
        const auto denominator = childValuation(arguments[1]);
        if (numerator.state == Valuation::State::Zero)
            return Valuation::zero();
        if (numerator.state != Valuation::State::Finite
            || denominator.state != Valuation::State::Finite)
            return Valuation::unknown();
        return Valuation::finite(numerator.exponent - denominator.exponent);
    }
    case BuiltinId::Power: {
        if (arguments.size() != 2) return Valuation::unknown();
        const auto exponent = exactInt64(arguments[1]);
        if (exponent) {
            if (*exponent == 0) return Valuation::finite(0);
            const auto base = childValuation(arguments[0]);
            if (base.state == Valuation::State::Zero)
                return *exponent > 0 ? Valuation::zero() : Valuation::unknown();
            if (base.state != Valuation::State::Finite)
                return Valuation::unknown();
            if (base.exponent != 0
                && std::abs(*exponent) > kMaximumInternalExponent / std::abs(base.exponent))
                return Valuation::unknown();
            return Valuation::finite(base.exponent * *exponent);
        }

        if (!exactRational(arguments[1])) return Valuation::unknown();
        Expr centerBase = simplify(
            substituteSymbol(arguments[0], variable, center),
            builtins, mathematics, angles, assumptions);
        if (!principalPowerAnalyticAt(centerBase, builtins, mathematics, assumptions))
            return Valuation::unknown();
        return Valuation::finite(0);
    }
    case BuiltinId::Sqrt: {
        if (arguments.size() != 1) return Valuation::unknown();
        Expr centerBase = simplify(
            substituteSymbol(arguments[0], variable, center),
            builtins, mathematics, angles, assumptions);
        if (!principalPowerAnalyticAt(centerBase, builtins, mathematics, assumptions))
            return Valuation::unknown();
        return Valuation::finite(0);
    }
    case BuiltinId::Expm1:
    case BuiltinId::Log1p: {
        if (arguments.size() != 1) return Valuation::unknown();
        Expr centerArgument = simplify(
            substituteSymbol(arguments[0], variable, center),
            builtins, mathematics, angles, assumptions);
        if (isZero(centerArgument)) {
            Expr variation = subtract(
                arguments[0], centerArgument,
                builtins, mathematics, angles, assumptions);
            return structuralValuation(
                variation, variable, center, builtins, mathematics, angles, assumptions);
        }
        return Valuation::unknown();
    }
    case BuiltinId::Exp:
    case BuiltinId::Log:
    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Sinh:
    case BuiltinId::Cosh: {
        if (*id == BuiltinId::Log && arguments.size() == 2) {
            auto rewritten = lowCostSeriesRewrite(
                *id, arguments, builtins, mathematics, angles, assumptions);
            if (!rewritten) return Valuation::unknown();
            return structuralValuation(
                *rewritten, variable, center, builtins, mathematics, angles, assumptions);
        }
        if (arguments.size() != 1) return Valuation::unknown();

        const bool directTrig = *id == BuiltinId::Sin || *id == BuiltinId::Cos;
        Expr analyticArgument = directTrig
            ? radianTrigArgument(arguments[0], builtins, mathematics, angles, assumptions)
            : arguments[0];
        Expr centerArgument = simplify(
            substituteSymbol(analyticArgument, variable, center),
            builtins, mathematics, angles, assumptions);

        if (*id == BuiltinId::Log) {
            const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
            const auto facts = knowledge.facts(centerArgument);
            if (!(facts.provablyNonReal || facts.sign == mathematics::RealSign::Positive))
                return Valuation::unknown();
        }

        Expr centerValue = integer(0);
        if (directTrig)
            centerValue = radianTrigValue(
                *id, centerArgument, builtins, mathematics, angles, assumptions);
        else
            centerValue = simplify(
                call(*id, {centerArgument}, builtins),
                builtins, mathematics, angles, assumptions);

        if (provablyNonZero(centerValue, builtins, mathematics, assumptions))
            return Valuation::finite(0);
        if (!isZero(centerValue))
            return Valuation::unknown();

        Expr derivativeValue = integer(0);
        switch (*id) {
        case BuiltinId::Exp:
            derivativeValue = centerValue;
            break;
        case BuiltinId::Log:
            derivativeValue = divide(
                integer(1), centerArgument,
                builtins, mathematics, angles, assumptions);
            break;
        case BuiltinId::Sin:
            derivativeValue = radianTrigValue(
                BuiltinId::Cos, centerArgument,
                builtins, mathematics, angles, assumptions);
            break;
        case BuiltinId::Cos:
            derivativeValue = negate(
                radianTrigValue(BuiltinId::Sin, centerArgument,
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
            break;
        case BuiltinId::Sinh:
            derivativeValue = simplify(
                call(BuiltinId::Cosh, {centerArgument}, builtins),
                builtins, mathematics, angles, assumptions);
            break;
        case BuiltinId::Cosh:
            derivativeValue = simplify(
                call(BuiltinId::Sinh, {centerArgument}, builtins),
                builtins, mathematics, angles, assumptions);
            break;
        default:
            return Valuation::unknown();
        }
        if (!provablyNonZero(derivativeValue, builtins, mathematics, assumptions))
            return Valuation::unknown();

        Expr variation = subtract(
            analyticArgument, centerArgument,
            builtins, mathematics, angles, assumptions);
        return structuralValuation(
            variation, variable, center, builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Gamma:
    case BuiltinId::LogGamma:
    case BuiltinId::Digamma:
    case BuiltinId::Trigamma: {
        if (arguments.size() != 1) return Valuation::unknown();
        Expr centerArgument = simplify(
            substituteSymbol(arguments[0], variable, center),
            builtins, mathematics, angles, assumptions);
        const auto plan = gammaSeriesPlan(centerArgument);
        if (!plan) return Valuation::unknown();
        if (*id == BuiltinId::Gamma || *id == BuiltinId::Digamma
            || *id == BuiltinId::Trigamma)
            return Valuation::finite(0);
        if (!logGammaZeroAtSupportedCenter(centerArgument))
            return Valuation::finite(0);

        Expr variation = subtract(
            arguments[0], centerArgument,
            builtins, mathematics, angles, assumptions);
        return structuralValuation(
            variation, variable, center, builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::LambertW: {
        const auto principal = lambertWPrincipalBranch(arguments);
        if (!principal) return Valuation::unknown();
        const Expr& argument = arguments.back();
        Expr centerArgument = simplify(
            substituteSymbol(argument, variable, center),
            builtins, mathematics, angles, assumptions);
        const auto branch = lambertWBranchIndex(arguments);
        if (!branch) return Valuation::unknown();
        if ((*branch == 0 || *branch == -1)
            && isLambertWBranchPoint(
                centerArgument, builtins, mathematics, angles, assumptions))
            return Valuation::finite(0);
        if (!lambertWAnalyticAt(
                *principal, centerArgument, builtins, mathematics, angles, assumptions))
            return Valuation::unknown();

        if (!*principal)
            return Valuation::finite(0);
        if (isZero(centerArgument)) {
            Expr variation = subtract(
                argument, centerArgument,
                builtins, mathematics, angles, assumptions);
            return structuralValuation(
                variation, variable, center, builtins, mathematics, angles, assumptions);
        }
        if (provablyNonZero(centerArgument, builtins, mathematics, assumptions))
            return Valuation::finite(0);
        return Valuation::unknown();
    }
    case BuiltinId::Asin:
    case BuiltinId::Acos:
    case BuiltinId::Atan: {
        if (arguments.size() != 1) return Valuation::unknown();
        Expr centerArgument = simplify(
            substituteSymbol(arguments[0], variable, center),
            builtins, mathematics, angles, assumptions);
        if (!inverseTrigAnalyticAt(
                *id, centerArgument, builtins, mathematics, angles, assumptions))
            return Valuation::unknown();

        Expr centerValue = simplify(
            call(*id, {centerArgument}, builtins),
            builtins, mathematics, angles, assumptions);
        if (provablyNonZero(centerValue, builtins, mathematics, assumptions))
            return Valuation::finite(0);
        if (!isZero(centerValue))
            return Valuation::unknown();

        Expr variation = subtract(
            arguments[0], centerArgument,
            builtins, mathematics, angles, assumptions);
        return structuralValuation(
            variation, variable, center, builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Polylog: {
        if (arguments.size() != 2 || containsSymbol(arguments[0], variable))
            return Valuation::unknown();
        Expr centerArgument = simplify(
            substituteSymbol(arguments[1], variable, center),
            builtins, mathematics, angles, assumptions);
        if (isZero(centerArgument)) {
            Expr variation = subtract(
                arguments[1], centerArgument,
                builtins, mathematics, angles, assumptions);
            return structuralValuation(
                variation, variable, center, builtins, mathematics, angles, assumptions);
        }
        if (!polylogRegularCenterAnalyticAt(
                arguments[0], centerArgument, builtins, mathematics, assumptions))
            return Valuation::unknown();

        Expr centerValue = simplify(
            call(BuiltinId::Polylog, {arguments[0], centerArgument}, builtins),
            builtins, mathematics, angles, assumptions);
        if (provablyNonZero(centerValue, builtins, mathematics, assumptions)
            || polylogPositiveRealCenter(
                arguments[0], centerArgument, builtins, mathematics, assumptions))
            return Valuation::finite(0);
        return Valuation::unknown();
    }
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi: {
        if (arguments.size() != 1) return Valuation::unknown();

        Expr centerArgument = simplify(
            substituteSymbol(arguments[0], variable, center),
            builtins, mathematics, angles, assumptions);
        if ((*id == BuiltinId::ExponentialIntegralEi
                || *id == BuiltinId::CosineIntegralCi)
            && !principalPowerAnalyticAt(
                centerArgument, builtins, mathematics, assumptions))
            return Valuation::unknown();
        if (*id == BuiltinId::LogarithmicIntegralLi
            && !logarithmicIntegralAnalyticAt(
                centerArgument, builtins, mathematics, assumptions))
            return Valuation::unknown();

        Expr centerValue = integer(0);
        if (isZero(centerArgument)) {
            if (*id == BuiltinId::Erfc)
                centerValue = integer(1);
            else if (*id == BuiltinId::ExponentialIntegralEi
                || *id == BuiltinId::CosineIntegralCi
                || *id == BuiltinId::LogarithmicIntegralLi)
                return Valuation::unknown();
        }
        else {
            centerValue = simplify(
                call(*id, {centerArgument}, builtins),
                builtins, mathematics, angles, assumptions);
        }
        if (provablyNonZero(centerValue, builtins, mathematics, assumptions))
            return Valuation::finite(0);
        if (!isZero(centerValue))
            return Valuation::unknown();

        // 現在exactに0へ簡約される起点では，既知の最低非零Taylor次数だけを使う。
        // FresnelSはS(z)=Pi*z^3/6+O[z^7]，その他の0起点函数は単純零点である。
        std::int64_t zeroMultiplicity = 1;
        if (*id == BuiltinId::FresnelS)
            zeroMultiplicity = 3;
        else if (*id == BuiltinId::Erfc
            || *id == BuiltinId::ExponentialIntegralEi
            || *id == BuiltinId::CosineIntegralCi
            || *id == BuiltinId::LogarithmicIntegralLi)
            return Valuation::unknown();

        Expr variation = subtract(
            arguments[0], centerArgument,
            builtins, mathematics, angles, assumptions);
        const auto argumentValuation = structuralValuation(
            variation, variable, center, builtins, mathematics, angles, assumptions);
        if (argumentValuation.state != Valuation::State::Finite)
            return argumentValuation;
        if (argumentValuation.exponent != 0
            && std::abs(argumentValuation.exponent)
                > kMaximumInternalExponent / zeroMultiplicity)
            return Valuation::unknown();
        return Valuation::finite(argumentValuation.exponent * zeroMultiplicity);
    }
    default:
        return Valuation::unknown();
    }
}

struct LaurentSeries final {
    std::int64_t minimumExponent = 0;
    std::int64_t ceiling = 0;
    std::vector<Expr> coefficients;
    bool exactZero = false;
};

[[nodiscard]] LaurentSeries zeroSeries(std::int64_t ceiling) {
    return LaurentSeries{0, ceiling, {}, true};
}

[[nodiscard]] Expr coefficientAt(const LaurentSeries& series, std::int64_t exponent) {
    if (series.exactZero || exponent < series.minimumExponent || exponent > series.ceiling)
        return integer(0);
    const auto index = static_cast<std::size_t>(exponent - series.minimumExponent);
    if (index >= series.coefficients.size())
        return integer(0);
    return series.coefficients[index];
}

void trimLeadingZeros(LaurentSeries& series) {
    if (series.exactZero) return;
    std::size_t leading = 0;
    while (leading < series.coefficients.size() && isZero(series.coefficients[leading]))
        ++leading;
    if (leading == series.coefficients.size()) {
        series = zeroSeries(series.ceiling);
        return;
    }
    if (leading != 0) {
        series.minimumExponent += static_cast<std::int64_t>(leading);
        series.coefficients.erase(
            series.coefficients.begin(), series.coefficients.begin() + static_cast<std::ptrdiff_t>(leading));
    }
}

[[nodiscard]] LaurentSeries denseSeries(
    std::int64_t minimumExponent,
    std::int64_t ceiling) {
    if (minimumExponent > ceiling)
        return zeroSeries(ceiling);
    const auto count = static_cast<std::size_t>(ceiling - minimumExponent + 1);
    return LaurentSeries{minimumExponent, ceiling, std::vector<Expr>(count, integer(0)), false};
}

[[nodiscard]] LaurentSeries addSeries(
    const LaurentSeries& lhs,
    const LaurentSeries& rhs,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    bool subtractRight = false) {
    if (lhs.exactZero && rhs.exactZero)
        return zeroSeries(ceiling);
    const std::int64_t minimum = lhs.exactZero ? rhs.minimumExponent
        : rhs.exactZero ? lhs.minimumExponent
        : std::min(lhs.minimumExponent, rhs.minimumExponent);
    if (minimum > ceiling)
        return zeroSeries(ceiling);
    LaurentSeries result = denseSeries(minimum, ceiling);
    for (std::int64_t exponent = minimum; exponent <= ceiling; ++exponent) {
        const Expr left = coefficientAt(lhs, exponent);
        const Expr right = coefficientAt(rhs, exponent);
        result.coefficients[static_cast<std::size_t>(exponent - minimum)] = subtractRight
            ? subtract(left, right, builtins, mathematics, angles, assumptions)
            : add(left, right, builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] LaurentSeries multiplySeries(
    const LaurentSeries& lhs,
    const LaurentSeries& rhs,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (lhs.exactZero || rhs.exactZero)
        return zeroSeries(ceiling);
    const std::int64_t minimum = lhs.minimumExponent + rhs.minimumExponent;
    if (minimum > ceiling)
        return zeroSeries(ceiling);
    LaurentSeries result = denseSeries(minimum, ceiling);
    for (std::int64_t leftExponent = lhs.minimumExponent;
         leftExponent <= lhs.ceiling; ++leftExponent) {
        const Expr left = coefficientAt(lhs, leftExponent);
        if (isZero(left)) continue;
        const std::int64_t rightMaximum = std::min(rhs.ceiling, ceiling - leftExponent);
        for (std::int64_t rightExponent = rhs.minimumExponent;
             rightExponent <= rightMaximum; ++rightExponent) {
            const Expr right = coefficientAt(rhs, rightExponent);
            if (isZero(right)) continue;
            const std::int64_t exponent = leftExponent + rightExponent;
            Expr& destination = result.coefficients[
                static_cast<std::size_t>(exponent - minimum)];
            destination = add(
                destination,
                multiply(left, right, builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
        }
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> inverseSeries(
    const LaurentSeries& base,
    std::int64_t outputCeiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    bool leadingNonZeroKnown = false) {
    if (base.exactZero || base.minimumExponent > base.ceiling)
        return std::nullopt;
    const Expr leading = coefficientAt(base, base.minimumExponent);
    if (!leadingNonZeroKnown
        && !provablyNonZero(leading, builtins, mathematics, assumptions))
        return std::nullopt;

    const std::int64_t resultMinimum = -base.minimumExponent;
    if (resultMinimum > outputCeiling)
        return zeroSeries(outputCeiling);
    const std::int64_t relativeMaximum = outputCeiling - resultMinimum;
    LaurentSeries result = denseSeries(resultMinimum, outputCeiling);

    const Expr first = divide(
        integer(1), leading, builtins, mathematics, angles, assumptions);
    result.coefficients[0] = first;

    for (std::int64_t degree = 1; degree <= relativeMaximum; ++degree) {
        Expr sum = integer(0);
        for (std::int64_t k = 1; k <= degree; ++k) {
            const Expr baseCoefficient = coefficientAt(
                base, base.minimumExponent + k);
            if (isZero(baseCoefficient)) continue;
            const Expr inverseCoefficient = result.coefficients[
                static_cast<std::size_t>(degree - k)];
            sum = add(
                sum,
                multiply(baseCoefficient, inverseCoefficient,
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
        }
        result.coefficients[static_cast<std::size_t>(degree)] = divide(
            negate(sum, builtins, mathematics, angles, assumptions),
            leading,
            builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandExpSeries(
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!argument.exactZero && argument.minimumExponent < 0)
        return std::nullopt;
    if (ceiling < 0)
        return zeroSeries(ceiling);

    LaurentSeries result = denseSeries(0, ceiling);
    const Expr a0 = coefficientAt(argument, 0);
    result.coefficients[0] = simplify(
        call(BuiltinId::Exp, {a0}, builtins),
        builtins, mathematics, angles, assumptions);

    for (std::int64_t n = 1; n <= ceiling; ++n) {
        Expr sum = integer(0);
        for (std::int64_t k = 1; k <= n; ++k) {
            const Expr ak = coefficientAt(argument, k);
            if (isZero(ak)) continue;
            Expr weighted = multiply(
                integer(k), ak, builtins, mathematics, angles, assumptions);
            Expr term = multiply(
                weighted, result.coefficients[static_cast<std::size_t>(n - k)],
                builtins, mathematics, angles, assumptions);
            sum = add(sum, term, builtins, mathematics, angles, assumptions);
        }
        result.coefficients[static_cast<std::size_t>(n)] = divide(
            sum, integer(n), builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandLogSeries(
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (argument.exactZero || argument.minimumExponent != 0 || ceiling < 0)
        return std::nullopt;
    const Expr a0 = coefficientAt(argument, 0);
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const auto facts = knowledge.facts(a0);
    if (!(facts.provablyNonReal || facts.sign == mathematics::RealSign::Positive))
        return std::nullopt;

    auto inverse = inverseSeries(
        argument, std::max<std::int64_t>(0, ceiling - 1),
        builtins, mathematics, angles, assumptions);
    if (!inverse) return std::nullopt;

    LaurentSeries result = denseSeries(0, ceiling);
    result.coefficients[0] = simplify(
        call(BuiltinId::Log, {a0}, builtins),
        builtins, mathematics, angles, assumptions);
    for (std::int64_t n = 1; n <= ceiling; ++n) {
        Expr derivativeCoefficient = integer(0);
        for (std::int64_t k = 1; k <= n; ++k) {
            const Expr ak = coefficientAt(argument, k);
            if (isZero(ak)) continue;
            Expr weighted = multiply(
                integer(k), ak, builtins, mathematics, angles, assumptions);
            Expr term = multiply(
                weighted, coefficientAt(*inverse, n - k),
                builtins, mathematics, angles, assumptions);
            derivativeCoefficient = add(
                derivativeCoefficient, term,
                builtins, mathematics, angles, assumptions);
        }
        result.coefficients[static_cast<std::size_t>(n)] = divide(
            derivativeCoefficient, integer(n),
            builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<std::pair<LaurentSeries, LaurentSeries>> expandSinCosSeries(
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!argument.exactZero && argument.minimumExponent < 0)
        return std::nullopt;
    if (ceiling < 0)
        return std::pair{zeroSeries(ceiling), zeroSeries(ceiling)};

    LaurentSeries sine = denseSeries(0, ceiling);
    LaurentSeries cosine = denseSeries(0, ceiling);
    const Expr a0 = coefficientAt(argument, 0);
    sine.coefficients[0] = radianTrigValue(
        BuiltinId::Sin, a0, builtins, mathematics, angles, assumptions);
    cosine.coefficients[0] = radianTrigValue(
        BuiltinId::Cos, a0, builtins, mathematics, angles, assumptions);

    for (std::int64_t n = 1; n <= ceiling; ++n) {
        Expr sineSum = integer(0);
        Expr cosineSum = integer(0);
        for (std::int64_t k = 1; k <= n; ++k) {
            const Expr ak = coefficientAt(argument, k);
            if (isZero(ak)) continue;
            Expr weighted = multiply(
                integer(k), ak, builtins, mathematics, angles, assumptions);
            sineSum = add(
                sineSum,
                multiply(weighted, cosine.coefficients[static_cast<std::size_t>(n - k)],
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
            cosineSum = add(
                cosineSum,
                multiply(weighted, sine.coefficients[static_cast<std::size_t>(n - k)],
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
        }
        sine.coefficients[static_cast<std::size_t>(n)] = divide(
            sineSum, integer(n), builtins, mathematics, angles, assumptions);
        cosine.coefficients[static_cast<std::size_t>(n)] = negate(
            divide(cosineSum, integer(n), builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(sine);
    trimLeadingZeros(cosine);
    return std::pair{std::move(sine), std::move(cosine)};
}

[[nodiscard]] std::optional<std::pair<LaurentSeries, LaurentSeries>> expandSinhCoshSeries(
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!argument.exactZero && argument.minimumExponent < 0)
        return std::nullopt;
    if (ceiling < 0)
        return std::pair{zeroSeries(ceiling), zeroSeries(ceiling)};

    LaurentSeries sine = denseSeries(0, ceiling);
    LaurentSeries cosine = denseSeries(0, ceiling);
    const Expr a0 = coefficientAt(argument, 0);
    sine.coefficients[0] = simplify(
        call(BuiltinId::Sinh, {a0}, builtins),
        builtins, mathematics, angles, assumptions);
    cosine.coefficients[0] = simplify(
        call(BuiltinId::Cosh, {a0}, builtins),
        builtins, mathematics, angles, assumptions);

    for (std::int64_t n = 1; n <= ceiling; ++n) {
        Expr sineSum = integer(0);
        Expr cosineSum = integer(0);
        for (std::int64_t k = 1; k <= n; ++k) {
            const Expr ak = coefficientAt(argument, k);
            if (isZero(ak)) continue;
            Expr weighted = multiply(
                integer(k), ak, builtins, mathematics, angles, assumptions);
            sineSum = add(
                sineSum,
                multiply(weighted, cosine.coefficients[static_cast<std::size_t>(n - k)],
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
            cosineSum = add(
                cosineSum,
                multiply(weighted, sine.coefficients[static_cast<std::size_t>(n - k)],
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
        }
        sine.coefficients[static_cast<std::size_t>(n)] = divide(
            sineSum, integer(n), builtins, mathematics, angles, assumptions);
        cosine.coefficients[static_cast<std::size_t>(n)] = divide(
            cosineSum, integer(n), builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(sine);
    trimLeadingZeros(cosine);
    return std::pair{std::move(sine), std::move(cosine)};
}

[[nodiscard]] LaurentSeries differentiateFormalSeries(
    const LaurentSeries& source,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (source.exactZero)
        return zeroSeries(ceiling);

    const std::int64_t minimum = source.minimumExponent - 1;
    if (minimum > ceiling)
        return zeroSeries(ceiling);
    LaurentSeries result = denseSeries(minimum, ceiling);
    const std::int64_t maximumSourceExponent = std::min(source.ceiling, ceiling + 1);
    for (std::int64_t exponent = source.minimumExponent;
         exponent <= maximumSourceExponent; ++exponent) {
        if (exponent == 0) continue;
        const Expr coefficient = coefficientAt(source, exponent);
        if (isZero(coefficient)) continue;
        result.coefficients[static_cast<std::size_t>(exponent - 1 - minimum)] = multiply(
            integer(exponent), coefficient,
            builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> integrateFormalDerivative(
    const LaurentSeries& derivative,
    const Expr& constant,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (ceiling < 0)
        return zeroSeries(ceiling);

    if (!derivative.exactZero) {
        const Expr logarithmicCoefficient = coefficientAt(derivative, -1);
        if (!isZero(logarithmicCoefficient))
            return std::nullopt;
    }

    const std::int64_t minimum = derivative.exactZero
        ? 0
        : std::min<std::int64_t>(0, derivative.minimumExponent + 1);
    LaurentSeries result = denseSeries(minimum, ceiling);
    if (0 >= minimum && 0 <= ceiling)
        result.coefficients[static_cast<std::size_t>(-minimum)] = constant;

    if (!derivative.exactZero) {
        for (std::int64_t exponent = derivative.minimumExponent;
             exponent <= derivative.ceiling; ++exponent) {
            if (exponent == -1 || exponent + 1 > ceiling) continue;
            const Expr coefficient = coefficientAt(derivative, exponent);
            if (isZero(coefficient)) continue;
            const std::int64_t target = exponent + 1;
            if (target < minimum) continue;
            result.coefficients[static_cast<std::size_t>(target - minimum)] = divide(
                coefficient, integer(target),
                builtins, mathematics, angles, assumptions);
        }
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] bool principalIntegralAnalyticAt(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    // Ei/Ciはprincipal Logと同じ負実軸の分岐切断を持ち，0で特異である。
    // 正の実中心または非実と証明できる中心だけを正則近傍として扱う。
    return principalPowerAnalyticAt(value, builtins, mathematics, assumptions);
}

[[nodiscard]] std::optional<LaurentSeries> expandClassicalIntegralSeries(
    BuiltinId id,
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!argument.exactZero && argument.minimumExponent < 0)
        return std::nullopt;
    if (ceiling < 0)
        return zeroSeries(ceiling);

    const Expr a0 = coefficientAt(argument, 0);
    if ((id == BuiltinId::ExponentialIntegralEi || id == BuiltinId::CosineIntegralCi)
        && !principalIntegralAnalyticAt(a0, builtins, mathematics, assumptions))
        return std::nullopt;
    if (id == BuiltinId::LogarithmicIntegralLi
        && !logarithmicIntegralAnalyticAt(a0, builtins, mathematics, assumptions))
        return std::nullopt;

    const bool zeroAtOrigin =
        id == BuiltinId::Erf
        || id == BuiltinId::SineIntegralSi
        || id == BuiltinId::FresnelC
        || id == BuiltinId::FresnelS;
    const Expr constant = isZero(a0) && zeroAtOrigin
        ? integer(0)
        : isZero(a0) && id == BuiltinId::Erfc
            ? integer(1)
            : simplify(call(id, {a0}, builtins),
                builtins, mathematics, angles, assumptions);
    if (argument.exactZero) {
        LaurentSeries result = denseSeries(0, ceiling);
        if (!result.coefficients.empty()) result.coefficients[0] = constant;
        trimLeadingZeros(result);
        return result;
    }

    const std::int64_t derivativeCeiling = std::max<std::int64_t>(0, ceiling - 1);
    LaurentSeries argumentDerivative = differentiateFormalSeries(
        argument, derivativeCeiling,
        builtins, mathematics, angles, assumptions);

    std::optional<LaurentSeries> kernel;
    switch (id) {
    case BuiltinId::Erf:
    case BuiltinId::Erfc: {
        LaurentSeries square = multiplySeries(
            argument, argument, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        for (Expr& coefficient : square.coefficients)
            coefficient = negate(coefficient, builtins, mathematics, angles, assumptions);
        auto exponential = expandExpSeries(
            square, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        if (!exponential) return std::nullopt;
        Expr scale = divide(
            integer(2), simplify(call(BuiltinId::Sqrt, {pi(mathematics)}, builtins),
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions);
        if (id == BuiltinId::Erfc)
            scale = negate(scale, builtins, mathematics, angles, assumptions);
        for (Expr& coefficient : exponential->coefficients)
            coefficient = multiply(
                scale, coefficient, builtins, mathematics, angles, assumptions);
        kernel = std::move(*exponential);
        break;
    }
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS: {
        LaurentSeries phase = multiplySeries(
            argument, argument, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        const Expr scale = divide(
            pi(mathematics), integer(2),
            builtins, mathematics, angles, assumptions);
        for (Expr& coefficient : phase.coefficients)
            coefficient = multiply(
                scale, coefficient, builtins, mathematics, angles, assumptions);
        auto sinCos = expandSinCosSeries(
            phase, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        if (!sinCos) return std::nullopt;
        kernel = id == BuiltinId::FresnelC
            ? std::move(sinCos->second)
            : std::move(sinCos->first);
        break;
    }
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi: {
        const std::int64_t quotientWorkCeiling = id == BuiltinId::SineIntegralSi
                && argument.minimumExponent > 0
            ? derivativeCeiling + argument.minimumExponent
            : derivativeCeiling;
        if (quotientWorkCeiling > kMaximumInternalExponent)
            return std::nullopt;
        auto sinCos = expandSinCosSeries(
            argument, quotientWorkCeiling,
            builtins, mathematics, angles, assumptions);
        if (!sinCos) return std::nullopt;
        auto inverse = inverseSeries(
            argument, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        if (!inverse) return std::nullopt;
        const LaurentSeries& numerator = id == BuiltinId::SineIntegralSi
            ? sinCos->first : sinCos->second;
        kernel = multiplySeries(
            numerator, *inverse, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        break;
    }
    case BuiltinId::ExponentialIntegralEi: {
        auto exponential = expandExpSeries(
            argument, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        auto inverse = inverseSeries(
            argument, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        if (!exponential || !inverse) return std::nullopt;
        kernel = multiplySeries(
            *exponential, *inverse, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        break;
    }
    case BuiltinId::LogarithmicIntegralLi: {
        auto logarithm = expandLogSeries(
            argument, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        if (!logarithm) return std::nullopt;
        auto inverse = inverseSeries(
            *logarithm, derivativeCeiling,
            builtins, mathematics, angles, assumptions, true);
        if (!inverse) return std::nullopt;
        kernel = std::move(*inverse);
        break;
    }
    default:
        return std::nullopt;
    }

    LaurentSeries derivative = multiplySeries(
        *kernel, argumentDerivative, derivativeCeiling,
        builtins, mathematics, angles, assumptions);
    return integrateFormalDerivative(
        derivative, constant, ceiling,
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<LaurentSeries> expandRationalPowerSeries(
    const LaurentSeries& base,
    const Expr& exponent,
    std::int64_t ceiling,
    BuiltinId resultHead,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    bool centerAnalyticKnown);

[[nodiscard]] std::optional<LaurentSeries> expandInverseTrigSeries(
    BuiltinId id,
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!argument.exactZero && argument.minimumExponent < 0)
        return std::nullopt;
    if (ceiling < 0)
        return zeroSeries(ceiling);

    const Expr a0 = coefficientAt(argument, 0);
    if (!inverseTrigAnalyticAt(id, a0, builtins, mathematics, angles, assumptions))
        return std::nullopt;

    const Expr constant = simplify(
        call(id, {a0}, builtins), builtins, mathematics, angles, assumptions);
    if (argument.exactZero) {
        LaurentSeries result = denseSeries(0, ceiling);
        if (!result.coefficients.empty()) result.coefficients[0] = constant;
        trimLeadingZeros(result);
        return result;
    }

    const std::int64_t derivativeCeiling = std::max<std::int64_t>(0, ceiling - 1);
    LaurentSeries square = multiplySeries(
        argument, argument, derivativeCeiling,
        builtins, mathematics, angles, assumptions);
    LaurentSeries one = denseSeries(0, derivativeCeiling);
    if (one.coefficients.empty())
        return std::nullopt;
    one.coefficients[0] = integer(1);

    std::optional<LaurentSeries> kernel;
    if (id == BuiltinId::Asin || id == BuiltinId::Acos) {
        LaurentSeries oneMinusSquare = addSeries(
            one, square, derivativeCeiling,
            builtins, mathematics, angles, assumptions, true);
        kernel = expandRationalPowerSeries(
            oneMinusSquare, rational(-1, 2), derivativeCeiling, BuiltinId::Power,
            builtins, mathematics, angles, assumptions, true);
    }
    else if (id == BuiltinId::Atan) {
        LaurentSeries onePlusSquare = addSeries(
            one, square, derivativeCeiling,
            builtins, mathematics, angles, assumptions);
        kernel = inverseSeries(
            onePlusSquare, derivativeCeiling,
            builtins, mathematics, angles, assumptions, true);
    }
    if (!kernel)
        return std::nullopt;

    Expr scale = inverseTrigScale(
        builtins, mathematics, angles, assumptions);
    if (id == BuiltinId::Acos)
        scale = negate(scale, builtins, mathematics, angles, assumptions);
    for (Expr& coefficient : kernel->coefficients)
        coefficient = multiply(
            scale, coefficient, builtins, mathematics, angles, assumptions);

    LaurentSeries argumentDerivative = differentiateFormalSeries(
        argument, derivativeCeiling,
        builtins, mathematics, angles, assumptions);
    LaurentSeries derivative = multiplySeries(
        *kernel, argumentDerivative, derivativeCeiling,
        builtins, mathematics, angles, assumptions);
    return integrateFormalDerivative(
        derivative, constant, ceiling,
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] LaurentSeries constantFormalSeries(
    const Expr& value,
    std::int64_t ceiling) {
    if (ceiling < 0)
        return zeroSeries(ceiling);
    LaurentSeries result = denseSeries(0, ceiling);
    if (!result.coefficients.empty())
        result.coefficients[0] = value;
    trimLeadingZeros(result);
    return result;
}

void scaleFormalSeries(
    LaurentSeries& series,
    const Expr& scale,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (series.exactZero || isOne(scale)) return;
    for (Expr& coefficient : series.coefficients)
        coefficient = multiply(
            scale, coefficient, builtins, mathematics, angles, assumptions);
    trimLeadingZeros(series);
}

[[nodiscard]] bool formalSeriesProvablyReal(
    const LaurentSeries& series,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AssumptionSet& assumptions) {
    if (series.exactZero) return true;
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    for (const Expr& coefficient : series.coefficients) {
        if (isZero(coefficient)) continue;
        if (!knowledge.facts(coefficient).isProvablyReal())
            return false;
    }
    return true;
}

[[nodiscard]] Expr gammaLogCoefficient(
    GammaSeriesBase base,
    std::int64_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (order == 1) {
        Expr coefficient = simplify(
            call(BuiltinId::Digamma, {integer(1)}, builtins),
            builtins, mathematics, angles, assumptions);
        if (base == GammaSeriesBase::Half) {
            Expr logTwo = simplify(
                call(BuiltinId::Log, {integer(2)}, builtins),
                builtins, mathematics, angles, assumptions);
            coefficient = subtract(
                coefficient,
                multiply(integer(2), logTwo,
                    builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
        }
        return coefficient;
    }

    Expr zeta = specialFunctionValue(
        BuiltinId::Zeta, integer(order),
        builtins, mathematics, angles, assumptions);
    Expr numerator = zeta;
    if (base == GammaSeriesBase::Half) {
        Expr powerOfTwo = simplify(
            call(BuiltinId::Power, {integer(2), integer(order)}, builtins),
            builtins, mathematics, angles, assumptions);
        Expr factor = subtract(
            powerOfTwo, integer(1), builtins, mathematics, angles, assumptions);
        numerator = multiply(
            factor, numerator, builtins, mathematics, angles, assumptions);
    }
    Expr coefficient = divide(
        numerator, integer(order), builtins, mathematics, angles, assumptions);
    if ((order & 1) != 0)
        coefficient = negate(
            coefficient, builtins, mathematics, angles, assumptions);
    return coefficient;
}

[[nodiscard]] std::optional<LaurentSeries> composeGammaLogIncrement(
    GammaSeriesBase base,
    const LaurentSeries& variation,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (variation.exactZero || ceiling < 1)
        return zeroSeries(ceiling);
    if (variation.minimumExponent <= 0)
        return std::nullopt;

    const std::int64_t maximumPower = ceiling / variation.minimumExponent;
    if (maximumPower > 256)
        return std::nullopt;

    LaurentSeries result = zeroSeries(ceiling);
    LaurentSeries power = variation;
    for (std::int64_t order = 1; order <= maximumPower; ++order) {
        LaurentSeries term = power;
        const Expr coefficient = gammaLogCoefficient(
            base, order, builtins, mathematics, angles, assumptions);
        scaleFormalSeries(
            term, coefficient, builtins, mathematics, angles, assumptions);
        result = addSeries(
            result, term, ceiling,
            builtins, mathematics, angles, assumptions);
        if (order != maximumPower)
            power = multiplySeries(
                power, variation, ceiling,
                builtins, mathematics, angles, assumptions);
    }
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandGammaFamilySeries(
    BuiltinId id,
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (id != BuiltinId::Gamma && id != BuiltinId::LogGamma)
        return std::nullopt;
    if (!argument.exactZero && argument.minimumExponent < 0)
        return std::nullopt;
    if (ceiling < 0)
        return zeroSeries(ceiling);

    const Expr center = coefficientAt(argument, 0);
    const auto plan = gammaSeriesPlan(center);
    if (!plan)
        return std::nullopt;
    if (id == BuiltinId::LogGamma
        && !formalSeriesProvablyReal(argument, builtins, mathematics, assumptions))
        return std::nullopt;

    LaurentSeries variation = addSeries(
        argument, constantFormalSeries(
            negate(center, builtins, mathematics, angles, assumptions), ceiling),
        ceiling, builtins, mathematics, angles, assumptions);
    auto deltaLog = composeGammaLogIncrement(
        plan->base, variation, ceiling,
        builtins, mathematics, angles, assumptions);
    if (!deltaLog)
        return std::nullopt;

    // 函数方程式Gamma(z+1)=z*Gamma(z)の対数微分だけを加える。
    // 各logの定数項は捨て，最後に目的中心のexact値を一度だけ戻す。
    const auto addLogFactor = [&](LaurentSeries factor, bool subtractFactor)
        -> bool {
        auto logarithm = expandLogSeries(
            factor, ceiling, builtins, mathematics, angles, assumptions);
        if (!logarithm) return false;
        if (!logarithm->exactZero && logarithm->minimumExponent <= 0
            && logarithm->ceiling >= 0) {
            const auto index = static_cast<std::size_t>(-logarithm->minimumExponent);
            if (index < logarithm->coefficients.size())
                logarithm->coefficients[index] = integer(0);
            trimLeadingZeros(*logarithm);
        }
        *deltaLog = addSeries(
            *deltaLog, *logarithm, ceiling,
            builtins, mathematics, angles, assumptions, subtractFactor);
        return true;
    };

    if (plan->shift > 0) {
        for (std::int64_t j = 1; j <= plan->shift; ++j) {
            LaurentSeries factor = addSeries(
                argument, constantFormalSeries(integer(-j), ceiling), ceiling,
                builtins, mathematics, angles, assumptions);
            if (!addLogFactor(std::move(factor), false))
                return std::nullopt;
        }
    }
    else if (plan->shift < 0) {
        for (std::int64_t j = 0; j < -plan->shift; ++j) {
            LaurentSeries factor = addSeries(
                argument, constantFormalSeries(integer(j), ceiling), ceiling,
                builtins, mathematics, angles, assumptions);
            // 負の半整数中心ではfactor<0なので，log[-factor]で同じ対数導函数を得る。
            for (Expr& coefficient : factor.coefficients)
                coefficient = negate(
                    coefficient, builtins, mathematics, angles, assumptions);
            if (!addLogFactor(std::move(factor), true))
                return std::nullopt;
        }
    }

    if (id == BuiltinId::LogGamma) {
        Expr constant = specialFunctionValue(
            BuiltinId::LogGamma, center,
            builtins, mathematics, angles, assumptions);
        if (deltaLog->exactZero)
            *deltaLog = constantFormalSeries(constant, ceiling);
        else if (deltaLog->minimumExponent > 0) {
            *deltaLog = addSeries(
                *deltaLog, constantFormalSeries(constant, ceiling), ceiling,
                builtins, mathematics, angles, assumptions);
        }
        else {
            const auto index = static_cast<std::size_t>(-deltaLog->minimumExponent);
            if (index < deltaLog->coefficients.size())
                deltaLog->coefficients[index] = constant;
        }
        trimLeadingZeros(*deltaLog);
        return deltaLog;
    }

    auto gamma = expandExpSeries(
        *deltaLog, ceiling,
        builtins, mathematics, angles, assumptions);
    if (!gamma)
        return std::nullopt;
    Expr constant = specialFunctionValue(
        BuiltinId::Gamma, center,
        builtins, mathematics, angles, assumptions);
    scaleFormalSeries(
        *gamma, constant, builtins, mathematics, angles, assumptions);
    return gamma;
}

[[nodiscard]] Expr gammaLogDerivativeCoefficient(
    GammaSeriesBase base,
    BuiltinId id,
    std::int64_t power,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const std::int64_t derivativeOrder = id == BuiltinId::Digamma ? 1 : 2;
    const std::int64_t logOrder = power + derivativeOrder;
    Expr coefficient = gammaLogCoefficient(
        base, logOrder, builtins, mathematics, angles, assumptions);
    coefficient = multiply(
        integer(logOrder), coefficient,
        builtins, mathematics, angles, assumptions);
    if (derivativeOrder == 2)
        coefficient = multiply(
            integer(logOrder - 1), coefficient,
            builtins, mathematics, angles, assumptions);
    return coefficient;
}

[[nodiscard]] std::optional<LaurentSeries> composeGammaLogDerivative(
    GammaSeriesBase base,
    BuiltinId id,
    const LaurentSeries& variation,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (id != BuiltinId::Digamma && id != BuiltinId::Trigamma)
        return std::nullopt;
    if (!variation.exactZero && variation.minimumExponent <= 0)
        return std::nullopt;

    LaurentSeries result = constantFormalSeries(
        gammaLogDerivativeCoefficient(
            base, id, 0, builtins, mathematics, angles, assumptions),
        ceiling);
    if (variation.exactZero || ceiling < variation.minimumExponent)
        return result;

    const std::int64_t maximumPower = ceiling / variation.minimumExponent;
    if (maximumPower > 256)
        return std::nullopt;

    LaurentSeries power = variation;
    for (std::int64_t exponent = 1; exponent <= maximumPower; ++exponent) {
        LaurentSeries term = power;
        scaleFormalSeries(
            term,
            gammaLogDerivativeCoefficient(
                base, id, exponent,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions);
        result = addSeries(
            result, term, ceiling,
            builtins, mathematics, angles, assumptions);
        if (exponent != maximumPower)
            power = multiplySeries(
                power, variation, ceiling,
                builtins, mathematics, angles, assumptions);
    }
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandPsiFamilySeries(
    BuiltinId id,
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (id != BuiltinId::Digamma && id != BuiltinId::Trigamma)
        return std::nullopt;
    if (!argument.exactZero && argument.minimumExponent < 0)
        return std::nullopt;
    if (ceiling < 0)
        return zeroSeries(ceiling);

    const Expr center = coefficientAt(argument, 0);
    const auto plan = gammaSeriesPlan(center);
    if (!plan)
        return std::nullopt;

    LaurentSeries variation = addSeries(
        argument, constantFormalSeries(
            negate(center, builtins, mathematics, angles, assumptions), ceiling),
        ceiling, builtins, mathematics, angles, assumptions);
    auto result = composeGammaLogDerivative(
        plan->base, id, variation, ceiling,
        builtins, mathematics, angles, assumptions);
    if (!result)
        return std::nullopt;

    // ψ(z+1)=ψ(z)+1/zとその微分を使い，基底中心から目的中心へ移送する。
    const auto addReciprocalFactor = [&](LaurentSeries factor, bool subtractFactor)
        -> bool {
        auto inverse = inverseSeries(
            factor, ceiling, builtins, mathematics, angles, assumptions);
        if (!inverse) return false;
        if (id == BuiltinId::Digamma) {
            *result = addSeries(
                *result, *inverse, ceiling,
                builtins, mathematics, angles, assumptions, subtractFactor);
            return true;
        }

        LaurentSeries square = multiplySeries(
            *inverse, *inverse, ceiling,
            builtins, mathematics, angles, assumptions);
        *result = addSeries(
            *result, square, ceiling,
            builtins, mathematics, angles, assumptions, !subtractFactor);
        return true;
    };

    if (plan->shift > 0) {
        for (std::int64_t j = 1; j <= plan->shift; ++j) {
            LaurentSeries factor = addSeries(
                argument, constantFormalSeries(integer(-j), ceiling), ceiling,
                builtins, mathematics, angles, assumptions);
            if (!addReciprocalFactor(std::move(factor), false))
                return std::nullopt;
        }
    }
    else if (plan->shift < 0) {
        for (std::int64_t j = 0; j < -plan->shift; ++j) {
            LaurentSeries factor = addSeries(
                argument, constantFormalSeries(integer(j), ceiling), ceiling,
                builtins, mathematics, angles, assumptions);
            if (!addReciprocalFactor(std::move(factor), true))
                return std::nullopt;
        }
    }

    trimLeadingZeros(*result);
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandPolylogOriginSeries(
    const Expr& order,
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (argument.exactZero)
        return zeroSeries(ceiling);
    if (argument.minimumExponent <= 0 || ceiling < argument.minimumExponent)
        return std::nullopt;

    const std::int64_t maximumPower = ceiling / argument.minimumExponent;
    if (maximumPower <= 0 || maximumPower > 256)
        return std::nullopt;

    LaurentSeries result = zeroSeries(ceiling);
    LaurentSeries power = argument;
    const Expr negativeOrder = negate(
        order, builtins, mathematics, angles, assumptions);
    for (std::int64_t n = 1; n <= maximumPower; ++n) {
        Expr coefficient = integer(1);
        if (n != 1)
            coefficient = simplify(
                call(BuiltinId::Power, {integer(n), negativeOrder}, builtins),
                builtins, mathematics, angles, assumptions);
        LaurentSeries term = power;
        scaleFormalSeries(
            term, coefficient, builtins, mathematics, angles, assumptions);
        result = addSeries(
            result, term, ceiling,
            builtins, mathematics, angles, assumptions);
        if (n != maximumPower)
            power = multiplySeries(
                power, argument, ceiling,
                builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<Expr> nonpositiveIntegerPolylogValue(
    std::int64_t order,
    const Expr& center,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (order > 0)
        return std::nullopt;

    const Expr oneMinusCenter = subtract(
        integer(1), center, builtins, mathematics, angles, assumptions);
    if (order == 0)
        return divide(
            center, oneMinusCenter,
            builtins, mathematics, angles, assumptions);

    constexpr std::int64_t kMaximumEulerianOrder = 64;
    const std::int64_t degree = -order;
    if (degree > kMaximumEulerianOrder)
        return std::nullopt;

    // Li_{-m}(z)=z*A_m(z)/(1-z)^(m+1) のEulerian多項式をexact整数係数で生成する。
    std::vector<numeric::BigInt> eulerian{numeric::BigInt{1}};
    for (std::int64_t n = 2; n <= degree; ++n) {
        std::vector<numeric::BigInt> next(
            static_cast<std::size_t>(n), numeric::BigInt{0});
        for (std::int64_t k = 0; k < n; ++k) {
            if (k < n - 1)
                next[static_cast<std::size_t>(k)] +=
                    numeric::BigInt{k + 1} * eulerian[static_cast<std::size_t>(k)];
            if (k > 0)
                next[static_cast<std::size_t>(k)] +=
                    numeric::BigInt{n - k} * eulerian[static_cast<std::size_t>(k - 1)];
        }
        eulerian = std::move(next);
    }

    Expr polynomial = integer(0);
    Expr centerPower = integer(1);
    for (std::size_t k = 0; k < eulerian.size(); ++k) {
        Expr term = centerPower;
        if (!(eulerian[k] == numeric::BigInt{1}))
            term = multiply(
                Expr{numeric::Number{eulerian[k]}}, term,
                builtins, mathematics, angles, assumptions);
        polynomial = add(
            polynomial, term,
            builtins, mathematics, angles, assumptions);
        if (k + 1 != eulerian.size())
            centerPower = multiply(
                centerPower, center,
                builtins, mathematics, angles, assumptions);
    }

    Expr numerator = multiply(
        center, polynomial,
        builtins, mathematics, angles, assumptions);
    Expr denominator = integer(1);
    for (std::int64_t n = 0; n <= degree; ++n)
        denominator = multiply(
            denominator, oneMinusCenter,
            builtins, mathematics, angles, assumptions);
    return divide(
        numerator, denominator,
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<Expr> integerOrderPolylogCenterValue(
    std::int64_t order,
    const Expr& center,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (order <= 0)
        return nonpositiveIntegerPolylogValue(
            order, center, builtins, mathematics, angles, assumptions);

    const std::array<Expr, 2> arguments{integer(order), center};
    return simplify(
        builtins::evaluateSpecialFunction(
            BuiltinId::Polylog, arguments,
            builtins, mathematics, angles),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<LaurentSeries> expandPolylogRegularCenterSeries(
    const Expr& order,
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto integerOrder = exactInt64(order);
    if (!integerOrder || *integerOrder <= 0 || *integerOrder > 64
        || argument.exactZero || argument.minimumExponent < 0)
        return std::nullopt;

    const Expr center = coefficientAt(argument, 0);
    if (isZero(center)
        || !provablyNonZero(center, builtins, mathematics, assumptions)
        || !polylogRegularCenterAnalyticAt(
            order, center, builtins, mathematics, assumptions))
        return std::nullopt;

    LaurentSeries variation = argument;
    if (variation.minimumExponent == 0 && !variation.coefficients.empty())
        variation.coefficients[0] = integer(0);
    trimLeadingZeros(variation);

    auto constant = integerOrderPolylogCenterValue(
        *integerOrder, center,
        builtins, mathematics, angles, assumptions);
    if (!constant)
        return std::nullopt;
    LaurentSeries result = constantFormalSeries(*constant, ceiling);
    if (variation.exactZero || ceiling < variation.minimumExponent)
        return result;
    if (variation.minimumExponent <= 0)
        return std::nullopt;

    const std::int64_t maximumPower = ceiling / variation.minimumExponent;
    if (maximumPower <= 0 || maximumPower > 256
        || maximumPower - *integerOrder > 64)
        return std::nullopt;

    // D^n Li_s(z)=z^-n sum_k s(n,k) Li_(s-k)(z) を使い，
    // signed Stirling係数を次数ごとに更新してTaylor係数を直接構成する。
    std::vector<numeric::BigInt> stirling(
        static_cast<std::size_t>(maximumPower + 1), numeric::BigInt{0});
    stirling[0] = numeric::BigInt{1};
    numeric::BigInt factorial{1};
    Expr centerPower = integer(1);
    LaurentSeries power = variation;

    for (std::int64_t n = 1; n <= maximumPower; ++n) {
        std::vector<numeric::BigInt> next(
            static_cast<std::size_t>(maximumPower + 1), numeric::BigInt{0});
        for (std::int64_t k = 1; k <= n; ++k)
            next[static_cast<std::size_t>(k)] =
                stirling[static_cast<std::size_t>(k - 1)]
                - numeric::BigInt{n - 1} * stirling[static_cast<std::size_t>(k)];
        stirling = std::move(next);
        factorial *= numeric::BigInt{n};
        centerPower = multiply(
            centerPower, center,
            builtins, mathematics, angles, assumptions);

        Expr derivativeNumerator = integer(0);
        for (std::int64_t k = 1; k <= n; ++k) {
            const numeric::BigInt& stirlingCoefficient =
                stirling[static_cast<std::size_t>(k)];
            if (stirlingCoefficient.isZero())
                continue;
            auto value = integerOrderPolylogCenterValue(
                *integerOrder - k, center,
                builtins, mathematics, angles, assumptions);
            if (!value)
                return std::nullopt;
            Expr term = *value;
            if (!(stirlingCoefficient == numeric::BigInt{1}))
                term = multiply(
                    Expr{numeric::Number{stirlingCoefficient}}, term,
                    builtins, mathematics, angles, assumptions);
            derivativeNumerator = add(
                derivativeNumerator, term,
                builtins, mathematics, angles, assumptions);
        }

        Expr denominator = multiply(
            centerPower, Expr{numeric::Number{factorial}},
            builtins, mathematics, angles, assumptions);
        Expr coefficient = divide(
            derivativeNumerator, denominator,
            builtins, mathematics, angles, assumptions);

        LaurentSeries term = power;
        scaleFormalSeries(
            term, coefficient, builtins, mathematics, angles, assumptions);
        result = addSeries(
            result, term, ceiling,
            builtins, mathematics, angles, assumptions);
        if (n != maximumPower)
            power = multiplySeries(
                power, variation, ceiling,
                builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandSeries(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions);

[[nodiscard]] std::optional<LaurentSeries> expandProduct(
    std::span<const Expr> factors,
    const expression::Symbol& variable,
    const Expr& center,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (factors.empty()) {
        LaurentSeries result = denseSeries(0, ceiling);
        if (!result.exactZero && !result.coefficients.empty())
            result.coefficients[0] = integer(1);
        return result;
    }

    std::vector<std::int64_t> valuations;
    valuations.reserve(factors.size());
    std::int64_t totalValuation = 0;
    for (const Expr& factor : factors) {
        const auto valuation = structuralValuation(factor, variable, center, builtins, mathematics, angles, assumptions);
        if (valuation.state == Valuation::State::Zero)
            return zeroSeries(ceiling);
        if (valuation.state != Valuation::State::Finite)
            return std::nullopt;
        valuations.push_back(valuation.exponent);
        totalValuation += valuation.exponent;
    }
    if (totalValuation > ceiling)
        return zeroSeries(ceiling);

    LaurentSeries result = denseSeries(0, ceiling - totalValuation);
    if (result.exactZero || result.coefficients.empty())
        return std::nullopt;
    result.minimumExponent = 0;
    result.ceiling = ceiling - totalValuation;
    result.coefficients.assign(
        static_cast<std::size_t>(result.ceiling + 1), integer(0));
    result.coefficients[0] = integer(1);

    std::int64_t accumulatedValuation = 0;
    for (std::size_t i = 0; i < factors.size(); ++i) {
        const std::int64_t remainingValuation = totalValuation - accumulatedValuation - valuations[i];
        const std::int64_t requiredCeiling = ceiling - accumulatedValuation - remainingValuation;
        auto child = expandSeries(
            factors[i], variable, center, requiredCeiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        result = multiplySeries(
            result, *child, ceiling - remainingValuation,
            builtins, mathematics, angles, assumptions);
        accumulatedValuation += valuations[i];
    }
    if (!result.exactZero)
        result.ceiling = ceiling;
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandPower(
    const Expr& baseExpression,
    std::int64_t exponent,
    const expression::Symbol& variable,
    const Expr& center,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (exponent == 0) {
        LaurentSeries result = denseSeries(0, ceiling);
        if (!result.exactZero && !result.coefficients.empty())
            result.coefficients[0] = integer(1);
        return result;
    }

    const auto baseValuation = structuralValuation(baseExpression, variable, center, builtins, mathematics, angles, assumptions);
    if (baseValuation.state != Valuation::State::Finite)
        return std::nullopt;
    const std::int64_t outputValuation = baseValuation.exponent * exponent;
    if (outputValuation > ceiling)
        return zeroSeries(ceiling);

    if (exponent < 0) {
        const std::uint64_t magnitude = static_cast<std::uint64_t>(-(exponent + 1)) + 1U;
        if (magnitude > 4096)
            return std::nullopt;
        const std::int64_t positiveValuation =
            baseValuation.exponent * static_cast<std::int64_t>(magnitude);
        const std::int64_t requiredPositiveCeiling = std::max(
            positiveValuation,
            ceiling + 2 * positiveValuation);
        if (requiredPositiveCeiling > kMaximumInternalExponent)
            return std::nullopt;

        auto positivePower = expandPower(
            baseExpression, static_cast<std::int64_t>(magnitude),
            variable, center, requiredPositiveCeiling,
            builtins, mathematics, angles, assumptions);
        if (!positivePower)
            return std::nullopt;
        const bool leadingNonZeroKnown = baseValuation.exponent == 0
            && seriesCenterValueProvablyNonZero(
                baseExpression, variable, center,
                builtins, mathematics, angles, assumptions);
        return inverseSeries(
            *positivePower, ceiling,
            builtins, mathematics, angles, assumptions, leadingNonZeroKnown);
    }

    const std::uint64_t magnitude = static_cast<std::uint64_t>(exponent);
    if (magnitude > 4096)
        return std::nullopt;
    std::vector<Expr> factors(static_cast<std::size_t>(magnitude), baseExpression);
    return expandProduct(
        factors, variable, center, ceiling,
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<LaurentSeries> expandRationalPowerSeries(
    const LaurentSeries& base,
    const Expr& exponent,
    std::int64_t ceiling,
    BuiltinId resultHead,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    bool centerAnalyticKnown = false) {
    if (base.exactZero || base.minimumExponent != 0 || ceiling < 0)
        return std::nullopt;

    const Expr a0 = coefficientAt(base, 0);
    if (!centerAnalyticKnown
        && !principalPowerAnalyticAt(a0, builtins, mathematics, assumptions))
        return std::nullopt;

    LaurentSeries result = denseSeries(0, ceiling);
    if (result.exactZero || result.coefficients.empty())
        return std::nullopt;

    if (isOne(a0)) {
        result.coefficients[0] = integer(1);
    }
    else if (resultHead == BuiltinId::Sqrt) {
        result.coefficients[0] = simplify(
            call(BuiltinId::Sqrt, {a0}, builtins),
            builtins, mathematics, angles, assumptions);
    }
    else {
        result.coefficients[0] = simplify(
            call(BuiltinId::Power, {a0, exponent}, builtins),
            builtins, mathematics, angles, assumptions);
    }

    const Expr exponentPlusOne = add(
        exponent, integer(1), builtins, mathematics, angles, assumptions);
    for (std::int64_t n = 1; n <= ceiling; ++n) {
        Expr sum = integer(0);
        for (std::int64_t k = 1; k <= n; ++k) {
            const Expr ak = coefficientAt(base, k);
            if (isZero(ak)) continue;

            Expr weight = subtract(
                multiply(exponentPlusOne, integer(k),
                    builtins, mathematics, angles, assumptions),
                integer(n), builtins, mathematics, angles, assumptions);
            Expr term = multiply(
                multiply(weight, ak, builtins, mathematics, angles, assumptions),
                result.coefficients[static_cast<std::size_t>(n - k)],
                builtins, mathematics, angles, assumptions);
            sum = add(sum, term, builtins, mathematics, angles, assumptions);
        }

        Expr denominator = multiply(
            integer(n), a0, builtins, mathematics, angles, assumptions);
        result.coefficients[static_cast<std::size_t>(n)] = divide(
            sum, denominator, builtins, mathematics, angles, assumptions);
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::vector<numeric::BigInt> nextLambertDerivativePolynomial(
    const std::vector<numeric::BigInt>& previous,
    std::uint64_t n) {
    std::vector<numeric::BigInt> next(previous.size() + 1, numeric::BigInt{0});
    const numeric::BigInt order = numeric::BigInt::fromUnsigned(n);
    const numeric::BigInt constantFactor = numeric::BigInt{1} - numeric::BigInt{3} * order;
    for (std::size_t i = 0; i < previous.size(); ++i) {
        next[i] += constantFactor * previous[i];
        next[i + 1] -= order * previous[i];
        if (i == 0) continue;
        const numeric::BigInt derivativeFactor = numeric::BigInt::fromUnsigned(i);
        const numeric::BigInt derivativeTerm = derivativeFactor * previous[i];
        next[i - 1] += derivativeTerm;
        next[i] += derivativeTerm;
    }
    return next;
}

[[nodiscard]] Expr evaluateIntegerPolynomial(
    const std::vector<numeric::BigInt>& coefficients,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    Expr result = integer(0);
    for (auto iterator = coefficients.rbegin(); iterator != coefficients.rend(); ++iterator) {
        result = add(
            multiply(result, value, builtins, mathematics, angles, assumptions),
            Expr{numeric::Number{*iterator}},
            builtins, mathematics, angles, assumptions);
    }
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandLambertWSeries(
    std::span<const Expr> originalArguments,
    const LaurentSeries& argument,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto principal = lambertWPrincipalBranch(originalArguments);
    if (!principal || (!argument.exactZero && argument.minimumExponent < 0))
        return std::nullopt;
    if (ceiling < 0)
        return zeroSeries(ceiling);

    const Expr a0 = coefficientAt(argument, 0);
    if (!lambertWAnalyticAt(
            *principal, a0, builtins, mathematics, angles, assumptions))
        return std::nullopt;

    std::vector<Expr> centerArguments;
    if (originalArguments.size() == 2)
        centerArguments.push_back(originalArguments[0]);
    centerArguments.push_back(a0);
    const Expr constant = simplify(
        builtins::evaluateSpecialFunction(
            BuiltinId::LambertW, centerArguments, builtins, mathematics, angles),
        builtins, mathematics, angles, assumptions);

    LaurentSeries result = denseSeries(0, ceiling);
    if (result.exactZero || result.coefficients.empty())
        return std::nullopt;
    result.coefficients[0] = constant;

    LaurentSeries delta = argument;
    if (!delta.exactZero && delta.minimumExponent <= 0 && delta.ceiling >= 0)
        delta.coefficients[static_cast<std::size_t>(-delta.minimumExponent)] = integer(0);
    trimLeadingZeros(delta);
    if (delta.exactZero) {
        trimLeadingZeros(result);
        return result;
    }

    // DLMF 4.13.4_1, 4.13.4_2の導函数多項式を生成し，局所Taylor係数をTPSA合成する。
    std::vector<numeric::BigInt> derivativePolynomial{numeric::BigInt{1}};
    LaurentSeries deltaPower = delta;
    for (std::int64_t n = 1; n <= ceiling && !deltaPower.exactZero; ++n) {
        if (deltaPower.minimumExponent > ceiling)
            break;

        const Expr polynomial = evaluateIntegerPolynomial(
            derivativePolynomial, constant,
            builtins, mathematics, angles, assumptions);
        Expr exponential = integer(1);
        if (!isZero(constant)) {
            const Expr exponent = negate(
                multiply(integer(n), constant, builtins, mathematics, angles, assumptions),
                builtins, mathematics, angles, assumptions);
            exponential = simplify(
                call(BuiltinId::Exp, {exponent}, builtins),
                builtins, mathematics, angles, assumptions);
        }
        const Expr onePlusW = add(
            integer(1), constant, builtins, mathematics, angles, assumptions);
        const Expr denominatorPower = simplify(
            call(BuiltinId::Power,
                {onePlusW, integer(2 * n - 1)}, builtins),
            builtins, mathematics, angles, assumptions);
        const Expr factorial = Expr{numeric::Number{
            numeric::factorial(static_cast<std::uint64_t>(n))}};
        const Expr taylorCoefficient = divide(
            multiply(exponential, polynomial, builtins, mathematics, angles, assumptions),
            multiply(factorial, denominatorPower,
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions);

        LaurentSeries term = deltaPower;
        for (Expr& coefficient : term.coefficients)
            coefficient = multiply(
                taylorCoefficient, coefficient,
                builtins, mathematics, angles, assumptions);
        result = addSeries(
            result, term, ceiling,
            builtins, mathematics, angles, assumptions);

        if (n < ceiling) {
            derivativePolynomial = nextLambertDerivativePolynomial(
                derivativePolynomial, static_cast<std::uint64_t>(n));
            deltaPower = multiplySeries(
                deltaPower, delta, ceiling,
                builtins, mathematics, angles, assumptions);
        }
    }
    trimLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<LaurentSeries> expandSeries(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::int64_t ceiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (ceiling > kMaximumInternalExponent || ceiling < -kMaximumInternalExponent)
        return std::nullopt;

    if (!containsSymbol(expression, variable)) {
        if (isZero(expression))
            return zeroSeries(ceiling);
        if (0 > ceiling)
            return zeroSeries(ceiling);
        LaurentSeries result = denseSeries(0, ceiling);
        if (result.exactZero || result.coefficients.empty()) return std::nullopt;
        result.coefficients[0] = expression;
        return result;
    }

    if (expression.isSymbol() && expression.asSymbol() == variable) {
        const std::int64_t minimum = isZero(center) ? 1 : 0;
        if (minimum > ceiling)
            return zeroSeries(ceiling);
        LaurentSeries result = denseSeries(minimum, ceiling);
        if (minimum == 0) {
            result.coefficients[0] = center;
            if (ceiling >= 1)
                result.coefficients[1] = integer(1);
        }
        else {
            result.coefficients[0] = integer(1);
        }
        trimLeadingZeros(result);
        return result;
    }

    if (!expression.isCall())
        return std::nullopt;
    const auto id = builtinId(expression, builtins);
    if (!id) return std::nullopt;
    const auto& arguments = expression.asCall().arguments;

    switch (*id) {
    case BuiltinId::Add: {
        LaurentSeries result = zeroSeries(ceiling);
        for (const Expr& argument : arguments) {
            const auto child = expandSeries(
                argument, variable, center, ceiling,
                builtins, mathematics, angles, assumptions);
            if (!child) return std::nullopt;
            result = addSeries(
                result, *child, ceiling,
                builtins, mathematics, angles, assumptions);
        }
        return result;
    }
    case BuiltinId::Subtract: {
        if (arguments.size() != 2) return std::nullopt;
        const auto lhs = expandSeries(
            arguments[0], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        const auto rhs = expandSeries(
            arguments[1], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!lhs || !rhs) return std::nullopt;
        return addSeries(
            *lhs, *rhs, ceiling,
            builtins, mathematics, angles, assumptions, true);
    }
    case BuiltinId::Negate: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandSeries(
            arguments[0], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        for (Expr& coefficient : child->coefficients)
            coefficient = negate(coefficient, builtins, mathematics, angles, assumptions);
        return child;
    }
    case BuiltinId::Multiply:
        return expandProduct(
            arguments, variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
    case BuiltinId::Divide: {
        if (arguments.size() != 2) return std::nullopt;
        const auto numeratorValuation = structuralValuation(arguments[0], variable, center, builtins, mathematics, angles, assumptions);
        const auto denominatorValuation = structuralValuation(arguments[1], variable, center, builtins, mathematics, angles, assumptions);
        if (numeratorValuation.state == Valuation::State::Zero)
            return zeroSeries(ceiling);
        if (numeratorValuation.state != Valuation::State::Finite
            || denominatorValuation.state != Valuation::State::Finite)
            return std::nullopt;

        const std::int64_t quotientValuation =
            numeratorValuation.exponent - denominatorValuation.exponent;
        if (quotientValuation > ceiling)
            return zeroSeries(ceiling);

        const std::int64_t numeratorCeiling = ceiling + denominatorValuation.exponent;
        const std::int64_t inverseCeiling = ceiling - numeratorValuation.exponent;
        const std::int64_t denominatorCeiling = std::max(
            denominatorValuation.exponent,
            inverseCeiling + 2 * denominatorValuation.exponent);
        if (numeratorCeiling > kMaximumInternalExponent
            || denominatorCeiling > kMaximumInternalExponent)
            return std::nullopt;

        auto numerator = expandSeries(
            arguments[0], variable, center, numeratorCeiling,
            builtins, mathematics, angles, assumptions);
        auto denominator = expandSeries(
            arguments[1], variable, center, denominatorCeiling,
            builtins, mathematics, angles, assumptions);
        if (!numerator || !denominator) return std::nullopt;
        const bool denominatorLeadingNonZeroKnown = denominatorValuation.exponent == 0
            && seriesCenterValueProvablyNonZero(
                arguments[1], variable, center,
                builtins, mathematics, angles, assumptions);
        auto inverse = inverseSeries(
            *denominator, inverseCeiling,
            builtins, mathematics, angles, assumptions,
            denominatorLeadingNonZeroKnown);
        if (!inverse) return std::nullopt;
        return multiplySeries(
            *numerator, *inverse, ceiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Power: {
        if (arguments.size() != 2) return std::nullopt;
        const auto exponent = exactInt64(arguments[1]);
        if (exponent)
            return expandPower(
                arguments[0], *exponent, variable, center, ceiling,
                builtins, mathematics, angles, assumptions);

        if (!exactRational(arguments[1])) return std::nullopt;
        auto base = expandSeries(
            arguments[0], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!base) return std::nullopt;
        return expandRationalPowerSeries(
            *base, arguments[1], ceiling, BuiltinId::Power,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Sqrt: {
        if (arguments.size() != 1) return std::nullopt;
        auto base = expandSeries(
            arguments[0], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!base) return std::nullopt;
        return expandRationalPowerSeries(
            *base, rational(1, 2), ceiling, BuiltinId::Sqrt,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc:
    case BuiltinId::Tanh:
    case BuiltinId::Coth:
    case BuiltinId::Sech:
    case BuiltinId::Csch:
    case BuiltinId::Expm1:
    case BuiltinId::Log1p:
    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc:
    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc:
    case BuiltinId::Log2:
    case BuiltinId::Log10: {
        auto rewritten = lowCostSeriesRewrite(
            *id, arguments, builtins, mathematics, angles, assumptions);
        if (!rewritten) return std::nullopt;
        return expandSeries(
            *rewritten, variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Exp:
    case BuiltinId::Log:
    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Sinh:
    case BuiltinId::Cosh: {
        if (*id == BuiltinId::Log && arguments.size() == 2) {
            auto rewritten = lowCostSeriesRewrite(
                *id, arguments, builtins, mathematics, angles, assumptions);
            if (!rewritten) return std::nullopt;
            return expandSeries(
                *rewritten, variable, center, ceiling,
                builtins, mathematics, angles, assumptions);
        }
        if (arguments.size() != 1) return std::nullopt;
        Expr analyticArgument = (*id == BuiltinId::Sin || *id == BuiltinId::Cos)
            ? radianTrigArgument(arguments[0], builtins, mathematics, angles, assumptions)
            : arguments[0];
        auto child = expandSeries(
            analyticArgument, variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;

        switch (*id) {
        case BuiltinId::Exp:
            return expandExpSeries(
                *child, ceiling, builtins, mathematics, angles, assumptions);
        case BuiltinId::Log:
            return expandLogSeries(
                *child, ceiling, builtins, mathematics, angles, assumptions);
        case BuiltinId::Sin:
        case BuiltinId::Cos: {
            auto pair = expandSinCosSeries(
                *child, ceiling, builtins, mathematics, angles, assumptions);
            if (!pair) return std::nullopt;
            return *id == BuiltinId::Sin
                ? std::move(pair->first) : std::move(pair->second);
        }
        case BuiltinId::Sinh:
        case BuiltinId::Cosh: {
            auto pair = expandSinhCoshSeries(
                *child, ceiling, builtins, mathematics, angles, assumptions);
            if (!pair) return std::nullopt;
            return *id == BuiltinId::Sinh
                ? std::move(pair->first) : std::move(pair->second);
        }
        default:
            return std::nullopt;
        }
    }
    case BuiltinId::Polylog: {
        if (arguments.size() != 2 || containsSymbol(arguments[0], variable))
            return std::nullopt;
        auto child = expandSeries(
            arguments[1], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        if (child->exactZero || child->minimumExponent > 0)
            return expandPolylogOriginSeries(
                arguments[0], *child, ceiling,
                builtins, mathematics, angles, assumptions);
        return expandPolylogRegularCenterSeries(
            arguments[0], *child, ceiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Gamma:
    case BuiltinId::LogGamma: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandSeries(
            arguments[0], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        return expandGammaFamilySeries(
            *id, *child, ceiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Digamma:
    case BuiltinId::Trigamma: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandSeries(
            arguments[0], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        return expandPsiFamilySeries(
            *id, *child, ceiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::LambertW: {
        const auto principal = lambertWPrincipalBranch(arguments);
        if (!principal) return std::nullopt;
        auto child = expandSeries(
            arguments.back(), variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        return expandLambertWSeries(
            arguments, *child, ceiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Asin:
    case BuiltinId::Acos:
    case BuiltinId::Atan: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandSeries(
            arguments[0], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        return expandInverseTrigSeries(
            *id, *child, ceiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandSeries(
            arguments[0], variable, center, ceiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        return expandClassicalIntegralSeries(
            *id, *child, ceiling,
            builtins, mathematics, angles, assumptions);
    }
    default:
        return std::nullopt;
    }
}


struct PuiseuxSeries final {
    std::uint32_t denominator = 1;
    LaurentSeries ticks;
};

[[nodiscard]] std::optional<std::uint32_t> commonPuiseuxDenominator(
    std::uint32_t lhs,
    std::uint32_t rhs) {
    const std::uint64_t gcd = std::gcd(lhs, rhs);
    const std::uint64_t lcm = (static_cast<std::uint64_t>(lhs) / gcd) * rhs;
    if (lcm == 0 || lcm > 64)
        return std::nullopt;
    return static_cast<std::uint32_t>(lcm);
}

[[nodiscard]] std::optional<PuiseuxSeries> rescalePuiseux(
    const PuiseuxSeries& source,
    std::uint32_t denominator,
    std::int64_t ceiling) {
    if (denominator == 0 || denominator % source.denominator != 0)
        return std::nullopt;
    const std::int64_t scale = static_cast<std::int64_t>(denominator / source.denominator);
    if (source.ticks.exactZero)
        return PuiseuxSeries{denominator, zeroSeries(ceiling)};
    if (source.ticks.minimumExponent != 0
        && std::abs(source.ticks.minimumExponent) > kMaximumInternalExponent / scale)
        return std::nullopt;
    const std::int64_t minimum = source.ticks.minimumExponent * scale;
    if (minimum > ceiling)
        return PuiseuxSeries{denominator, zeroSeries(ceiling)};
    LaurentSeries result = denseSeries(minimum, ceiling);
    for (std::int64_t exponent = source.ticks.minimumExponent;
         exponent <= source.ticks.ceiling; ++exponent) {
        const Expr coefficient = coefficientAt(source.ticks, exponent);
        if (isZero(coefficient)) continue;
        const std::int64_t target = exponent * scale;
        if (target < minimum || target > ceiling) continue;
        result.coefficients[static_cast<std::size_t>(target - minimum)] = coefficient;
    }
    trimLeadingZeros(result);
    return PuiseuxSeries{denominator, std::move(result)};
}

[[nodiscard]] std::optional<PuiseuxSeries> addPuiseux(
    const PuiseuxSeries& lhs,
    const PuiseuxSeries& rhs,
    std::int64_t integerCeiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    bool subtractRight = false) {
    const auto denominator = commonPuiseuxDenominator(lhs.denominator, rhs.denominator);
    if (!denominator) return std::nullopt;
    const std::int64_t ceiling = integerCeiling * static_cast<std::int64_t>(*denominator);
    if (ceiling > kMaximumInternalExponent) return std::nullopt;
    auto left = rescalePuiseux(lhs, *denominator, ceiling);
    auto right = rescalePuiseux(rhs, *denominator, ceiling);
    if (!left || !right) return std::nullopt;
    return PuiseuxSeries{*denominator, addSeries(
        left->ticks, right->ticks, ceiling,
        builtins, mathematics, angles, assumptions, subtractRight)};
}

[[nodiscard]] std::optional<PuiseuxSeries> multiplyPuiseux(
    const PuiseuxSeries& lhs,
    const PuiseuxSeries& rhs,
    std::int64_t integerCeiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto denominator = commonPuiseuxDenominator(lhs.denominator, rhs.denominator);
    if (!denominator) return std::nullopt;
    const std::int64_t ceiling = integerCeiling * static_cast<std::int64_t>(*denominator);
    if (ceiling > kMaximumInternalExponent) return std::nullopt;
    auto left = rescalePuiseux(lhs, *denominator, ceiling);
    auto right = rescalePuiseux(rhs, *denominator, ceiling);
    if (!left || !right) return std::nullopt;
    return PuiseuxSeries{*denominator, multiplySeries(
        left->ticks, right->ticks, ceiling,
        builtins, mathematics, angles, assumptions)};
}

[[nodiscard]] std::optional<PuiseuxSeries> inversePuiseux(
    const PuiseuxSeries& source,
    std::int64_t integerCeiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const std::int64_t ceiling = integerCeiling * static_cast<std::int64_t>(source.denominator);
    if (ceiling > kMaximumInternalExponent) return std::nullopt;
    auto result = inverseSeries(
        source.ticks, ceiling,
        builtins, mathematics, angles, assumptions);
    if (!result) return std::nullopt;
    return PuiseuxSeries{source.denominator, std::move(*result)};
}

[[nodiscard]] std::optional<PuiseuxSeries> expandPuiseuxSeries(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::int64_t integerCeiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions);

[[nodiscard]] std::optional<PuiseuxSeries> branchPointPowerPuiseux(
    const PuiseuxSeries& base,
    const Expr& exponentExpression,
    std::int64_t integerCeiling,
    BuiltinId resultHead,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto exponent = smallRational(exponentExpression);
    if (!exponent || exponent->denominator == 1 || exponent->numerator == 0
        || base.ticks.exactZero)
        return std::nullopt;

    // ordinary analytic baseの高重複零点をprincipal powerへ
    // 無条件に潰さない。simple zero/pole，または既にbranchを持つPuiseux baseだけを扱う。
    if (base.denominator == 1 && std::abs(base.ticks.minimumExponent) != 1)
        return std::nullopt;

    const Expr leading = coefficientAt(base.ticks, base.ticks.minimumExponent);
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.facts(leading).sign != mathematics::RealSign::Positive)
        return std::nullopt;

    const std::uint64_t rawDenominator =
        static_cast<std::uint64_t>(base.denominator) * exponent->denominator;
    if (rawDenominator > 64)
        return std::nullopt;
    const auto targetDenominator = static_cast<std::uint32_t>(rawDenominator);
    const std::int64_t targetCeiling =
        integerCeiling * static_cast<std::int64_t>(targetDenominator);
    if (targetCeiling > kMaximumInternalExponent)
        return std::nullopt;

    const std::int64_t shift =
        base.ticks.minimumExponent * exponent->numerator;
    const std::int64_t unitCeiling = std::max<std::int64_t>(
        0, (targetCeiling - shift) / static_cast<std::int64_t>(exponent->denominator));
    LaurentSeries unit = denseSeries(0, unitCeiling);
    if (unit.exactZero) return std::nullopt;
    for (std::int64_t tick = 0; tick <= unitCeiling; ++tick)
        unit.coefficients[static_cast<std::size_t>(tick)] =
            coefficientAt(base.ticks, base.ticks.minimumExponent + tick);
    trimLeadingZeros(unit);
    if (unit.exactZero || unit.minimumExponent != 0)
        return std::nullopt;

    auto poweredUnit = expandRationalPowerSeries(
        unit, exponentExpression, unitCeiling, resultHead,
        builtins, mathematics, angles, assumptions);
    if (!poweredUnit) return std::nullopt;

    PuiseuxSeries unitSeries{base.denominator, std::move(*poweredUnit)};
    auto scaled = rescalePuiseux(unitSeries, targetDenominator, targetCeiling - shift);
    if (!scaled) return std::nullopt;
    if (scaled->ticks.exactZero)
        return scaled;
    scaled->ticks.minimumExponent += shift;
    scaled->ticks.ceiling += shift;
    if (scaled->ticks.minimumExponent > targetCeiling)
        scaled->ticks = zeroSeries(targetCeiling);
    else if (scaled->ticks.ceiling > targetCeiling) {
        const auto newSize = static_cast<std::size_t>(
            targetCeiling - scaled->ticks.minimumExponent + 1);
        if (newSize < scaled->ticks.coefficients.size())
            scaled->ticks.coefficients.erase(
                scaled->ticks.coefficients.begin() + static_cast<std::ptrdiff_t>(newSize),
                scaled->ticks.coefficients.end());
        scaled->ticks.ceiling = targetCeiling;
    }
    return scaled;
}

[[nodiscard]] std::optional<PuiseuxSeries> rationalPowerPuiseux(
    const PuiseuxSeries& base,
    const Expr& exponent,
    std::int64_t integerCeiling,
    BuiltinId resultHead,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const std::int64_t tickCeiling =
        integerCeiling * static_cast<std::int64_t>(base.denominator);
    if (tickCeiling > kMaximumInternalExponent)
        return std::nullopt;
    if (!base.ticks.exactZero && base.ticks.minimumExponent == 0) {
        auto result = expandRationalPowerSeries(
            base.ticks, exponent, tickCeiling, resultHead,
            builtins, mathematics, angles, assumptions);
        if (!result) return std::nullopt;
        return PuiseuxSeries{base.denominator, std::move(*result)};
    }
    return branchPointPowerPuiseux(
        base, exponent, integerCeiling, resultHead,
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<PuiseuxSeries> expandLambertWBranchPointPuiseux(
    std::span<const Expr> arguments,
    const expression::Symbol& variable,
    const Expr& center,
    std::int64_t integerCeiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto branch = lambertWBranchIndex(arguments);
    if (!branch || (*branch != 0 && *branch != -1))
        return std::nullopt;

    const Expr& argument = arguments.back();
    const Expr centerArgument = simplify(
        substituteSymbol(argument, variable, center),
        builtins, mathematics, angles, assumptions);
    if (!isLambertWBranchPoint(
            centerArgument, builtins, mathematics, angles, assumptions))
        return std::nullopt;

    // DLMF 4.13.9_1の局所変数s=sqrt[E z+1]を既存Puiseux演算へ接続する。
    Expr localBase = add(
        multiply(e(mathematics), argument, builtins, mathematics, angles, assumptions),
        integer(1), builtins, mathematics, angles, assumptions);
    Expr sExpression = call(BuiltinId::Sqrt, {localBase}, builtins);
    if (*branch == -1)
        sExpression = negate(
            sExpression, builtins, mathematics, angles, assumptions);

    auto local = expandPuiseuxSeries(
        sExpression, variable, center, integerCeiling + 2,
        builtins, mathematics, angles, assumptions);
    if (!local || local->ticks.exactZero)
        return std::nullopt;

    const std::int64_t tickCeiling =
        integerCeiling * static_cast<std::int64_t>(local->denominator);
    if (tickCeiling > kMaximumInternalExponent || tickCeiling > 512)
        return std::nullopt;

    PuiseuxSeries result{1, denseSeries(0, integerCeiling)};
    if (result.ticks.exactZero || result.ticks.coefficients.empty())
        return std::nullopt;
    result.ticks.coefficients[0] = integer(-1);

    // d_n=q_n*(sqrt[2])^(n mod 2)として，有理係数q_nだけを漸化式で更新する。
    std::vector<numeric::Rational> rationalCoefficients(
        static_cast<std::size_t>(tickCeiling + 1), numeric::Rational{});
    rationalCoefficients[0] = numeric::Rational{numeric::BigInt{-1}};
    if (tickCeiling >= 1)
        rationalCoefficients[1] = numeric::Rational{numeric::BigInt{1}};
    const numeric::Rational two{numeric::BigInt{2}};

    const auto productFactor = [&](std::int64_t lhs, std::int64_t rhs) {
        numeric::Rational value =
            rationalCoefficients[static_cast<std::size_t>(lhs)]
            * rationalCoefficients[static_cast<std::size_t>(rhs)];
        if ((lhs & 1) != 0 && (rhs & 1) != 0)
            value *= two;
        return value;
    };

    PuiseuxSeries localPower = *local;
    for (std::int64_t termIndex = 1;
         termIndex <= tickCeiling && !localPower.ticks.exactZero; ++termIndex) {
        if (termIndex >= 2) {
            const std::int64_t n = termIndex - 1;
            numeric::Rational firstSum;
            numeric::Rational secondSum;
            for (std::int64_t k = 1; k <= n - 1; ++k) {
                firstSum += productFactor(k, n - k);
                secondSum += productFactor(k + 1, n - k + 1);
            }

            numeric::Rational numerator =
                -two * rationalCoefficients[static_cast<std::size_t>(n)]
                + numeric::Rational{numeric::BigInt{n}, numeric::BigInt{2}} * firstSum
                - numeric::Rational{numeric::BigInt{n + 2}, numeric::BigInt{2}} * secondSum;
            const numeric::Rational denominator{numeric::BigInt{n + 2}};
            rationalCoefficients[static_cast<std::size_t>(termIndex)] = (n & 1) != 0
                ? numerator / denominator
                : numerator / (two * denominator);
        }

        Expr coefficient{numeric::Number{
            rationalCoefficients[static_cast<std::size_t>(termIndex)]}};
        if ((termIndex & 1) != 0) {
            const Expr sqrtTwo = simplify(
                call(BuiltinId::Sqrt, {integer(2)}, builtins),
                builtins, mathematics, angles, assumptions);
            coefficient = multiply(
                coefficient, sqrtTwo, builtins, mathematics, angles, assumptions);
        }
        PuiseuxSeries term = localPower;
        for (Expr& value : term.ticks.coefficients)
            value = multiply(
                coefficient, value, builtins, mathematics, angles, assumptions);
        auto next = addPuiseux(
            result, term, integerCeiling,
            builtins, mathematics, angles, assumptions);
        if (!next) return std::nullopt;
        result = std::move(*next);

        if (termIndex == tickCeiling)
            break;
        auto nextPower = multiplyPuiseux(
            localPower, *local, integerCeiling,
            builtins, mathematics, angles, assumptions);
        if (!nextPower) return std::nullopt;
        localPower = std::move(*nextPower);
    }
    return result;
}

[[nodiscard]] std::optional<PuiseuxSeries> expandPuiseuxSeries(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::int64_t integerCeiling,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (integerCeiling < 0 || integerCeiling > 1024)
        return std::nullopt;

    if (auto ordinary = expandSeries(
            expression, variable, center, integerCeiling,
            builtins, mathematics, angles, assumptions))
        return PuiseuxSeries{1, std::move(*ordinary)};

    if (!expression.isCall())
        return std::nullopt;
    const auto id = builtinId(expression, builtins);
    if (!id) return std::nullopt;
    const auto& arguments = expression.asCall().arguments;

    switch (*id) {
    case BuiltinId::Add: {
        PuiseuxSeries result{1, zeroSeries(integerCeiling)};
        for (const Expr& argument : arguments) {
            auto child = expandPuiseuxSeries(
                argument, variable, center, integerCeiling,
                builtins, mathematics, angles, assumptions);
            if (!child) return std::nullopt;
            auto next = addPuiseux(
                result, *child, integerCeiling,
                builtins, mathematics, angles, assumptions);
            if (!next) return std::nullopt;
            result = std::move(*next);
        }
        return result;
    }
    case BuiltinId::Subtract: {
        if (arguments.size() != 2) return std::nullopt;
        auto lhs = expandPuiseuxSeries(arguments[0], variable, center, integerCeiling,
            builtins, mathematics, angles, assumptions);
        auto rhs = expandPuiseuxSeries(arguments[1], variable, center, integerCeiling,
            builtins, mathematics, angles, assumptions);
        if (!lhs || !rhs) return std::nullopt;
        return addPuiseux(*lhs, *rhs, integerCeiling,
            builtins, mathematics, angles, assumptions, true);
    }
    case BuiltinId::Negate: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandPuiseuxSeries(arguments[0], variable, center, integerCeiling,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        for (Expr& coefficient : child->ticks.coefficients)
            coefficient = negate(coefficient, builtins, mathematics, angles, assumptions);
        return child;
    }
    case BuiltinId::Multiply: {
        PuiseuxSeries result{1, denseSeries(0, integerCeiling)};
        if (result.ticks.coefficients.empty()) return std::nullopt;
        result.ticks.coefficients[0] = integer(1);
        for (const Expr& factor : arguments) {
            auto child = expandPuiseuxSeries(factor, variable, center, integerCeiling,
                builtins, mathematics, angles, assumptions);
            if (!child) return std::nullopt;
            auto next = multiplyPuiseux(result, *child, integerCeiling,
                builtins, mathematics, angles, assumptions);
            if (!next) return std::nullopt;
            result = std::move(*next);
        }
        return result;
    }
    case BuiltinId::Divide: {
        if (arguments.size() != 2) return std::nullopt;
        auto lhs = expandPuiseuxSeries(arguments[0], variable, center, integerCeiling,
            builtins, mathematics, angles, assumptions);
        auto rhs = expandPuiseuxSeries(arguments[1], variable, center, integerCeiling + 2,
            builtins, mathematics, angles, assumptions);
        if (!lhs || !rhs) return std::nullopt;
        auto inv = inversePuiseux(*rhs, integerCeiling + 1,
            builtins, mathematics, angles, assumptions);
        if (!inv) return std::nullopt;
        return multiplyPuiseux(*lhs, *inv, integerCeiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Power: {
        if (arguments.size() != 2) return std::nullopt;
        if (const auto integerExponent = exactInt64(arguments[1])) {
            if (*integerExponent == 0)
                return expandPuiseuxSeries(integer(1), variable, center, integerCeiling,
                    builtins, mathematics, angles, assumptions);
            auto base = expandPuiseuxSeries(arguments[0], variable, center, integerCeiling + 2,
                builtins, mathematics, angles, assumptions);
            if (!base) return std::nullopt;
            if (*integerExponent < 0) {
                auto inv = inversePuiseux(*base, integerCeiling + 1,
                    builtins, mathematics, angles, assumptions);
                if (!inv) return std::nullopt;
                base = std::move(inv);
            }
            const std::uint64_t magnitude = static_cast<std::uint64_t>(
                *integerExponent < 0 ? -(*integerExponent + 1) + 1 : *integerExponent);
            PuiseuxSeries result{1, denseSeries(0, integerCeiling)};
            if (result.ticks.coefficients.empty()) return std::nullopt;
            result.ticks.coefficients[0] = integer(1);
            for (std::uint64_t i = 0; i < magnitude; ++i) {
                auto next = multiplyPuiseux(result, *base, integerCeiling,
                    builtins, mathematics, angles, assumptions);
                if (!next) return std::nullopt;
                result = std::move(*next);
            }
            return result;
        }
        if (!smallRational(arguments[1])) return std::nullopt;
        auto base = expandPuiseuxSeries(arguments[0], variable, center, integerCeiling + 4,
            builtins, mathematics, angles, assumptions);
        if (!base) return std::nullopt;
        return rationalPowerPuiseux(*base, arguments[1], integerCeiling,
            BuiltinId::Power, builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Sqrt: {
        if (arguments.size() != 1) return std::nullopt;
        auto base = expandPuiseuxSeries(arguments[0], variable, center, integerCeiling + 4,
            builtins, mathematics, angles, assumptions);
        if (!base) return std::nullopt;
        return rationalPowerPuiseux(*base, rational(1, 2), integerCeiling,
            BuiltinId::Sqrt, builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc:
    case BuiltinId::Tanh:
    case BuiltinId::Coth:
    case BuiltinId::Sech:
    case BuiltinId::Csch:
    case BuiltinId::Expm1:
    case BuiltinId::Log1p:
    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc:
    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc:
    case BuiltinId::Log2:
    case BuiltinId::Log10: {
        auto rewritten = lowCostSeriesRewrite(
            *id, arguments, builtins, mathematics, angles, assumptions);
        if (!rewritten) return std::nullopt;
        return expandPuiseuxSeries(
            *rewritten, variable, center, integerCeiling,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Exp:
    case BuiltinId::Log:
    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Sinh:
    case BuiltinId::Cosh: {
        if (*id == BuiltinId::Log && arguments.size() == 2) {
            auto rewritten = lowCostSeriesRewrite(
                *id, arguments, builtins, mathematics, angles, assumptions);
            if (!rewritten) return std::nullopt;
            return expandPuiseuxSeries(
                *rewritten, variable, center, integerCeiling,
                builtins, mathematics, angles, assumptions);
        }
        if (arguments.size() != 1) return std::nullopt;
        Expr analyticArgument = (*id == BuiltinId::Sin || *id == BuiltinId::Cos)
            ? radianTrigArgument(arguments[0], builtins, mathematics, angles, assumptions)
            : arguments[0];
        auto child = expandPuiseuxSeries(analyticArgument, variable, center, integerCeiling + 2,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        const std::int64_t tickCeiling =
            integerCeiling * static_cast<std::int64_t>(child->denominator);
        if (tickCeiling > kMaximumInternalExponent) return std::nullopt;
        std::optional<LaurentSeries> result;
        switch (*id) {
        case BuiltinId::Exp:
            result = expandExpSeries(child->ticks, tickCeiling,
                builtins, mathematics, angles, assumptions); break;
        case BuiltinId::Log:
            result = expandLogSeries(child->ticks, tickCeiling,
                builtins, mathematics, angles, assumptions); break;
        case BuiltinId::Sin:
        case BuiltinId::Cos: {
            auto pair = expandSinCosSeries(child->ticks, tickCeiling,
                builtins, mathematics, angles, assumptions);
            if (!pair) return std::nullopt;
            result = *id == BuiltinId::Sin ? std::move(pair->first) : std::move(pair->second);
            break;
        }
        case BuiltinId::Sinh:
        case BuiltinId::Cosh: {
            auto pair = expandSinhCoshSeries(child->ticks, tickCeiling,
                builtins, mathematics, angles, assumptions);
            if (!pair) return std::nullopt;
            result = *id == BuiltinId::Sinh ? std::move(pair->first) : std::move(pair->second);
            break;
        }
        default: break;
        }
        if (!result) return std::nullopt;
        return PuiseuxSeries{child->denominator, std::move(*result)};
    }
    case BuiltinId::Polylog: {
        if (arguments.size() != 2 || containsSymbol(arguments[0], variable))
            return std::nullopt;
        auto child = expandPuiseuxSeries(
            arguments[1], variable, center, integerCeiling + 2,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        const std::int64_t tickCeiling =
            integerCeiling * static_cast<std::int64_t>(child->denominator);
        if (tickCeiling > kMaximumInternalExponent)
            return std::nullopt;
        auto result = child->ticks.exactZero || child->ticks.minimumExponent > 0
            ? expandPolylogOriginSeries(
                arguments[0], child->ticks, tickCeiling,
                builtins, mathematics, angles, assumptions)
            : expandPolylogRegularCenterSeries(
                arguments[0], child->ticks, tickCeiling,
                builtins, mathematics, angles, assumptions);
        if (!result) return std::nullopt;
        return PuiseuxSeries{child->denominator, std::move(*result)};
    }
    case BuiltinId::Gamma:
    case BuiltinId::LogGamma: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandPuiseuxSeries(
            arguments[0], variable, center, integerCeiling + 2,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        const std::int64_t tickCeiling =
            integerCeiling * static_cast<std::int64_t>(child->denominator);
        if (tickCeiling > kMaximumInternalExponent) return std::nullopt;
        auto result = expandGammaFamilySeries(
            *id, child->ticks, tickCeiling,
            builtins, mathematics, angles, assumptions);
        if (!result) return std::nullopt;
        return PuiseuxSeries{child->denominator, std::move(*result)};
    }
    case BuiltinId::Digamma:
    case BuiltinId::Trigamma: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandPuiseuxSeries(
            arguments[0], variable, center, integerCeiling + 2,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        const std::int64_t tickCeiling =
            integerCeiling * static_cast<std::int64_t>(child->denominator);
        if (tickCeiling > kMaximumInternalExponent) return std::nullopt;
        auto result = expandPsiFamilySeries(
            *id, child->ticks, tickCeiling,
            builtins, mathematics, angles, assumptions);
        if (!result) return std::nullopt;
        return PuiseuxSeries{child->denominator, std::move(*result)};
    }
    case BuiltinId::LambertW: {
        const auto principal = lambertWPrincipalBranch(arguments);
        if (!principal) return std::nullopt;
        if (auto branchPoint = expandLambertWBranchPointPuiseux(
                arguments, variable, center, integerCeiling,
                builtins, mathematics, angles, assumptions))
            return branchPoint;
        auto child = expandPuiseuxSeries(
            arguments.back(), variable, center, integerCeiling + 2,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        const std::int64_t tickCeiling =
            integerCeiling * static_cast<std::int64_t>(child->denominator);
        if (tickCeiling > kMaximumInternalExponent) return std::nullopt;
        auto result = expandLambertWSeries(
            arguments, child->ticks, tickCeiling,
            builtins, mathematics, angles, assumptions);
        if (!result) return std::nullopt;
        return PuiseuxSeries{child->denominator, std::move(*result)};
    }
    case BuiltinId::Asin:
    case BuiltinId::Acos:
    case BuiltinId::Atan: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandPuiseuxSeries(
            arguments[0], variable, center, integerCeiling + 2,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        const std::int64_t tickCeiling =
            integerCeiling * static_cast<std::int64_t>(child->denominator);
        if (tickCeiling > kMaximumInternalExponent) return std::nullopt;
        auto result = expandInverseTrigSeries(
            *id, child->ticks, tickCeiling,
            builtins, mathematics, angles, assumptions);
        if (!result) return std::nullopt;
        return PuiseuxSeries{child->denominator, std::move(*result)};
    }
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi: {
        if (arguments.size() != 1) return std::nullopt;
        auto child = expandPuiseuxSeries(
            arguments[0], variable, center, integerCeiling + 2,
            builtins, mathematics, angles, assumptions);
        if (!child) return std::nullopt;
        const std::int64_t tickCeiling =
            integerCeiling * static_cast<std::int64_t>(child->denominator);
        if (tickCeiling > kMaximumInternalExponent) return std::nullopt;
        auto result = expandClassicalIntegralSeries(
            *id, child->ticks, tickCeiling,
            builtins, mathematics, angles, assumptions);
        if (!result) return std::nullopt;
        return PuiseuxSeries{child->denominator, std::move(*result)};
    }
    default:
        return std::nullopt;
    }
}


} // namespace

std::optional<SeriesData> parseSeriesData(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!expression.isCall()) return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition || definition->id != BuiltinId::SeriesData)
        return std::nullopt;

    const auto& arguments = expression.asCall().arguments;
    if ((arguments.size() != 6 && arguments.size() != 7) || !arguments[0].isSymbol())
        return std::nullopt;
    auto coefficients = braceElements(arguments[2]);
    const auto minimumExponent = exactInt64(arguments[3]);
    const auto orderNumerator = exactInt64(arguments[4]);
    const auto denominator = exactUint32(arguments[5]);
    if (!coefficients || !minimumExponent || !orderNumerator || !denominator)
        return std::nullopt;
    if (*denominator == 0 || *denominator > 64
        || *orderNumerator <= *minimumExponent)
        return std::nullopt;
    const auto expectedSize = static_cast<std::uint64_t>(
        *orderNumerator - *minimumExponent);
    if (expectedSize != static_cast<std::uint64_t>(coefficients->size()))
        return std::nullopt;

    std::vector<std::vector<Expr>> logarithmicCoefficients;
    if (arguments.size() == 7) {
        if (arguments[6].isArray()) {
            const auto& array = arguments[6].asArray();
            if (array.shape.size() != 2
                || array.shape[1] != static_cast<std::size_t>(expectedSize))
                return std::nullopt;
            const auto flat = array.materialize();
            logarithmicCoefficients.reserve(array.shape[0]);
            for (std::size_t row = 0; row < array.shape[0]; ++row) {
                const auto begin = flat.begin() + static_cast<std::ptrdiff_t>(row * array.shape[1]);
                logarithmicCoefficients.emplace_back(begin, begin + static_cast<std::ptrdiff_t>(array.shape[1]));
            }
        }
        else if (arguments[6].isList()) {
            logarithmicCoefficients.reserve(arguments[6].asList().elements.size());
            for (const Expr& layerExpression : arguments[6].asList().elements) {
                auto layer = braceElements(layerExpression);
                if (!layer || expectedSize != static_cast<std::uint64_t>(layer->size()))
                    return std::nullopt;
                logarithmicCoefficients.push_back(std::move(*layer));
            }
        }
        else {
            return std::nullopt;
        }
    }

    return SeriesData{
        arguments[0].asSymbol(), arguments[1], std::move(*coefficients),
        *minimumExponent, *orderNumerator, *denominator,
        std::move(logarithmicCoefficients)};
}

[[nodiscard]] bool coefficientLayerIsZero(const std::vector<Expr>& layer) noexcept {
    return std::all_of(layer.begin(), layer.end(), [](const Expr& coefficient) {
        return isZero(coefficient);
    });
}

void trimTrailingLogarithmicLayers(SeriesData& series) {
    while (!series.logarithmicCoefficients.empty()
        && coefficientLayerIsZero(series.logarithmicCoefficients.back()))
        series.logarithmicCoefficients.pop_back();
}

Expr makeSeriesData(
    SeriesData data,
    const evaluation::BuiltinRegistry& builtins) {
    trimTrailingLogarithmicLayers(data);
    std::vector<Expr> coefficients = std::move(data.coefficients);
    std::vector<Expr> arguments{
        Expr{data.variable},
        std::move(data.center),
        expression::braceValue(std::move(coefficients)),
        integer(data.minimumExponent),
        integer(data.orderNumerator),
        integer(static_cast<std::int64_t>(data.exponentDenominator))
    };
    if (!data.logarithmicCoefficients.empty()) {
        std::vector<Expr> layers;
        layers.reserve(data.logarithmicCoefficients.size());
        for (auto& layer : data.logarithmicCoefficients)
            layers.push_back(expression::braceValue(std::move(layer)));
        arguments.push_back(expression::braceValue(std::move(layers)));
    }
    return Expr::call(builtins.symbol(BuiltinId::SeriesData), std::move(arguments));
}

void trimSeriesDataLeadingZeros(SeriesData& series);

[[nodiscard]] std::optional<SeriesData> directLogarithmicOriginSeriesData(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!(center.isNumber() && center.asNumber().isZero()) || !expression.isCall())
        return std::nullopt;

    const auto* definition = builtins.find(expression.asCall().head);
    const auto& arguments = expression.asCall().arguments;
    if (!definition || arguments.size() != 1 || !arguments[0].isSymbol()
        || arguments[0].asSymbol() != variable
        || (definition->id != BuiltinId::Log
            && definition->id != BuiltinId::ExponentialIntegralEi
            && definition->id != BuiltinId::CosineIntegralCi))
        return std::nullopt;

    std::vector<Expr> coefficients(order + 1, integer(0));
    std::vector<Expr> logarithmic(order + 1, integer(0));
    logarithmic[0] = integer(1);

    if (definition->id == BuiltinId::ExponentialIntegralEi
        || definition->id == BuiltinId::CosineIntegralCi) {
        coefficients[0] = negate(
            simplify(call(BuiltinId::Digamma, {integer(1)}, builtins),
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions);
    }

    if (definition->id == BuiltinId::ExponentialIntegralEi && order >= 1) {
        Expr coefficient = integer(1);
        coefficients[1] = coefficient;
        for (std::size_t n = 2; n <= order; ++n) {
            const auto ni = static_cast<std::int64_t>(n);
            coefficient = multiply(
                coefficient, rational(ni - 1, ni * ni),
                builtins, mathematics, angles, assumptions);
            coefficients[n] = coefficient;
        }
    }
    else if (definition->id == BuiltinId::CosineIntegralCi && order >= 2) {
        Expr coefficient = rational(-1, 4);
        coefficients[2] = coefficient;
        for (std::size_t m = 2; 2 * m <= order; ++m) {
            const auto twoM = static_cast<std::int64_t>(2 * m);
            coefficient = multiply(
                coefficient,
                rational(-(twoM - 2), twoM * twoM * (twoM - 1)),
                builtins, mathematics, angles, assumptions);
            coefficients[2 * m] = coefficient;
        }
    }

    return SeriesData{
        variable, center, std::move(coefficients), 0,
        static_cast<std::int64_t>(order + 1), 1,
        {std::move(logarithmic)}};
}

[[nodiscard]] std::optional<SeriesData> pureSeriesData(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const std::int64_t ceiling = static_cast<std::int64_t>(order);
    if (auto expanded = expandSeries(
            expression, variable, center, ceiling,
            builtins, mathematics, angles, assumptions)) {
        const std::int64_t orderExclusive = static_cast<std::int64_t>(order + 1);
        if (expanded->exactZero)
            return SeriesData{
                variable, center, std::vector<Expr>(order + 1, integer(0)),
                0, orderExclusive, 1, {}};

        for (Expr& coefficient : expanded->coefficients)
            coefficient = simplify(
                std::move(coefficient), builtins, mathematics, angles, assumptions);
        trimLeadingZeros(*expanded);
        const std::int64_t minimum = expanded->minimumExponent;
        if (minimum >= orderExclusive)
            return SeriesData{
                variable, center, std::vector<Expr>(order + 1, integer(0)),
                0, orderExclusive, 1, {}};

        std::vector<Expr> coefficients;
        coefficients.reserve(static_cast<std::size_t>(orderExclusive - minimum));
        for (std::int64_t exponent = minimum; exponent < orderExclusive; ++exponent)
            coefficients.push_back(coefficientAt(*expanded, exponent));
        return SeriesData{
            variable, center, std::move(coefficients), minimum, orderExclusive, 1, {}};
    }

    const std::int64_t guardedOrder = std::min<std::int64_t>(
        1024, static_cast<std::int64_t>(order) + 4);
    auto puiseux = expandPuiseuxSeries(
        expression, variable, center, guardedOrder,
        builtins, mathematics, angles, assumptions);
    if (!puiseux || puiseux->denominator <= 1)
        return std::nullopt;

    for (Expr& coefficient : puiseux->ticks.coefficients)
        coefficient = simplify(
            std::move(coefficient), builtins, mathematics, angles, assumptions);
    trimLeadingZeros(puiseux->ticks);

    const std::int64_t denominator = static_cast<std::int64_t>(puiseux->denominator);
    const std::int64_t orderExclusive = static_cast<std::int64_t>(order) * denominator + 1;
    const std::int64_t minimum = puiseux->ticks.exactZero
        ? 0 : puiseux->ticks.minimumExponent;
    if (minimum >= orderExclusive)
        return SeriesData{
            variable, center,
            std::vector<Expr>(static_cast<std::size_t>(orderExclusive), integer(0)),
            0, orderExclusive, puiseux->denominator, {}};

    std::vector<Expr> coefficients;
    coefficients.reserve(static_cast<std::size_t>(orderExclusive - minimum));
    for (std::int64_t tick = minimum; tick < orderExclusive; ++tick)
        coefficients.push_back(coefficientAt(puiseux->ticks, tick));
    return SeriesData{
        variable, center, std::move(coefficients), minimum,
        orderExclusive, puiseux->denominator, {}};
}

[[nodiscard]] std::size_t maximumLogDegree(const SeriesData& series) noexcept {
    return series.logarithmicCoefficients.size();
}

[[nodiscard]] const Expr& seriesDataCoefficient(
    const SeriesData& series,
    std::size_t logDegree,
    std::size_t index) {
    return logDegree == 0
        ? series.coefficients[index]
        : series.logarithmicCoefficients[logDegree - 1][index];
}

void addSeriesDataContribution(
    SeriesData& series,
    std::size_t logDegree,
    std::size_t index,
    const Expr& contribution,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (isZero(contribution)) return;
    Expr& destination = logDegree == 0
        ? series.coefficients[index]
        : series.logarithmicCoefficients[logDegree - 1][index];
    destination = add(
        destination, contribution,
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] std::optional<std::uint32_t> commonSeriesDenominator(
    std::uint32_t lhs,
    std::uint32_t rhs) noexcept {
    const std::uint64_t divisor = std::gcd(lhs, rhs);
    const std::uint64_t value = static_cast<std::uint64_t>(lhs) / divisor * rhs;
    if (value == 0 || value > 64) return std::nullopt;
    return static_cast<std::uint32_t>(value);
}

[[nodiscard]] std::optional<SeriesData> addSeriesData(
    const SeriesData& lhs,
    const SeriesData& rhs,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto denominator = commonSeriesDenominator(
        lhs.exponentDenominator, rhs.exponentDenominator);
    if (!denominator) return std::nullopt;
    const std::int64_t lhsScale = *denominator / lhs.exponentDenominator;
    const std::int64_t rhsScale = *denominator / rhs.exponentDenominator;
    const std::int64_t minimum = std::min(
        lhs.minimumExponent * lhsScale,
        rhs.minimumExponent * rhsScale);
    const std::int64_t orderExclusive =
        static_cast<std::int64_t>(order) * static_cast<std::int64_t>(*denominator) + 1;
    if (orderExclusive <= minimum
        || orderExclusive - minimum > kMaximumInternalExponent)
        return std::nullopt;

    const std::size_t count = static_cast<std::size_t>(orderExclusive - minimum);
    const std::size_t logDegree = std::max(maximumLogDegree(lhs), maximumLogDegree(rhs));
    SeriesData result{
        lhs.variable, lhs.center, std::vector<Expr>(count, integer(0)),
        minimum, orderExclusive, *denominator,
        std::vector<std::vector<Expr>>(logDegree, std::vector<Expr>(count, integer(0)))};

    const auto accumulate = [&](const SeriesData& source, std::int64_t scale) {
        for (std::size_t i = 0; i < source.coefficients.size(); ++i) {
            const std::int64_t tick =
                (source.minimumExponent + static_cast<std::int64_t>(i)) * scale;
            if (tick < minimum || tick >= orderExclusive) continue;
            const std::size_t destination = static_cast<std::size_t>(tick - minimum);
            for (std::size_t degree = 0; degree <= maximumLogDegree(source); ++degree)
                addSeriesDataContribution(
                    result, degree, destination,
                    seriesDataCoefficient(source, degree, i),
                    builtins, mathematics, angles, assumptions);
        }
    };
    accumulate(lhs, lhsScale);
    accumulate(rhs, rhsScale);
    trimSeriesDataLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<SeriesData> multiplySeriesData(
    const SeriesData& lhs,
    const SeriesData& rhs,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    const auto denominator = commonSeriesDenominator(
        lhs.exponentDenominator, rhs.exponentDenominator);
    if (!denominator) return std::nullopt;
    const std::int64_t lhsScale = *denominator / lhs.exponentDenominator;
    const std::int64_t rhsScale = *denominator / rhs.exponentDenominator;
    const std::int64_t minimum = lhs.minimumExponent * lhsScale
        + rhs.minimumExponent * rhsScale;
    const std::int64_t orderExclusive =
        static_cast<std::int64_t>(order) * static_cast<std::int64_t>(*denominator) + 1;
    if (orderExclusive <= minimum
        || orderExclusive - minimum > kMaximumInternalExponent)
        return std::nullopt;

    const std::size_t count = static_cast<std::size_t>(orderExclusive - minimum);
    const std::size_t logDegree = maximumLogDegree(lhs) + maximumLogDegree(rhs);
    SeriesData result{
        lhs.variable, lhs.center, std::vector<Expr>(count, integer(0)),
        minimum, orderExclusive, *denominator,
        std::vector<std::vector<Expr>>(logDegree, std::vector<Expr>(count, integer(0)))};

    for (std::size_t i = 0; i < lhs.coefficients.size(); ++i) {
        const std::int64_t lhsTick =
            (lhs.minimumExponent + static_cast<std::int64_t>(i)) * lhsScale;
        for (std::size_t j = 0; j < rhs.coefficients.size(); ++j) {
            const std::int64_t tick = lhsTick
                + (rhs.minimumExponent + static_cast<std::int64_t>(j)) * rhsScale;
            if (tick < minimum || tick >= orderExclusive) continue;
            const std::size_t destination = static_cast<std::size_t>(tick - minimum);
            for (std::size_t lhsDegree = 0; lhsDegree <= maximumLogDegree(lhs); ++lhsDegree) {
                const Expr& lhsCoefficient = seriesDataCoefficient(lhs, lhsDegree, i);
                if (isZero(lhsCoefficient)) continue;
                for (std::size_t rhsDegree = 0; rhsDegree <= maximumLogDegree(rhs); ++rhsDegree) {
                    const Expr& rhsCoefficient = seriesDataCoefficient(rhs, rhsDegree, j);
                    if (isZero(rhsCoefficient)) continue;
                    addSeriesDataContribution(
                        result, lhsDegree + rhsDegree, destination,
                        multiply(lhsCoefficient, rhsCoefficient,
                            builtins, mathematics, angles, assumptions),
                        builtins, mathematics, angles, assumptions);
                }
            }
        }
    }
    trimSeriesDataLeadingZeros(result);
    return result;
}

[[nodiscard]] std::optional<SeriesData> composedLogarithmicOriginSeriesData(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!(center.isNumber() && center.asNumber().isZero()) || !expression.isCall())
        return std::nullopt;

    const auto* definition = builtins.find(expression.asCall().head);
    const auto& arguments = expression.asCall().arguments;
    if (!definition || arguments.size() != 1
        || (definition->id != BuiltinId::Log
            && definition->id != BuiltinId::ExponentialIntegralEi
            && definition->id != BuiltinId::CosineIntegralCi))
        return std::nullopt;

    // A(t)=c*t^r*(1+h)で0<r<=1かつc>0を証明できる場合だけ，
    // principal log A=log c+r log t+log(1+h)としてlog係数層へ分解する。
    const std::int64_t probeCeiling = std::min<std::int64_t>(
        1024, static_cast<std::int64_t>(order) + 8);
    auto argument = expandPuiseuxSeries(
        arguments[0], variable, center, probeCeiling,
        builtins, mathematics, angles, assumptions);
    if (!argument || argument->ticks.exactZero || argument->ticks.minimumExponent <= 0)
        return std::nullopt;

    const std::int64_t denominator = static_cast<std::int64_t>(argument->denominator);
    const std::int64_t leadingExponent = argument->ticks.minimumExponent;
    if (leadingExponent > denominator)
        return std::nullopt;
    const Expr leading = coefficientAt(argument->ticks, leadingExponent);
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.facts(leading).sign != mathematics::RealSign::Positive)
        return std::nullopt;

    const std::int64_t targetTick =
        static_cast<std::int64_t>(order) * denominator;
    if (targetTick > kMaximumInternalExponent
        || leadingExponent > kMaximumInternalExponent - targetTick)
        return std::nullopt;
    const std::int64_t requiredTick = leadingExponent + targetTick;
    const std::int64_t requiredIntegerCeiling =
        (requiredTick + denominator - 1) / denominator;
    if (requiredIntegerCeiling > 1024)
        return std::nullopt;
    if (argument->ticks.ceiling < requiredTick) {
        argument = expandPuiseuxSeries(
            arguments[0], variable, center, requiredIntegerCeiling,
            builtins, mathematics, angles, assumptions);
        if (!argument || argument->ticks.exactZero
            || argument->ticks.minimumExponent != leadingExponent
            || argument->denominator != static_cast<std::uint32_t>(denominator))
            return std::nullopt;
    }

    LaurentSeries unit = denseSeries(0, targetTick);
    for (std::int64_t tick = 0; tick <= targetTick; ++tick) {
        unit.coefficients[static_cast<std::size_t>(tick)] = divide(
            coefficientAt(argument->ticks, leadingExponent + tick), leading,
            builtins, mathematics, angles, assumptions);
    }
    auto regularLog = expandLogSeries(
        unit, targetTick, builtins, mathematics, angles, assumptions);
    if (!regularLog)
        return std::nullopt;

    const std::size_t count = static_cast<std::size_t>(targetTick + 1);
    SeriesData result{
        variable, center, std::vector<Expr>(count, integer(0)),
        0, targetTick + 1, argument->denominator,
        {std::vector<Expr>(count, integer(0))}};
    for (std::int64_t tick = 0; tick <= targetTick; ++tick)
        result.coefficients[static_cast<std::size_t>(tick)] =
            coefficientAt(*regularLog, tick);
    result.coefficients[0] = add(
        result.coefficients[0],
        simplify(call(BuiltinId::Log, {leading}, builtins),
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    result.logarithmicCoefficients[0][0] =
        rational(leadingExponent, denominator);

    if (definition->id == BuiltinId::Log)
        return result;

    result.coefficients[0] = add(
        result.coefficients[0],
        negate(
            simplify(call(BuiltinId::Digamma, {integer(1)}, builtins),
                builtins, mathematics, angles, assumptions),
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    if (targetTick == 0)
        return result;

    const std::int64_t derivativeCeiling = targetTick - 1;
    const std::int64_t kernelWorkCeiling = derivativeCeiling + leadingExponent;
    LaurentSeries numerator = zeroSeries(kernelWorkCeiling);
    if (definition->id == BuiltinId::ExponentialIntegralEi) {
        auto exponential = expandExpSeries(
            argument->ticks, kernelWorkCeiling,
            builtins, mathematics, angles, assumptions);
        if (!exponential) return std::nullopt;
        numerator = std::move(*exponential);
    }
    else {
        auto sinCos = expandSinCosSeries(
            argument->ticks, kernelWorkCeiling,
            builtins, mathematics, angles, assumptions);
        if (!sinCos) return std::nullopt;
        numerator = std::move(sinCos->second);
    }
    if (numerator.exactZero || numerator.minimumExponent > 0)
        return std::nullopt;
    const std::size_t constantIndex = static_cast<std::size_t>(-numerator.minimumExponent);
    if (constantIndex >= numerator.coefficients.size())
        return std::nullopt;
    numerator.coefficients[constantIndex] = subtract(
        numerator.coefficients[constantIndex], integer(1),
        builtins, mathematics, angles, assumptions);
    trimLeadingZeros(numerator);

    auto inverse = inverseSeries(
        argument->ticks, derivativeCeiling,
        builtins, mathematics, angles, assumptions, true);
    if (!inverse) return std::nullopt;
    LaurentSeries kernel = multiplySeries(
        numerator, *inverse, derivativeCeiling,
        builtins, mathematics, angles, assumptions);
    LaurentSeries argumentDerivative = differentiateFormalSeries(
        argument->ticks, derivativeCeiling,
        builtins, mathematics, angles, assumptions);
    LaurentSeries derivative = multiplySeries(
        kernel, argumentDerivative, derivativeCeiling,
        builtins, mathematics, angles, assumptions);
    auto regularCorrection = integrateFormalDerivative(
        derivative, integer(0), targetTick,
        builtins, mathematics, angles, assumptions);
    if (!regularCorrection) return std::nullopt;
    for (std::int64_t tick = 0; tick <= targetTick; ++tick) {
        const Expr correction = coefficientAt(*regularCorrection, tick);
        if (isZero(correction)) continue;
        result.coefficients[static_cast<std::size_t>(tick)] = add(
            result.coefficients[static_cast<std::size_t>(tick)], correction,
            builtins, mathematics, angles, assumptions);
    }
    return result;
}

[[nodiscard]] bool containsLogHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins);

[[nodiscard]] std::optional<SeriesData> logarithmicSeriesData(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::size_t depth = 0) {
    if (depth > 32) return std::nullopt;
    if (expression.isCall())
        if (const auto* rewriteDefinition = builtins.find(expression.asCall().head))
            if (auto rewritten = lowCostSeriesRewrite(
                    rewriteDefinition->id, expression.asCall().arguments,
                    builtins, mathematics, angles, assumptions))
                return logarithmicSeriesData(
                    *rewritten, variable, center, order,
                    builtins, mathematics, angles, assumptions, depth + 1);
    if (auto direct = directLogarithmicOriginSeriesData(
            expression, variable, center, order,
            builtins, mathematics, angles, assumptions))
        return direct;
    if (auto pure = pureSeriesData(
            expression, variable, center, order,
            builtins, mathematics, angles, assumptions))
        return pure;
    if (auto composed = composedLogarithmicOriginSeriesData(
            expression, variable, center, order,
            builtins, mathematics, angles, assumptions))
        return composed;
    if (!expression.isCall()) return std::nullopt;

    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition) return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    switch (definition->id) {
    case BuiltinId::Negate: {
        if (arguments.size() != 1) return std::nullopt;
        auto result = logarithmicSeriesData(
            arguments[0], variable, center, order,
            builtins, mathematics, angles, assumptions, depth + 1);
        if (!result) return std::nullopt;
        for (Expr& coefficient : result->coefficients)
            coefficient = negate(
                coefficient, builtins, mathematics, angles, assumptions);
        for (auto& layer : result->logarithmicCoefficients)
            for (Expr& coefficient : layer)
                coefficient = negate(
                    coefficient, builtins, mathematics, angles, assumptions);
        return result;
    }
    case BuiltinId::Add: {
        if (arguments.empty()) return std::nullopt;
        auto result = logarithmicSeriesData(
            arguments[0], variable, center, order,
            builtins, mathematics, angles, assumptions, depth + 1);
        if (!result) return std::nullopt;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            auto term = logarithmicSeriesData(
                arguments[i], variable, center, order,
                builtins, mathematics, angles, assumptions, depth + 1);
            if (!term) return std::nullopt;
            auto sum = addSeriesData(
                *result, *term, order,
                builtins, mathematics, angles, assumptions);
            if (!sum) return std::nullopt;
            result = std::move(*sum);
        }
        return result;
    }
    case BuiltinId::Subtract: {
        if (arguments.size() != 2) return std::nullopt;
        auto lhs = logarithmicSeriesData(
            arguments[0], variable, center, order,
            builtins, mathematics, angles, assumptions, depth + 1);
        auto rhs = logarithmicSeriesData(
            arguments[1], variable, center, order,
            builtins, mathematics, angles, assumptions, depth + 1);
        if (!lhs || !rhs) return std::nullopt;
        for (Expr& coefficient : rhs->coefficients)
            coefficient = negate(
                coefficient, builtins, mathematics, angles, assumptions);
        for (auto& layer : rhs->logarithmicCoefficients)
            for (Expr& coefficient : layer)
                coefficient = negate(
                    coefficient, builtins, mathematics, angles, assumptions);
        return addSeriesData(
            *lhs, *rhs, order, builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Multiply: {
        if (arguments.empty()) return std::nullopt;
        auto result = logarithmicSeriesData(
            arguments[0], variable, center, order,
            builtins, mathematics, angles, assumptions, depth + 1);
        if (!result) return std::nullopt;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            auto factor = logarithmicSeriesData(
                arguments[i], variable, center, order,
                builtins, mathematics, angles, assumptions, depth + 1);
            if (!factor) return std::nullopt;
            auto product = multiplySeriesData(
                *result, *factor, order,
                builtins, mathematics, angles, assumptions);
            if (!product) return std::nullopt;
            result = std::move(*product);
        }
        return result;
    }
    case BuiltinId::Divide: {
        if (arguments.size() != 2) return std::nullopt;
        auto lhs = logarithmicSeriesData(
            arguments[0], variable, center, order,
            builtins, mathematics, angles, assumptions, depth + 1);
        if (!lhs) return std::nullopt;

        if (!containsSymbol(arguments[1], variable)) {
            if (!seriesCenterValueProvablyNonZero(
                    arguments[1], variable, center,
                    builtins, mathematics, angles, assumptions))
                return std::nullopt;
            for (Expr& coefficient : lhs->coefficients)
                coefficient = divide(
                    coefficient, arguments[1],
                    builtins, mathematics, angles, assumptions);
            for (auto& layer : lhs->logarithmicCoefficients)
                for (Expr& coefficient : layer)
                    coefficient = divide(
                        coefficient, arguments[1],
                        builtins, mathematics, angles, assumptions);
            return lhs;
        }

        if (containsLogHead(arguments[1], builtins)) return std::nullopt;
        const Expr reciprocal = call(
            BuiltinId::Divide, {integer(1), arguments[1]}, builtins);
        auto rhsInverse = pureSeriesData(
            reciprocal, variable, center, order,
            builtins, mathematics, angles, assumptions);
        if (!rhsInverse) return std::nullopt;
        return multiplySeriesData(
            *lhs, *rhsInverse, order,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Power: {
        if (arguments.size() != 2) return std::nullopt;
        const auto exponent = exactInt64(arguments[1]);
        if (!exponent || *exponent < 0) return std::nullopt;
        auto base = logarithmicSeriesData(
            arguments[0], variable, center, order,
            builtins, mathematics, angles, assumptions, depth + 1);
        if (!base) return std::nullopt;
        SeriesData result{
            variable, center, {integer(1)}, 0,
            static_cast<std::int64_t>(order + 1), 1, {}};
        std::int64_t power = *exponent;
        while (power > 0) {
            if ((power & 1) != 0) {
                auto product = multiplySeriesData(
                    result, *base, order,
                    builtins, mathematics, angles, assumptions);
                if (!product) return std::nullopt;
                result = std::move(*product);
            }
            power >>= 1;
            if (power == 0) break;
            auto square = multiplySeriesData(
                *base, *base, order,
                builtins, mathematics, angles, assumptions);
            if (!square) return std::nullopt;
            base = std::move(*square);
        }
        return result;
    }
    default:
        return std::nullopt;
    }
}

[[nodiscard]] bool containsLogHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (expression.isCall()) {
        if (const auto* definition = builtins.find(expression.asCall().head); definition) {
            switch (definition->id) {
            case BuiltinId::Log:
            case BuiltinId::Log1p:
            case BuiltinId::Log2:
            case BuiltinId::Log10:
                return true;
            default:
                break;
            }
        }
        for (const Expr& argument : expression.asCall().arguments)
            if (containsLogHead(argument, builtins)) return true;
    }
    if (expression.isList()) {
        for (const Expr& element : expression.asList().elements)
            if (containsLogHead(element, builtins)) return true;
    }
    return false;
}

[[nodiscard]] std::optional<SeriesData> positiveInfinityLogSeriesData(
    const Expr& expression,
    const expression::Symbol& variable,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (!expression.isCall()) return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    const auto& arguments = expression.asCall().arguments;
    if (!definition || definition->id != BuiltinId::Log || arguments.size() != 1
        || containsLogHead(arguments[0], builtins))
        return std::nullopt;

    auto argumentExpression = seriesExpression(
        arguments[0], variable, Expr{expression::Symbol{"Infinity"}}, order,
        builtins, mathematics, angles, assumptions);
    if (!argumentExpression) return std::nullopt;
    auto argument = parseSeriesData(*argumentExpression, builtins);
    if (!argument || !argument->logarithmicCoefficients.empty()
        || argument->coefficients.empty())
        return std::nullopt;
    trimSeriesDataLeadingZeros(*argument);
    if (argument->coefficients.empty() || isZero(argument->coefficients[0]))
        return std::nullopt;

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    const Expr leading = argument->coefficients[0];
    if (knowledge.facts(leading).sign != mathematics::RealSign::Positive)
        return std::nullopt;

    const std::int64_t denominator = static_cast<std::int64_t>(argument->exponentDenominator);
    if (denominator <= 0 || order > static_cast<std::size_t>(kMaximumInternalExponent / denominator))
        return std::nullopt;
    const std::int64_t targetTick = static_cast<std::int64_t>(order) * denominator;
    LaurentSeries unit = denseSeries(0, targetTick);
    for (std::int64_t tick = 0; tick <= targetTick; ++tick) {
        const auto index = static_cast<std::size_t>(tick);
        const Expr coefficient = index < argument->coefficients.size()
            ? argument->coefficients[index] : integer(0);
        unit.coefficients[index] = divide(
            coefficient, leading, builtins, mathematics, angles, assumptions);
    }
    auto regularLog = expandLogSeries(
        unit, targetTick, builtins, mathematics, angles, assumptions);
    if (!regularLog) return std::nullopt;

    const std::size_t count = static_cast<std::size_t>(targetTick + 1);
    SeriesData result{
        variable, Expr{expression::Symbol{"Infinity"}},
        std::vector<Expr>(count, integer(0)), 0, targetTick + 1,
        argument->exponentDenominator,
        {std::vector<Expr>(count, integer(0))}};
    for (std::int64_t tick = 0; tick <= targetTick; ++tick)
        result.coefficients[static_cast<std::size_t>(tick)] = coefficientAt(*regularLog, tick);
    result.coefficients[0] = add(
        result.coefficients[0],
        simplify(call(BuiltinId::Log, {leading}, builtins),
            builtins, mathematics, angles, assumptions),
        builtins, mathematics, angles, assumptions);
    result.logarithmicCoefficients[0][0] = rational(
        argument->minimumExponent, denominator);
    trimTrailingLogarithmicLayers(result);
    return result;
}

[[nodiscard]] std::optional<SeriesData> positiveInfinityLogarithmicSeriesData(
    const Expr& expression,
    const expression::Symbol& variable,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions,
    std::size_t depth = 0) {
    if (depth > 32) return std::nullopt;
    if (expression.isCall())
        if (const auto* rewriteDefinition = builtins.find(expression.asCall().head))
            if (auto rewritten = lowCostSeriesRewrite(
                    rewriteDefinition->id, expression.asCall().arguments,
                    builtins, mathematics, angles, assumptions))
                return positiveInfinityLogarithmicSeriesData(
                    *rewritten, variable, order,
                    builtins, mathematics, angles, assumptions, depth + 1);
    if (auto logarithm = positiveInfinityLogSeriesData(
            expression, variable, order,
            builtins, mathematics, angles, assumptions))
        return logarithm;
    if (!expression.isCall())
        return std::nullopt;

    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition) return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    const Expr infinity{expression::Symbol{"Infinity"}};
    const auto pure = [&](const Expr& value) -> std::optional<SeriesData> {
        if (containsLogHead(value, builtins)) return std::nullopt;
        auto expanded = seriesExpression(
            value, variable, infinity, order,
            builtins, mathematics, angles, assumptions);
        if (!expanded) return std::nullopt;
        return parseSeriesData(*expanded, builtins);
    };
    const auto child = [&](const Expr& value) -> std::optional<SeriesData> {
        if (!containsLogHead(value, builtins)) return pure(value);
        return positiveInfinityLogarithmicSeriesData(
            value, variable, order,
            builtins, mathematics, angles, assumptions, depth + 1);
    };

    switch (definition->id) {
    case BuiltinId::Negate: {
        if (arguments.size() != 1) return std::nullopt;
        auto result = child(arguments[0]);
        if (!result) return std::nullopt;
        for (Expr& coefficient : result->coefficients)
            coefficient = negate(coefficient, builtins, mathematics, angles, assumptions);
        for (auto& layer : result->logarithmicCoefficients)
            for (Expr& coefficient : layer)
                coefficient = negate(coefficient, builtins, mathematics, angles, assumptions);
        return result;
    }
    case BuiltinId::Add: {
        if (arguments.empty()) return std::nullopt;
        auto result = child(arguments[0]);
        if (!result) return std::nullopt;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            auto term = child(arguments[i]);
            if (!term) return std::nullopt;
            auto sum = addSeriesData(
                *result, *term, order,
                builtins, mathematics, angles, assumptions);
            if (!sum) return std::nullopt;
            result = std::move(*sum);
        }
        return result;
    }
    case BuiltinId::Subtract: {
        if (arguments.size() != 2) return std::nullopt;
        auto lhs = child(arguments[0]);
        auto rhs = child(arguments[1]);
        if (!lhs || !rhs) return std::nullopt;
        for (Expr& coefficient : rhs->coefficients)
            coefficient = negate(coefficient, builtins, mathematics, angles, assumptions);
        for (auto& layer : rhs->logarithmicCoefficients)
            for (Expr& coefficient : layer)
                coefficient = negate(coefficient, builtins, mathematics, angles, assumptions);
        return addSeriesData(
            *lhs, *rhs, order, builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Multiply: {
        if (arguments.empty()) return std::nullopt;
        auto result = child(arguments[0]);
        if (!result) return std::nullopt;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            auto factor = child(arguments[i]);
            if (!factor) return std::nullopt;
            auto product = multiplySeriesData(
                *result, *factor, order,
                builtins, mathematics, angles, assumptions);
            if (!product) return std::nullopt;
            result = std::move(*product);
        }
        return result;
    }
    case BuiltinId::Divide: {
        if (arguments.size() != 2
            || (containsLogHead(arguments[1], builtins)
                && containsSymbol(arguments[1], variable)))
            return std::nullopt;
        auto lhs = child(arguments[0]);
        if (!lhs) return std::nullopt;

        if (!containsSymbol(arguments[1], variable)) {
            if (!seriesCenterValueProvablyNonZero(
                    arguments[1], variable, infinity,
                    builtins, mathematics, angles, assumptions))
                return std::nullopt;
            for (Expr& coefficient : lhs->coefficients)
                coefficient = divide(
                    coefficient, arguments[1],
                    builtins, mathematics, angles, assumptions);
            for (auto& layer : lhs->logarithmicCoefficients)
                for (Expr& coefficient : layer)
                    coefficient = divide(
                        coefficient, arguments[1],
                        builtins, mathematics, angles, assumptions);
            return lhs;
        }

        const Expr reciprocal = call(
            BuiltinId::Divide, {integer(1), arguments[1]}, builtins);
        auto rhsInverse = pure(reciprocal);
        if (!rhsInverse) return std::nullopt;
        return multiplySeriesData(
            *lhs, *rhsInverse, order,
            builtins, mathematics, angles, assumptions);
    }
    case BuiltinId::Power: {
        if (arguments.size() != 2) return std::nullopt;
        const auto exponent = exactInt64(arguments[1]);
        if (!exponent) return std::nullopt;
        if (*exponent < 0)
            return containsLogHead(arguments[0], builtins) ? std::nullopt : pure(expression);
        auto base = child(arguments[0]);
        if (!base) return std::nullopt;
        SeriesData result{
            variable, infinity, {integer(1)}, 0,
            static_cast<std::int64_t>(order + 1), 1, {}};
        std::int64_t power = *exponent;
        while (power > 0) {
            if ((power & 1) != 0) {
                auto product = multiplySeriesData(
                    result, *base, order,
                    builtins, mathematics, angles, assumptions);
                if (!product) return std::nullopt;
                result = std::move(*product);
            }
            power >>= 1;
            if (power == 0) break;
            auto square = multiplySeriesData(
                *base, *base, order,
                builtins, mathematics, angles, assumptions);
            if (!square) return std::nullopt;
            base = std::move(*square);
        }
        return result;
    }
    default:
        return std::nullopt;
    }
}

std::optional<Expr> seriesExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const Expr& center,
    std::size_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const mathematics::AssumptionSet& assumptions) {
    if (containsSymbol(center, variable))
        return std::nullopt;
    if (order > 1024)
        return std::nullopt;

    if (isPositiveInfinityCenter(center)) {
        // +Infinityではt->0+の方向が固定されるため，正の先頭方向のlogを安全に因数分解できる。
        if (containsLogHead(expression, builtins)) {
            if (auto logarithmic = positiveInfinityLogarithmicSeriesData(
                    expression, variable, order,
                    builtins, mathematics, angles, assumptions))
                return makeSeriesData(std::move(*logarithmic), builtins);
            return std::nullopt;
        }
        // x->+Infinityをt=1/x，t->0+の局所展開へ写し，既存TPSAを再利用する。
        const expression::Symbol temporary = temporarySeriesVariable(expression);
        const Expr reciprocal = divide(
            integer(1), Expr{temporary},
            builtins, mathematics, angles, assumptions);
        const Expr transformed = simplify(
            canonicalizeInfinityReciprocals(
                substituteSymbol(expression, variable, reciprocal), builtins),
            builtins, mathematics, angles, assumptions);
        auto local = seriesExpression(
            transformed, temporary, integer(0), order,
            builtins, mathematics, angles, assumptions);
        if (!local) return std::nullopt;
        auto data = parseSeriesData(*local, builtins);
        if (!data) return std::nullopt;
        data->variable = variable;
        data->center = center;
        return makeSeriesData(std::move(*data), builtins);
    }

    const Expr prepared = simplify(
        expression, builtins, mathematics, angles, assumptions);

    // 原点の対数特異性は通常のLaurent/Puiseux係数だけでは表せないため，
    // DLMF 6.6.1/6.6.6の局所形をlog係数層へ直接格納する。
    if (auto logarithmic = directLogarithmicOriginSeriesData(
            prepared, variable, center, order,
            builtins, mathematics, angles, assumptions))
        return makeSeriesData(std::move(*logarithmic), builtins);

    const std::int64_t ceiling = static_cast<std::int64_t>(order);
    if (auto expanded = expandSeries(
            prepared, variable, center, ceiling,
            builtins, mathematics, angles, assumptions)) {
        if (expanded->exactZero) {
            std::vector<Expr> coefficients(order + 1, integer(0));
            return makeSeriesData(SeriesData{
                variable, center, std::move(coefficients),
                0, static_cast<std::int64_t>(order + 1), 1, {}}, builtins);
        }

        for (Expr& coefficient : expanded->coefficients)
            coefficient = simplify(
                std::move(coefficient), builtins, mathematics, angles, assumptions);
        trimLeadingZeros(*expanded);

        const std::int64_t minimum = expanded->exactZero ? 0 : expanded->minimumExponent;
        const std::int64_t orderExclusive = static_cast<std::int64_t>(order + 1);
        if (minimum >= orderExclusive) {
            std::vector<Expr> coefficients(order + 1, integer(0));
            return makeSeriesData(SeriesData{
                variable, center, std::move(coefficients),
                0, orderExclusive, 1, {}}, builtins);
        }

        std::vector<Expr> coefficients;
        coefficients.reserve(static_cast<std::size_t>(orderExclusive - minimum));
        for (std::int64_t exponent = minimum; exponent < orderExclusive; ++exponent)
            coefficients.push_back(coefficientAt(*expanded, exponent));

        return makeSeriesData(SeriesData{
            variable, center, std::move(coefficients),
            minimum, orderExclusive, 1, {}}, builtins);
    }

    const std::int64_t guardedOrder = std::min<std::int64_t>(
        1024, static_cast<std::int64_t>(order) + 4);
    auto puiseux = expandPuiseuxSeries(
        prepared, variable, center, guardedOrder,
        builtins, mathematics, angles, assumptions);
    if (!puiseux || puiseux->denominator <= 1) {
        if (auto logarithmic = logarithmicSeriesData(
                prepared, variable, center, order,
                builtins, mathematics, angles, assumptions))
            return makeSeriesData(std::move(*logarithmic), builtins);
        return std::nullopt;
    }

    for (Expr& coefficient : puiseux->ticks.coefficients)
        coefficient = simplify(
            std::move(coefficient), builtins, mathematics, angles, assumptions);
    trimLeadingZeros(puiseux->ticks);

    const std::int64_t denominator = static_cast<std::int64_t>(puiseux->denominator);
    const std::int64_t orderExclusive =
        static_cast<std::int64_t>(order) * denominator + 1;
    const std::int64_t minimum = puiseux->ticks.exactZero
        ? 0 : puiseux->ticks.minimumExponent;
    if (minimum >= orderExclusive) {
        std::vector<Expr> coefficients(static_cast<std::size_t>(orderExclusive), integer(0));
        return makeSeriesData(SeriesData{
            variable, center, std::move(coefficients),
            0, orderExclusive, puiseux->denominator, {}}, builtins);
    }

    std::vector<Expr> coefficients;
    coefficients.reserve(static_cast<std::size_t>(orderExclusive - minimum));
    for (std::int64_t tick = minimum; tick < orderExclusive; ++tick)
        coefficients.push_back(coefficientAt(puiseux->ticks, tick));

    return makeSeriesData(SeriesData{
        variable, center, std::move(coefficients),
        minimum, orderExclusive, puiseux->denominator, {}}, builtins);
}

Expr normalSeriesExpression(
    const SeriesData& series,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const bool atInfinity = isPositiveInfinityCenter(series.center);
    Expr base = atInfinity
        ? Expr{series.variable}
        : (series.center.isNumber() && series.center.asNumber().isZero()
            ? Expr{series.variable}
            : builtins::evaluateSubtract(
                std::array<Expr, 2>{Expr{series.variable}, series.center}, builtins));
    const Expr logarithm = simplify(
        call(BuiltinId::Log, {base}, builtins), builtins, mathematics, angles, {});

    const auto powerFactor = [&](std::int64_t exponentNumerator) {
        if (exponentNumerator == 0)
            return integer(1);
        const std::int64_t displayedExponent = atInfinity
            ? -exponentNumerator : exponentNumerator;
        return builtins::evaluatePower(
            std::array<Expr, 2>{
                base,
                rational(displayedExponent,
                    static_cast<std::int64_t>(series.exponentDenominator))},
            builtins, mathematics);
    };
    const auto logFactor = [&](std::size_t degree) {
        if (degree == 0) return integer(1);
        Expr factor = degree == 1
            ? logarithm
            : builtins::evaluatePower(
                std::array<Expr, 2>{logarithm, integer(static_cast<std::int64_t>(degree))},
                builtins, mathematics);
        if (atInfinity && (degree & 1U) != 0)
            factor = negate(factor, builtins, mathematics, angles, {});
        return factor;
    };

    Expr result = integer(0);
    for (std::size_t i = 0; i < series.coefficients.size(); ++i) {
        const std::int64_t exponentNumerator =
            series.minimumExponent + static_cast<std::int64_t>(i);
        const Expr power = powerFactor(exponentNumerator);

        const auto appendTerm = [&](const Expr& coefficient, std::size_t logDegree) {
            if (isZero(coefficient)) return;
            Expr term = coefficient;
            if (exponentNumerator != 0)
                term = multiply(term, power, builtins, mathematics, angles, {});
            if (logDegree != 0)
                term = multiply(term, logFactor(logDegree),
                    builtins, mathematics, angles, {});
            result = add(result, term, builtins, mathematics, angles, {});
        };

        appendTerm(series.coefficients[i], 0);
        for (std::size_t layer = 0; layer < series.logarithmicCoefficients.size(); ++layer)
            appendTerm(series.logarithmicCoefficients[layer][i], layer + 1);
    }
    return simplify(std::move(result), builtins, mathematics, angles, {});
}

Expr toNormalExpression(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto series = parseSeriesData(expression, builtins))
        return toNormalExpression(
            normalSeriesExpression(*series, builtins, mathematics, angles),
            builtins, mathematics, angles);

    if (expression.isCall()) {
        const auto& source = expression.asCall();
        std::vector<Expr> arguments;
        arguments.reserve(source.arguments.size());
        bool changed = false;
        for (const Expr& argument : source.arguments) {
            Expr converted = toNormalExpression(argument, builtins, mathematics, angles);
            changed |= !(converted == argument);
            arguments.push_back(std::move(converted));
        }
        if (!changed) return expression;
        return Expr::rebuildCall(source, std::move(arguments));
    }

    if (expression.isList()) {
        const auto& source = expression.asList().elements;
        std::vector<Expr> elements;
        elements.reserve(source.size());
        bool changed = false;
        for (const Expr& element : source) {
            Expr converted = toNormalExpression(element, builtins, mathematics, angles);
            changed |= !(converted == element);
            elements.push_back(std::move(converted));
        }
        if (!changed) return expression;
        return Expr::list(std::move(elements));
    }

    if (expression.isArray() && expression.asArray().hasStoredExpressions()) {
        const auto entries = expression.asArray().expressionEntries();
        std::vector<std::size_t> indices;
        std::vector<Expr> values;
        indices.reserve(entries.size());
        values.reserve(entries.size());
        for (const auto& entry : entries) {
            Expr converted = toNormalExpression(
                entry.expression, builtins, mathematics, angles);
            if (converted == entry.expression) continue;
            indices.push_back(entry.index);
            values.push_back(std::move(converted));
        }
        if (indices.empty()) return expression;
        return Expr::array(expression.asArray().replacedExpressions(indices, std::move(values)));
    }

    if (expression.isSolutionSet()) {
        const solver::SolutionSet& solutions = expression.asSolutionSet();
        bool changed = false;
        const auto convertBranch = [&](const solver::SolutionBranch& source) {
            solver::SolutionBranch result = source;
            for (solver::SolutionBinding& binding : result.bindings) {
                Expr converted = toNormalExpression(
                    binding.value, builtins, mathematics, angles);
                changed |= !(converted == binding.value);
                binding.value = std::move(converted);
            }
            return result;
        };
        const auto variables = [&] {
            return std::vector<solver::SolverVariable>(
                solutions.variables().begin(), solutions.variables().end());
        };

        switch (solutions.kind()) {
        case solver::SolutionSetKind::Finite: {
            std::vector<solver::SolutionBranch> branches;
            branches.reserve(solutions.branches().size());
            for (const solver::SolutionBranch& branch : solutions.branches())
                branches.push_back(convertBranch(branch));
            if (!changed) return expression;
            return Expr::solutionSet(solver::SolutionSet::finite(variables(), std::move(branches)));
        }
        case solver::SolutionSetKind::Conditional: {
            std::vector<solver::SolutionCase> cases(
                solutions.cases().begin(), solutions.cases().end());
            for (solver::SolutionCase& solutionCase : cases)
                if (solutionCase.outcome == solver::SolutionSetKind::Finite)
                    for (solver::SolutionBranch& branch : solutionCase.branches)
                        branch = convertBranch(branch);
            if (!changed) return expression;
            return Expr::solutionSet(solver::SolutionSet::conditional(variables(), std::move(cases)));
        }
        case solver::SolutionSetKind::Empty:
        case solver::SolutionSetKind::Universal:
        case solver::SolutionSetKind::Unresolved:
            return expression;
        }
    }

    return expression;
}

[[nodiscard]] bool seriesIndexIsZero(const SeriesData& series, std::size_t index) noexcept {
    if (!isZero(series.coefficients[index])) return false;
    for (const auto& layer : series.logarithmicCoefficients)
        if (!isZero(layer[index])) return false;
    return true;
}

void trimSeriesDataLeadingZeros(SeriesData& series) {
    std::size_t leading = 0;
    while (leading < series.coefficients.size() && seriesIndexIsZero(series, leading))
        ++leading;
    if (leading == 0 || leading == series.coefficients.size()) {
        trimTrailingLogarithmicLayers(series);
        return;
    }
    series.minimumExponent += static_cast<std::int64_t>(leading);
    series.coefficients.erase(
        series.coefficients.begin(),
        series.coefficients.begin() + static_cast<std::ptrdiff_t>(leading));
    for (auto& layer : series.logarithmicCoefficients)
        layer.erase(layer.begin(), layer.begin() + static_cast<std::ptrdiff_t>(leading));
    trimTrailingLogarithmicLayers(series);
}

[[nodiscard]] Expr& logarithmicCoefficient(
    SeriesData& series,
    std::size_t logDegree,
    std::size_t index) {
    return logDegree == 0
        ? series.coefficients[index]
        : series.logarithmicCoefficients[logDegree - 1][index];
}

[[nodiscard]] const Expr& logarithmicCoefficient(
    const SeriesData& series,
    std::size_t logDegree,
    std::size_t index) {
    return logDegree == 0
        ? series.coefficients[index]
        : series.logarithmicCoefficients[logDegree - 1][index];
}

void addLogarithmicCoefficient(
    SeriesData& series,
    std::size_t logDegree,
    std::size_t index,
    const Expr& contribution,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (isZero(contribution)) return;
    Expr& destination = logarithmicCoefficient(series, logDegree, index);
    destination = simplify(
        add(destination, contribution, builtins, mathematics, angles, {}),
        builtins, mathematics, angles, {});
}

std::optional<Expr> differentiateSeriesExpression(
    const SeriesData& series,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (series.variable != variable || series.exponentDenominator == 0)
        return std::nullopt;

    const std::int64_t denominator = static_cast<std::int64_t>(series.exponentDenominator);
    const std::size_t coefficientCount = series.coefficients.size();
    const bool atInfinity = isPositiveInfinityCenter(series.center);
    const std::int64_t exponentShift = atInfinity ? denominator : -denominator;
    SeriesData result{
        series.variable,
        series.center,
        std::vector<Expr>(coefficientCount, integer(0)),
        series.minimumExponent + exponentShift,
        series.orderNumerator + exponentShift,
        series.exponentDenominator,
        std::vector<std::vector<Expr>>(
            series.logarithmicCoefficients.size(),
            std::vector<Expr>(coefficientCount, integer(0)))};

    const std::size_t maximumLogDegree = series.logarithmicCoefficients.size();
    for (std::size_t i = 0; i < coefficientCount; ++i) {
        const std::int64_t exponentNumerator =
            series.minimumExponent + static_cast<std::int64_t>(i);
        const Expr exponent = rational(exponentNumerator, denominator);
        for (std::size_t logDegree = 0; logDegree <= maximumLogDegree; ++logDegree) {
            const Expr& coefficient = logarithmicCoefficient(series, logDegree, i);
            if (isZero(coefficient)) continue;

            if (exponentNumerator != 0) {
                Expr contribution = multiply(
                    coefficient, exponent,
                    builtins, mathematics, angles, {});
                if (atInfinity)
                    contribution = negate(
                        contribution, builtins, mathematics, angles, {});
                addLogarithmicCoefficient(
                    result, logDegree, i, contribution,
                    builtins, mathematics, angles);
            }
            if (logDegree != 0) {
                Expr contribution = multiply(
                    coefficient, integer(static_cast<std::int64_t>(logDegree)),
                    builtins, mathematics, angles, {});
                if (atInfinity)
                    contribution = negate(
                        contribution, builtins, mathematics, angles, {});
                addLogarithmicCoefficient(
                    result, logDegree - 1, i, contribution,
                    builtins, mathematics, angles);
            }
        }
    }

    trimSeriesDataLeadingZeros(result);
    return makeSeriesData(std::move(result), builtins);
}

std::optional<Expr> integrateSeriesExpression(
    const SeriesData& series,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (series.variable != variable || series.exponentDenominator == 0)
        return std::nullopt;

    const std::int64_t denominator = static_cast<std::int64_t>(series.exponentDenominator);
    const std::size_t coefficientCount = series.coefficients.size();
    const std::size_t inputMaximumLogDegree = series.logarithmicCoefficients.size();
    const bool atInfinity = isPositiveInfinityCenter(series.center);
    // O(t)をxで積分するとlog(t)型の剰余を生み得るため，現在のO(t^r)表現では閉じない。
    if (atInfinity && series.orderNumerator == denominator)
        return std::nullopt;
    const std::int64_t exponentShift = atInfinity ? -denominator : denominator;
    SeriesData result{
        series.variable,
        series.center,
        std::vector<Expr>(coefficientCount, integer(0)),
        series.minimumExponent + exponentShift,
        series.orderNumerator + exponentShift,
        series.exponentDenominator,
        std::vector<std::vector<Expr>>(
            inputMaximumLogDegree + 1,
            std::vector<Expr>(coefficientCount, integer(0)))};

    for (std::size_t i = 0; i < coefficientCount; ++i) {
        const std::int64_t exponentNumerator =
            series.minimumExponent + static_cast<std::int64_t>(i);
        for (std::size_t logDegree = 0; logDegree <= inputMaximumLogDegree; ++logDegree) {
            const Expr& coefficient = logarithmicCoefficient(series, logDegree, i);
            if (isZero(coefficient)) continue;

            if (atInfinity) {
                // t=1/xではdx=-t^-2 dt。指数rは積分でr-1へ移る。
                if (exponentNumerator == denominator) {
                    addLogarithmicCoefficient(
                        result, logDegree + 1, i,
                        negate(
                            divide(coefficient,
                                integer(static_cast<std::int64_t>(logDegree + 1)),
                                builtins, mathematics, angles, {}),
                            builtins, mathematics, angles, {}),
                        builtins, mathematics, angles);
                    continue;
                }

                const Expr shiftedExponent = rational(
                    exponentNumerator - denominator, denominator);
                Expr current = negate(
                    divide(coefficient, shiftedExponent,
                        builtins, mathematics, angles, {}),
                    builtins, mathematics, angles, {});
                addLogarithmicCoefficient(
                    result, logDegree, i, current,
                    builtins, mathematics, angles);

                for (std::size_t degree = logDegree; degree != 0; --degree) {
                    current = multiply(
                        current,
                        divide(
                            integer(-static_cast<std::int64_t>(degree)), shiftedExponent,
                            builtins, mathematics, angles, {}),
                        builtins, mathematics, angles, {});
                    addLogarithmicCoefficient(
                        result, degree - 1, i, current,
                        builtins, mathematics, angles);
                }
                continue;
            }

            if (exponentNumerator == -denominator) {
                addLogarithmicCoefficient(
                    result, logDegree + 1, i,
                    divide(coefficient, integer(static_cast<std::int64_t>(logDegree + 1)),
                        builtins, mathematics, angles, {}),
                    builtins, mathematics, angles);
                continue;
            }

            const Expr shiftedExponent = rational(
                exponentNumerator + denominator, denominator);
            Expr current = divide(
                coefficient, shiftedExponent,
                builtins, mathematics, angles, {});
            addLogarithmicCoefficient(
                result, logDegree, i, current,
                builtins, mathematics, angles);

            for (std::size_t degree = logDegree; degree != 0; --degree) {
                current = multiply(
                    current,
                    divide(
                        integer(-static_cast<std::int64_t>(degree)), shiftedExponent,
                        builtins, mathematics, angles, {}),
                    builtins, mathematics, angles, {});
                addLogarithmicCoefficient(
                    result, degree - 1, i, current,
                    builtins, mathematics, angles);
            }
        }
    }

    trimSeriesDataLeadingZeros(result);
    return makeSeriesData(std::move(result), builtins);
}

} // namespace mmcal::symbolic
