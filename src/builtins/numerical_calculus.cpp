// 数値微分・保証付き数値積分
#include "numerical_calculus.hpp"

#include "approximation/approximation_context.hpp"
#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/expression_interval.hpp"
#include "approximation/complex_interval.hpp"
#include "approximation/real_interval.hpp"
#include "error/error_message.hpp"
#include "expression/exact_value.hpp"
#include "evaluation/iterator_spec.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/substitution.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using approximation::CertifiedBinding;
using approximation::CertifiedEvaluator;
using approximation::CertifiedValue;
using approximation::ComplexInterval;
using approximation::RealInterval;
using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

constexpr std::size_t defaultDigits = 16;
constexpr std::size_t maximumSubintervals = 1U << 18;
constexpr std::size_t newtonCotesDegree = 8;

[[nodiscard]] Rational rational(std::int64_t value) {
    return Rational{BigInt{value}};
}

[[nodiscard]] Rational rationalFromSize(std::size_t value) {
    return Rational{BigInt::parse(std::to_string(value))};
}

[[nodiscard]] Rational absRational(const Rational& value) {
    return value.numerator().isNegative() ? -value : value;
}

[[nodiscard]] Rational maxAbs(const RealInterval& interval) {
    const Rational lower = absRational(interval.lower().toRational());
    const Rational upper = absRational(interval.upper().toRational());
    return lower < upper ? upper : lower;
}



[[nodiscard]] Expr simplifyExpression(
    Expr expression,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{registry, mathematics, angles});
}

[[nodiscard]] Expr addExpression(
    Expr lhs, Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplifyExpression(
        Expr::call(registry.symbol(BuiltinId::Add), {std::move(lhs), std::move(rhs)}),
        registry, mathematics, angles);
}

[[nodiscard]] Expr subtractExpression(
    Expr lhs, Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplifyExpression(
        Expr::call(registry.symbol(BuiltinId::Subtract), {std::move(lhs), std::move(rhs)}),
        registry, mathematics, angles);
}

[[nodiscard]] Expr multiplyExpression(
    Expr lhs, Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplifyExpression(
        Expr::call(registry.symbol(BuiltinId::Multiply), {std::move(lhs), std::move(rhs)}),
        registry, mathematics, angles);
}

[[nodiscard]] bool containsDerivativeCall(const Expr& expression, const evaluation::BuiltinRegistry& registry) {
    if (expression.isCall()) {
        if (expression.asCall().head.sameIdentity(registry.symbol(BuiltinId::Derivative)))
            return true;
        for (const Expr& argument : expression.asCall().arguments)
            if (containsDerivativeCall(argument, registry))
                return true;
    }
    if (expression.isArray())
        for (const Expr& element : expression.asArray().storedExpressions())
            if (containsDerivativeCall(element, registry))
                return true;
    return false;
}

[[nodiscard]] CertifiedValue addValue(
    const CertifiedValue& lhs, const CertifiedValue& rhs, std::size_t bits) {
    if (lhs.isReal() && rhs.isReal())
        return CertifiedValue{approximation::add(lhs.asReal(), rhs.asReal(), bits)};
    return CertifiedValue{approximation::add(lhs.toComplex(), rhs.toComplex(), bits)};
}

[[nodiscard]] CertifiedValue scaleValue(
    const CertifiedValue& value, const Rational& factor, std::size_t bits) {
    const RealInterval scalar = RealInterval::fromRational(factor, bits);
    if (value.isReal())
        return CertifiedValue{approximation::multiply(value.asReal(), scalar, bits)};
    return CertifiedValue{approximation::multiply(
        value.asComplex(), ComplexInterval::fromReal(scalar), bits)};
}

[[nodiscard]] RealInterval inflate(const RealInterval& interval, const Rational& error, std::size_t bits) {
    return RealInterval::fromRationalBounds(
        interval.lower().toRational() - error,
        interval.upper().toRational() + error,
        bits);
}

[[nodiscard]] CertifiedValue inflate(
    const CertifiedValue& value,
    const Rational& realError,
    const Rational& imaginaryError,
    std::size_t bits) {
    if (value.isReal())
        return CertifiedValue{inflate(value.asReal(), realError, bits)};
    return CertifiedValue{ComplexInterval{
        inflate(value.asComplex().real(), realError, bits),
        inflate(value.asComplex().imaginary(), imaginaryError, bits)}};
}

[[nodiscard]] std::optional<CertifiedValue> encloseAt(
    const CertifiedEvaluator& evaluator,
    const Expr& expression,
    const expression::Symbol& variable,
    const Rational& point,
    std::size_t bits,
    CertifiedEvaluator::EnclosureKind enclosureKind) {
    const CertifiedBinding binding{
        variable,
        CertifiedValue{RealInterval::fromRational(point, bits)}};
    return evaluator.enclose(
        expression, bits, std::span<const CertifiedBinding>{&binding, 1}, enclosureKind);
}

[[nodiscard]] std::optional<CertifiedValue> encloseOn(
    const CertifiedEvaluator& evaluator,
    const Expr& expression,
    const expression::Symbol& variable,
    const Rational& lower,
    const Rational& upper,
    std::size_t bits,
    CertifiedEvaluator::EnclosureKind enclosureKind) {
    const CertifiedBinding binding{
        variable,
        CertifiedValue{RealInterval::fromRationalBounds(lower, upper, bits)}};
    return evaluator.enclose(
        expression, bits, std::span<const CertifiedBinding>{&binding, 1}, enclosureKind);
}

struct ErrorBound final {
    Rational real{BigInt{0}};
    Rational imaginary{BigInt{0}};
};

[[nodiscard]] std::vector<Rational> multiplyPolynomialByLinear(
    const std::vector<Rational>& polynomial,
    const Rational& constant,
    const Rational& linear) {
    std::vector<Rational> result(polynomial.size() + 1, Rational{BigInt{0}});
    for (std::size_t i = 0; i < polynomial.size(); ++i) {
        result[i] += polynomial[i] * constant;
        result[i + 1] += polynomial[i] * linear;
    }
    return result;
}

[[nodiscard]] Rational evaluatePolynomial(
    const std::vector<Rational>& coefficients,
    const Rational& x) {
    Rational result{BigInt{0}};
    for (std::size_t i = coefficients.size(); i-- > 0;)
        result = result * x + coefficients[i];
    return result;
}

[[nodiscard]] std::vector<Rational> antiderivative(
    const std::vector<Rational>& polynomial) {
    std::vector<Rational> result(polynomial.size() + 1, Rational{BigInt{0}});
    for (std::size_t i = 0; i < polynomial.size(); ++i)
        result[i + 1] = polynomial[i] / rationalFromSize(i + 1);
    return result;
}

struct NewtonCotesRule final {
    std::size_t degree = 0;
    std::vector<Rational> weights;
    Rational absoluteRemainderIntegral{BigInt{0}};
};

[[nodiscard]] NewtonCotesRule buildNewtonCotesRule(std::size_t degree) {
    NewtonCotesRule rule;
    rule.degree = degree;
    rule.weights.reserve(degree + 1);
    const Rational endpoint = rationalFromSize(degree);

    // Lagrange基底をexact Rational多項式として構築し、[0,m]で積分する。
    for (std::size_t i = 0; i <= degree; ++i) {
        std::vector<Rational> basis{Rational{BigInt{1}}};
        Rational denominator{BigInt{1}};
        for (std::size_t j = 0; j <= degree; ++j) {
            if (j == i)
                continue;
            basis = multiplyPolynomialByLinear(
                basis, -rationalFromSize(j), Rational{BigInt{1}});
            const BigInt difference = BigInt::parse(std::to_string(i))
                - BigInt::parse(std::to_string(j));
            denominator *= Rational{difference};
        }
        for (Rational& coefficient : basis)
            coefficient /= denominator;
        const auto integral = antiderivative(basis);
        rule.weights.push_back(evaluatePolynomial(integral, endpoint));
    }

    // 補間剰余の product(t-i) は各整数区間で符号一定。Integrate|product| dt を区間ごとにexactに積分して誤差定数を得る。
    std::vector<Rational> nodal{Rational{BigInt{1}}};
    for (std::size_t j = 0; j <= degree; ++j)
        nodal = multiplyPolynomialByLinear(
            nodal, -rationalFromSize(j), Rational{BigInt{1}});
    const auto nodalIntegral = antiderivative(nodal);
    Rational absoluteIntegral{BigInt{0}};
    for (std::size_t j = 0; j < degree; ++j) {
        const Rational a = rationalFromSize(j);
        const Rational b = rationalFromSize(j + 1);
        absoluteIntegral += absRational(
            evaluatePolynomial(nodalIntegral, b) - evaluatePolynomial(nodalIntegral, a));
    }
    rule.absoluteRemainderIntegral = std::move(absoluteIntegral);
    return rule;
}

[[nodiscard]] const NewtonCotesRule& integrationRule() {
    static const NewtonCotesRule rule = buildNewtonCotesRule(newtonCotesDegree);
    return rule;
}

[[nodiscard]] ErrorBound compositeNewtonCotesErrorBound(
    const CertifiedValue& derivativeRange,
    const Rational& totalWidth,
    std::size_t subintervals) {
    const NewtonCotesRule& rule = integrationRule();
    const Rational n = rationalFromSize(subintervals);
    const Rational h = absRational(totalWidth) / n;
    const Rational panelCount = n / rationalFromSize(rule.degree);
    const Rational coefficient = panelCount
        * numeric::pow(h, static_cast<std::uint64_t>(rule.degree + 2))
        * rule.absoluteRemainderIntegral
        / Rational{numeric::factorial(static_cast<std::uint64_t>(rule.degree + 1))};

    ErrorBound result;
    if (derivativeRange.isReal()) {
        result.real = coefficient * maxAbs(derivativeRange.asReal());
        return result;
    }
    result.real = coefficient * maxAbs(derivativeRange.asComplex().real());
    result.imaginary = coefficient * maxAbs(derivativeRange.asComplex().imaginary());
    return result;
}

[[nodiscard]] std::optional<CertifiedValue> integrateWithPanels(
    const CertifiedEvaluator& evaluator,
    const Expr& integrand,
    const expression::Symbol& variable,
    const Rational& lower,
    const Rational& upper,
    std::size_t subintervals,
    const ErrorBound& errorBound,
    std::size_t bits,
    CertifiedEvaluator::EnclosureKind enclosureKind) {
    const NewtonCotesRule& rule = integrationRule();
    const Rational h = (upper - lower) / rationalFromSize(subintervals);
    CertifiedValue total{RealInterval::fromRational(rational(0), bits)};
    for (std::size_t panel = 0; panel < subintervals; panel += rule.degree) {
        CertifiedValue panelValue{RealInterval::fromRational(rational(0), bits)};
        const Rational a = lower + h * rationalFromSize(panel);
        for (std::size_t node = 0; node <= rule.degree; ++node) {
            const Rational x = a + h * rationalFromSize(node);
            const auto value = encloseAt(
                evaluator, integrand, variable, x, bits, enclosureKind);
            if (!value)
                return std::nullopt;
            panelValue = addValue(
                panelValue,
                scaleValue(*value, rule.weights[node], bits),
                bits);
        }
        total = addValue(total, scaleValue(panelValue, h, bits), bits);
    }
    return inflate(total, errorBound.real, errorBound.imaginary, bits);
}

[[nodiscard]] Rational decimalTolerance(std::size_t digits) {
    std::string denominator = "1";
    denominator.append(digits + 2, '0');
    return Rational{BigInt{1}, BigInt::parse(denominator)};
}

[[nodiscard]] std::size_t nextGuardDigits(std::size_t current) {
    const std::size_t growth = std::max<std::size_t>(8, current / 2);
    if (growth > std::numeric_limits<std::size_t>::max() - current)
        error::throwCalcError(error::CalcErrorType::Overflow, "Numerical precision is too large");
    return current + growth;
}

[[nodiscard]] bool exceedsLocalGuardBudget(
    std::size_t guardDigits,
    std::size_t requestedDigits) noexcept {
    constexpr std::size_t extraGuardDigits = 256;
    return guardDigits > extraGuardDigits
        && guardDigits - extraGuardDigits > requestedDigits;
}

} // namespace

[[nodiscard]] std::size_t requestedSignificantDigits(const Expr& value) {
    if (value.isDecimalApproximation())
        return value.asDecimalApproximation().requestedSignificantDigits();
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return std::min(
            complex.real().requestedSignificantDigits(),
            complex.imaginary().requestedSignificantDigits());
    }
    return 0;
}

[[nodiscard]] std::optional<Expr> finalizeNumericalCalculusApproximation(
    const CertifiedValue& certified,
    const CertifiedValue& information,
    std::size_t fractionalDigits) {
    // まずNと同じ有効桁評価でfinite-input由来の上限を検査する。
    // 情報量が要求値を下回る場合はその表示を採用し，十分な情報がある場合だけ
    // diff/nintegrate従来契約の「小数部p桁」で確定する。
    const auto limited = approximation::finalizeCertifiedApproximation(
        certified, information, fractionalDigits);
    if (!limited)
        return std::nullopt;
    const std::size_t available = requestedSignificantDigits(*limited);
    if (available != 0 && available < fractionalDigits)
        return limited;
    return approximation::finalizeCertifiedApproximationFixed(
        certified, information, fractionalDigits);
}

Expr evaluateNumericDerivative(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() < 3 || arguments.size() > 4 || !arguments[1].isSymbol())
        error::throwCalcError(
            error::CalcErrorType::Type,
            "diff expects diff[expression, variable, point] with optional precision");
    const std::size_t digits = arguments.size() == 4
        ? expression::exact::positiveSize(arguments[3]).value_or(0)
        : defaultDigits;
    if (digits == 0)
        error::throwCalcError(error::CalcErrorType::Type, "diff precision must be a positive integer");
    evaluation::checkEvaluationRequestedPrecisionDigits(digits);

    const expression::Symbol variable = arguments[1].asSymbol();
    const Expr derivative = symbolic::differentiateExpression(
        arguments[0], variable, registry, mathematics, angles);
    if (containsDerivativeCall(derivative, registry))
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            "diff cannot numerically evaluate an unresolved symbolic derivative");

    const Expr pointExpression = arguments[2];
    const Expr atPoint = symbolic::substituteSymbol(derivative, variable, pointExpression);
    CertifiedEvaluator evaluator{registry, mathematics, angles};
    approximation::ApproximationContext context{digits};
    for (;;) {
        evaluation::consumeEvaluationBudget(evaluation::EvaluationResource::CertifiedRefinement);
        try {
            const std::size_t bits = context.workingBinaryBits();
            const auto information = evaluator.enclose(
                atPoint, bits, CertifiedEvaluator::EnclosureKind::Information);
            const auto enclosed = evaluator.enclose(
                atPoint, bits, CertifiedEvaluator::EnclosureKind::Certified);
            if (!information || !enclosed)
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "diff point or derivative is not numerically evaluable");
            if (const auto decimal = finalizeNumericalCalculusApproximation(
                    *enclosed, *information, digits))
                return *decimal;
        }
        catch (const approximation::PrecisionInsufficient& exception) {
            if (!exception.refinable())
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "diff cannot resolve finite-precision input information");
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            error::throwCalcError(
                error::CalcErrorType::Evaluation,
                "diff encountered an expression unsupported by certified evaluation");
        }
        context.setGuardDigits(nextGuardDigits(context.guardDigits()));
    }
}

Expr evaluateNumericIntegral(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments.size() < 2 || arguments.size() > 3)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "nintegrate expects nintegrate[expression, {variable, lower, upper}] with optional precision");

    const auto iterator = evaluation::parseRangeIteratorSpec(arguments[1]);
    if (!iterator)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "nintegrate iterator must be {variable, lower, upper}");

    const std::size_t digits = arguments.size() == 3
        ? expression::exact::positiveSize(arguments[2]).value_or(0)
        : defaultDigits;
    if (digits == 0)
        error::throwCalcError(
            error::CalcErrorType::Type,
            "nintegrate precision must be a positive integer");
    evaluation::checkEvaluationRequestedPrecisionDigits(digits);

    const expression::Symbol variable = iterator->variable;
    const Expr lowerExpression = iterator->lower;
    const Expr upperExpression = iterator->upper;

    // 任意の、実数であることが証明可能な有限の境界は、[0,1] 内のtに正規化される：
    //   x = a + (b-a)t,  dx = (b-a)dt。
    // これにより、Pi、E、根号などを使用しながらも、求積メッシュを有理数に保つことができる。
    // 逆の境界は、(b-a)の符号によって自動的に処理される。

    CertifiedEvaluator evaluator{registry, mathematics, angles};
    approximation::ApproximationContext context{digits};
    for (;;) {
        evaluation::consumeEvaluationBudget(evaluation::EvaluationResource::CertifiedRefinement);
        const std::size_t bits = context.workingBinaryBits();
        try {
            const auto lowerValue = evaluator.enclose(lowerExpression, bits);
            const auto upperValue = evaluator.enclose(upperExpression, bits);
            if (!lowerValue || !upperValue)
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "nintegrate bounds are not numerically evaluable");
            if (!lowerValue->isReal() || !upperValue->isReal())
                error::throwCalcError(
                    error::CalcErrorType::Type,
                    "nintegrate requires finite real bounds");
            break;
        }
        catch (const approximation::PrecisionInsufficient&) {
            context.setGuardDigits(nextGuardDigits(context.guardDigits()));
            if (exceedsLocalGuardBudget(context.guardDigits(), digits))
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "nintegrate could not certify real finite bounds");
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            error::throwCalcError(
                error::CalcErrorType::Evaluation,
                "nintegrate bounds are unsupported by certified evaluation");
        }
        catch (const std::domain_error&) {
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "nintegrate requires finite real bounds");
        }
    }

    const Expr delta = subtractExpression(
        upperExpression, lowerExpression, registry, mathematics, angles);
    if (delta.isNumber() && delta.asNumber().isZero())
        return Expr{Number{BigInt{0}}};

    // Internal-only symbol. Identity, rather than spelling, is used by differentiation/binding.
    const expression::Symbol parameter{"$__mmcal_nintegrate_parameter"};
    const Expr mappedPoint = addExpression(
        lowerExpression,
        multiplyExpression(delta, Expr{parameter}, registry, mathematics, angles),
        registry, mathematics, angles);
    const Expr transformedIntegrand = multiplyExpression(
        delta,
        symbolic::substituteSymbol(arguments[0], variable, mappedPoint),
        registry, mathematics, angles);

    // 高階導函数を構築する前に、元の被積分函数を区間全体で一度評価する。
    // 1/x のように区間内に明白な特異点がある場合、導函数ASTを巨大化させる前にDomainErrorへ落とす。
    // PrecisionInsufficient/未対応はここでは決め打ちせず、後段の導函数上界評価へ委ねる。
    try {
        static_cast<void>(encloseOn(
            evaluator, transformedIntegrand, parameter,
            Rational{BigInt{0}}, Rational{BigInt{1}},
            context.workingBinaryBits(),
            CertifiedEvaluator::EnclosureKind::Information));
    }
    catch (const approximation::PrecisionInsufficient&) {
    }
    catch (const approximation::CertifiedBackendUnsupported&) {
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            "nintegrate encountered an expression unsupported by certified evaluation");
    }
    catch (const std::domain_error&) {
        error::throwCalcError(
            error::CalcErrorType::Domain,
            "nintegrate could not certify the integrand over the interval (possible singularity)");
    }

    Expr errorDerivative = transformedIntegrand;
    for (std::size_t order = 0; order < integrationRule().degree + 1; ++order)
        errorDerivative = symbolic::differentiateExpression(
            errorDerivative, parameter, registry, mathematics, angles);
    if (containsDerivativeCall(errorDerivative, registry))
        error::throwCalcError(
            error::CalcErrorType::Evaluation,
            "nintegrate requires a certifiable derivative required by its Newton-Cotes error bound");

    // precision誤差とquadrature誤差を分離する。
    // 高階導函数を区間全体で一度だけ囲い、補間剰余のexact Rational上界が十分小さくなるsubinterval数を先に決めてから函数値を採る。
    // これによりnを増やすたび全点を再評価する無駄を避ける。
    const Rational tolerance = decimalTolerance(digits);
    for (;;) {
        evaluation::consumeEvaluationBudget(evaluation::EvaluationResource::CertifiedRefinement);
        const std::size_t bits = context.workingBinaryBits();
        try {
            const Rational lower{BigInt{0}};
            const Rational upper{BigInt{1}};
            const auto informationDerivativeRange = encloseOn(
                evaluator, errorDerivative, parameter, lower, upper, bits,
                CertifiedEvaluator::EnclosureKind::Information);
            const auto derivativeRange = encloseOn(
                evaluator, errorDerivative, parameter, lower, upper, bits,
                CertifiedEvaluator::EnclosureKind::Certified);
            if (!informationDerivativeRange || !derivativeRange)
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "nintegrate cannot certify the required derivative over the interval");

            std::size_t n = integrationRule().degree;
            ErrorBound bound = compositeNewtonCotesErrorBound(*derivativeRange, upper - lower, n);
            while ((bound.real > tolerance || bound.imaginary > tolerance)
                && n < maximumSubintervals) {
                n *= 2;
                bound = compositeNewtonCotesErrorBound(*derivativeRange, upper - lower, n);
            }
            if (bound.real > tolerance || bound.imaginary > tolerance)
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "nintegrate requires too many Newton-Cotes subintervals for the requested precision");

            const ErrorBound informationBound = compositeNewtonCotesErrorBound(
                *informationDerivativeRange, upper - lower, n);
            auto information = integrateWithPanels(
                evaluator, transformedIntegrand, parameter, lower, upper, n,
                informationBound, bits, CertifiedEvaluator::EnclosureKind::Information);
            auto enclosed = integrateWithPanels(
                evaluator, transformedIntegrand, parameter, lower, upper, n,
                bound, bits, CertifiedEvaluator::EnclosureKind::Certified);
            if (!information || !enclosed)
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "nintegrate encountered an expression unsupported by certified evaluation");
            if (const auto decimal = finalizeNumericalCalculusApproximation(
                    *enclosed, *information, digits))
                return *decimal;
        }
        catch (const approximation::PrecisionInsufficient& exception) {
            if (!exception.refinable())
                error::throwCalcError(
                    error::CalcErrorType::Evaluation,
                    "nintegrate cannot resolve finite-precision input information");
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            error::throwCalcError(
                error::CalcErrorType::Evaluation,
                "nintegrate encountered an expression unsupported by certified evaluation");
        }
        catch (const std::domain_error&) {
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "nintegrate could not certify the integrand over the interval (possible singularity)");
        }

        context.setGuardDigits(nextGuardDigits(context.guardDigits()));
        if (exceedsLocalGuardBudget(context.guardDigits(), digits))
            error::throwCalcError(
                error::CalcErrorType::Evaluation,
                "nintegrate could not certify the requested decimal precision");
    }
}

} // namespace mmcal::builtins
