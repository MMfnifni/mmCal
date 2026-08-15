// builtin dispatchを非再帰評価機械から分離する。
#include "evaluator.hpp"

#include "mathematics/assumption_parser.hpp"

#include "builtins/arithmetic.hpp"
#include "builtins/approximation_utilities.hpp"
#include "builtins/aggregate.hpp"
#include "builtins/statistics.hpp"
#include "builtins/signal_processing.hpp"
#include "builtins/array.hpp"
#include "builtins/array_vector.hpp"
#include "builtins/comparison.hpp"
#include "builtins/combinatorics.hpp"
#include "builtins/discrete_math.hpp"
#include "builtins/elementary_utilities.hpp"
#include "builtins/complex_functions.hpp"
#include "builtins/hyperbolic.hpp"
#include "builtins/linear_algebra.hpp"
#include "builtins/numerical_calculus.hpp"
#include "builtins/trigonometric.hpp"
#include "builtins/stable_elementary.hpp"
#include "builtins/special_functions.hpp"
#include "builtins/random_functions.hpp"
#include "builtins/transcendental.hpp"
#include "error/error_message.hpp"
#include "evaluation/iterator_spec.hpp"
#include "simplification/full_simplifier.hpp"
#include "simplification/simplifier.hpp"
#include "solver/polynomial_solver.hpp"
#include "solver/solve_constraints.hpp"
#include "solver/transcendental_solver.hpp"
#include "symbolic/algebra_transforms.hpp"
#include "symbolic/differentiation.hpp"
#include "symbolic/integration.hpp"
#include "symbolic/limit.hpp"
#include "numeric/integer_algorithms.hpp"

#include <algorithm>
#include <optional>
#include <span>
#include <vector>

namespace mmcal::evaluation {
namespace {

[[nodiscard]] bool requiresRectangularArray(BuiltinId id) noexcept {
    switch (id) {
    case BuiltinId::Transpose:
    case BuiltinId::ConjugateTranspose:
    case BuiltinId::MatrixAdd:
    case BuiltinId::MatrixMultiply:
    case BuiltinId::Determinant:
    case BuiltinId::Inverse:
    case BuiltinId::Rref:
    case BuiltinId::Rank:
    case BuiltinId::SolveLinear:
    case BuiltinId::NullSpace:
    case BuiltinId::LuDecomposition:
    case BuiltinId::QrDecomposition:
    case BuiltinId::SingularValueDecomposition:
    case BuiltinId::Eigenvalues:
    case BuiltinId::Eigenvectors:
    case BuiltinId::Eigensystem:
    case BuiltinId::Trace:
    case BuiltinId::Rows:
    case BuiltinId::Cols:
    case BuiltinId::Diag:
    case BuiltinId::VectorAdd:
    case BuiltinId::VectorSubtract:
    case BuiltinId::VectorScale:
    case BuiltinId::VectorDot:
    case BuiltinId::VectorCross:
    case BuiltinId::VectorNorm:
    case BuiltinId::VectorManhattan:
    case BuiltinId::VectorEuclidean:
    case BuiltinId::VectorNormalize:
    case BuiltinId::VectorProject:
    case BuiltinId::VectorAngle:
    case BuiltinId::VectorReflect:
    case BuiltinId::VectorReflectAxis:
    case BuiltinId::VectorSum:
        return true;
    default:
        return false;
    }
}

[[nodiscard]] bool containsBuiltinCall(
    const expression::Expr& root,
    const expression::Symbol& head) {
    std::vector<expression::Expr> pending{root};
    while (!pending.empty()) {
        expression::Expr current = std::move(pending.back());
        pending.pop_back();
        if (current.isCall()) {
            const auto& call = current.asCall();
            if (call.head.sameIdentity(head))
                return true;
            for (const expression::Expr& argument : call.arguments)
                pending.push_back(argument);
        }
        else if (current.isArray()) {
            for (const expression::Expr& element : current.asArray().storedExpressions())
                pending.push_back(element);
        }
        else if (current.isList()) {
            for (const expression::Expr& element : current.asList().elements)
                pending.push_back(element);
        }
    }
    return false;
}

[[nodiscard]] bool containsUnresolvedSolution(const solver::SolutionSet& solutions) {
    if (solutions.kind() == solver::SolutionSetKind::Unresolved)
        return true;
    if (solutions.kind() != solver::SolutionSetKind::Conditional)
        return false;
    return std::any_of(solutions.cases().begin(), solutions.cases().end(),
        [](const solver::SolutionCase& item) {
            return item.outcome == solver::SolutionSetKind::Unresolved;
        });
}

} // namespace

expression::Expr Evaluator::dispatchBuiltin(
    const BuiltinDefinition& definition,
    const expression::CallExpr& call,
    std::span<const expression::Expr> arguments) {
    if (requiresRectangularArray(definition.id)) {
        const bool nonRectangular = std::any_of(arguments.begin(), arguments.end(),
            [](const expression::Expr& value) { return value.isList(); });
        if (nonRectangular) {
            emitWarning("Array::nonRectangular",
                std::string{definition.name()}
                    + " requires a rectangular dense array; the brace value remains unevaluated");
            return expression::Expr::call(call.head,
                std::vector<expression::Expr>{arguments.begin(), arguments.end()});
        }
    }

    switch (definition.id) {
    case BuiltinId::Add:
        return builtins::evaluateAdd(arguments, registry_);
    case BuiltinId::Subtract:
        return builtins::evaluateSubtract(arguments, registry_);
    case BuiltinId::Multiply:
        return builtins::evaluateMultiply(arguments, registry_);
    case BuiltinId::Divide:
        return builtins::evaluateDivide(arguments, registry_);
    case BuiltinId::Power:
        return builtins::evaluatePower(arguments, registry_, mathematics_);
    case BuiltinId::Negate:
        return builtins::evaluateNegate(arguments, registry_);
    case BuiltinId::Factorial:
        return builtins::evaluateFactorial(arguments);
    case BuiltinId::Derivative: {
        if (arguments.size() < 2)
            error::throwCalcError(error::CalcErrorType::Type,
                "D expects an expression followed by one or more derivative specifications");

        expression::Expr result = arguments[0];
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            expression::Symbol variable;
            std::uint64_t order = 1;
            if (arguments[i].isSymbol()) {
                variable = arguments[i].asSymbol();
            }
            else if (arguments[i].isArray()) {
                const auto& spec = arguments[i].asArray();
                const expression::Expr specVariable = spec.size() > 0 ? spec.element(0) : arguments[i];
                const expression::Expr specOrder = spec.size() > 1 ? spec.element(1) : arguments[i];
                if (spec.shape.size() != 1 || spec.shape[0] != 2 || spec.size() != 2
                    || !specVariable.isSymbol() || !specOrder.isNumber()
                    || !specOrder.asNumber().isReal()
                    || !specOrder.asNumber().asReal().isInteger()) {
                    error::throwCalcError(error::CalcErrorType::Type,
                        "D derivative specification must be a symbol or {symbol, nonnegative integer}");
                }
                const auto parsed = numeric::tryToUint64(
                    specOrder.asNumber().asReal().asInteger());
                if (!parsed)
                    error::throwCalcError(error::CalcErrorType::Domain,
                        "D derivative order must be a nonnegative integer that fits in uint64");
                if (*parsed > 4096)
                    error::throwCalcError(error::CalcErrorType::Overflow,
                        "D derivative order is too large");
                variable = specVariable.asSymbol();
                order = *parsed;
            }
            else {
                error::throwCalcError(error::CalcErrorType::Type,
                    "D derivative specification must be a symbol or {symbol, nonnegative integer}");
            }

            for (std::uint64_t derivative = 0; derivative < order; ++derivative)
                result = symbolic::differentiateExpression(
                    result, variable, registry_, mathematics_, angleSemantics_);
        }

        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Derivative)))
            emitWarning("D::unevaluated",
                "D could not fully evaluate the derivative; unevaluated D[...] remains");
        return result;
    }
    case BuiltinId::SymbolicIntegral: {
        if (arguments.size() < 2 || arguments.size() > 3)
            error::throwCalcError(error::CalcErrorType::Type,
                "integrate expects integrate[expression, variable, optional assumptions] or integrate[expression, {variable, lower, upper}, optional assumptions]");

        mathematics::AssumptionSet assumptions;
        if (arguments.size() == 3)
            assumptions = mathematics::parseAssumptions(arguments[2], registry_, mathematics_);

        std::optional<expression::Expr> result;
        std::optional<symbolic::IntegrationDisposition> integrationDisposition;
        if (arguments[1].isSymbol()) {
            symbolic::IntegrationResult detailed = symbolic::integrateExpressionDetailed(
                arguments[0], arguments[1].asSymbol(), registry_, mathematics_,
                angleSemantics_, assumptions);
            integrationDisposition = detailed.disposition;
            result = std::move(detailed.expression);
        }
        else if (const auto iterator = parseRangeIteratorSpec(arguments[1])) {
            const auto* infinity = symbolRegistry_.find("Infinity");
            if (!infinity)
                error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
            result = symbolic::integrateExpression(
                arguments[0], iterator->variable, iterator->lower, iterator->upper,
                registry_, mathematics_, angleSemantics_, infinity->symbol, assumptions);
        }
        else {
            error::throwCalcError(error::CalcErrorType::Type,
                "integrate expects a symbol or {variable, lower, upper} as the second argument");
        }

        if (containsBuiltinCall(*result, registry_.symbol(BuiltinId::SymbolicIntegral))) {
            if (!integrationDisposition) {
                emitWarning("integrate::conditionsRequired",
                    "integrate kept the definite integral unevaluated because a safe symbolic result could not be established on the requested interval");
            }
            else {
                switch (*integrationDisposition) {
                case symbolic::IntegrationDisposition::Partial:
                    emitWarning("integrate::partial",
                        "integrate partially evaluated the expression; remaining subintegral(s) are outside the current symbolic rule set");
                    break;
                case symbolic::IntegrationDisposition::KnownNoFiniteClosedForm:
                    emitWarning("integrate::noKnownClosedForm",
                        "integrate recognized a family with no known finite closed form in mmCal's supported standard-function vocabulary; the integral remains unevaluated");
                    break;
                case symbolic::IntegrationDisposition::ConditionsRequired:
                    emitWarning("integrate::conditionsRequired",
                        "integrate needs additional domain or branch assumptions before it can choose a safe symbolic antiderivative");
                    break;
                case symbolic::IntegrationDisposition::UnsupportedByEngine:
                    emitWarning("integrate::unsupported",
                        "mmCal has no implemented symbolic integration rule for this expression; this does not imply that no closed form exists");
                    break;
                case symbolic::IntegrationDisposition::Solved:
                    emitWarning("integrate::unsupported",
                        "integrate left an unexpected unevaluated subintegral; this does not imply that no closed form exists");
                    break;
                }
            }
        }
        return *result;
    }
    case BuiltinId::Limit: {
        if (arguments.size() < 3 || arguments.size() > 4 || !arguments[1].isSymbol())
            error::throwCalcError(error::CalcErrorType::Type,
                "limit expects limit[expression, variable, point] or limit[expression, variable, point, direction]");

        symbolic::LimitDirection direction = symbolic::LimitDirection::TwoSided;
        if (arguments.size() == 4) {
            if (!arguments[3].isNumber() || !arguments[3].asNumber().isReal()
                || !arguments[3].asNumber().asReal().isInteger())
                error::throwCalcError(error::CalcErrorType::Type,
                    "limit direction must be -1 for left or 1 for right");
            const auto& value = arguments[3].asNumber().asReal().asInteger();
            if (value == numeric::BigInt{-1})
                direction = symbolic::LimitDirection::Left;
            else if (value == numeric::BigInt{1})
                direction = symbolic::LimitDirection::Right;
            else
                error::throwCalcError(error::CalcErrorType::Domain,
                    "limit direction must be -1 for left or 1 for right");
        }

        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
        expression::Expr result = symbolic::limitExpression(
            arguments[0], arguments[1].asSymbol(), arguments[2], direction,
            registry_, mathematics_, angleSemantics_, infinity->symbol);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Limit)))
            emitWarning("limit::unevaluated",
                "limit could not prove the requested limit; unevaluated limit[...] remains");
        return result;
    }
    case BuiltinId::Floor:
        return builtins::evaluateFloor(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Ceil:
        return builtins::evaluateCeil(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Trunc:
        return builtins::evaluateTrunc(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Round:
        return builtins::evaluateRound(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Frac:
        return builtins::evaluateFrac(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Gcd:
        return builtins::evaluateGcd(arguments, registry_);
    case BuiltinId::Lcm:
        return builtins::evaluateLcm(arguments, registry_);
    case BuiltinId::Mod:
        return builtins::evaluateMod(arguments, registry_);
    case BuiltinId::Rem:
        return builtins::evaluateRem(arguments, registry_);
    case BuiltinId::Quotient:
        return builtins::evaluateQuotient(arguments, registry_);
    case BuiltinId::Permutation:
        return builtins::evaluatePermutation(arguments, registry_);
    case BuiltinId::Combination:
        return builtins::evaluateCombination(arguments, registry_);
    case BuiltinId::Fibonacci:
        return builtins::evaluateFibonacci(arguments, registry_);
    case BuiltinId::DiscreteFourierTransform:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateDft(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateDft(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::FastFourierTransform:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateFft(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateFft(arguments, registry_, mathematics_, angleSemantics_, fourierTransformCache_);
    case BuiltinId::InverseFourierTransform:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateIfft(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateIfft(arguments, registry_, mathematics_, angleSemantics_, fourierTransformCache_);
    case BuiltinId::Convolution:
        return builtins::evaluateConvolution(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Transpose:
        return builtins::evaluateTranspose(arguments, registry_);
    case BuiltinId::ConjugateTranspose:
        return builtins::evaluateConjugateTranspose(
            arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::MatrixAdd:
        return builtins::evaluateMatrixAdd(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::MatrixMultiply:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateDot(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateDot(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Determinant:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateDeterminant(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateDeterminant(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Inverse:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateInverse(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateMatrixInverse(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Rref: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateRref(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateRref(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Rref)))
            emitWarning("rref::unevaluated",
                "rref could not determine the symbolic pivots; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::Rank: {
        // rankは不連続量なので、exact入力ではまずexact eliminationを優先する。
        // exactに決まらない場合だけcertified intervalへ降ろし、epsilon判定は導入しない。
        expression::Expr result = builtins::evaluateMatrixRank(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Rank))) {
            if (const auto* approximation = currentApproximationContext())
                if (const auto approximateResult = builtins::evaluateApproximateMatrixRank(
                    arguments, registry_, mathematics_, angleSemantics_, *approximation))
                    return *approximateResult;
            emitWarning("matrixRank::unevaluated",
                "matrixRank could not determine the symbolic pivots; the expression remains unevaluated");
        }
        return result;
    }
    case BuiltinId::SolveLinear: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateSolveLinear(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateSolveLinear(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::SolveLinear)))
            emitWarning("solveLinear::unevaluated",
                "solveLinear could not determine the symbolic pivots; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::NullSpace: {
        // nullSpaceもrankと同じくrank deficiencyに依存する不連続量なので，
        // exact入力ではまずexact eliminationを優先する。未解決時だけcertified intervalへ降ろす。
        expression::Expr result = builtins::evaluateNullSpace(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::NullSpace))) {
            if (const auto* approximation = currentApproximationContext())
                if (const auto approximateResult = builtins::evaluateApproximateNullSpace(
                    arguments, registry_, mathematics_, angleSemantics_, *approximation))
                    return *approximateResult;
            emitWarning("nullSpace::unevaluated",
                "nullSpace could not certify the pivot structure; the expression remains unevaluated");
        }
        return result;
    }
    case BuiltinId::LuDecomposition: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateLuDecomposition(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateLuDecomposition(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::LuDecomposition)))
            emitWarning("luDecomposition::unevaluated",
                "luDecomposition could not prove a required symbolic pivot nonzero; the expression remains unevaluated");
        return result;
    }
    case BuiltinId::QrDecomposition: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateQrDecomposition(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateQrDecomposition(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::QrDecomposition)))
            emitWarning("qrDecomposition::unevaluated",
                "qrDecomposition exact Householder expansion is unavailable for this matrix; use N[...] for the certified numerical backend");
        return result;
    }
    case BuiltinId::SingularValueDecomposition: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateSingularValueDecomposition(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateSingularValueDecomposition(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::SingularValueDecomposition)))
            emitWarning("svd::unevaluated",
                "svd exact form is only emitted for natural exact cases; use N[...] for the numerical SVD backend");
        return result;
    }
    case BuiltinId::Eigenvalues: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateEigenvalues(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateEigenvalues(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Eigenvalues)))
            emitWarning("eigenvalues::unevaluated",
                "eigenvalues exact form is unavailable for this matrix; use N[...] for the numerical Schur backend");
        return result;
    }
    case BuiltinId::Eigenvectors: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateEigenvectors(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateEigenvectors(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Eigenvectors)))
            emitWarning("eigenvectors::unevaluated",
                "eigenvectors exact form is only emitted for natural exact cases; use N[...] for the numerical Schur backend");
        return result;
    }
    case BuiltinId::Eigensystem: {
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateEigensystem(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        expression::Expr result = builtins::evaluateEigensystem(
            arguments, registry_, mathematics_, angleSemantics_);
        if (containsBuiltinCall(result, registry_.symbol(BuiltinId::Eigensystem)))
            emitWarning("eigensystem::unevaluated",
                "eigensystem exact form is only emitted for natural exact cases; use N[...] for the numerical Schur backend");
        return result;
    }
    case BuiltinId::NumericDerivative:
        return builtins::evaluateNumericDerivative(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::NumericIntegral:
        return builtins::evaluateNumericIntegral(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Cbrt:
        return builtins::evaluateCbrt(arguments, registry_, mathematics_);
    case BuiltinId::Hypot:
        return builtins::evaluateHypot(arguments, registry_);
    case BuiltinId::Cis:
        return builtins::evaluateCis(arguments, registry_);
    case BuiltinId::Polar:
        return builtins::evaluatePolar(arguments, registry_);
    case BuiltinId::NextPow2:
        return builtins::evaluateNextPow2(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::DegreeToRadian:
        return builtins::evaluateDegreeToRadian(arguments, registry_, mathematics_);
    case BuiltinId::DegreeToGradian:
        return builtins::evaluateDegreeToGradian(arguments, registry_);
    case BuiltinId::RadianToDegree:
        return builtins::evaluateRadianToDegree(arguments, registry_, mathematics_);
    case BuiltinId::RadianToGradian:
        return builtins::evaluateRadianToGradian(arguments, registry_, mathematics_);
    case BuiltinId::GradianToDegree:
        return builtins::evaluateGradianToDegree(arguments, registry_);
    case BuiltinId::GradianToRadian:
        return builtins::evaluateGradianToRadian(arguments, registry_, mathematics_);
    case BuiltinId::Sum:
        return builtins::evaluateSum(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Product:
        return builtins::evaluateProduct(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Min:
        return builtins::evaluateMin(arguments, registry_);
    case BuiltinId::Max:
        return builtins::evaluateMax(arguments, registry_);
    case BuiltinId::Mean:
        return builtins::evaluateMean(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Median:
    case BuiltinId::Mode:
    case BuiltinId::Quantile:
    case BuiltinId::Percentile:
    case BuiltinId::VariancePopulation:
    case BuiltinId::VarianceSample:
    case BuiltinId::StddevPopulation:
    case BuiltinId::StddevSample:
    case BuiltinId::GeometricMean:
    case BuiltinId::HarmonicMean:
    case BuiltinId::Rms:
    case BuiltinId::MedianAbsoluteDeviation:
    case BuiltinId::MeanAbsoluteDeviation:
    case BuiltinId::Skewness:
    case BuiltinId::KurtosisPopulation:
    case BuiltinId::KurtosisSample:
    case BuiltinId::CoefficientVariation:
    case BuiltinId::StandardError:
    case BuiltinId::ZScore:
    case BuiltinId::Iqr:
    case BuiltinId::TrimMean:
    case BuiltinId::WinsorMean:
    case BuiltinId::Winsorized:
    case BuiltinId::Covariance:
    case BuiltinId::Correlation:
    case BuiltinId::SpearmanCorrelation:
    case BuiltinId::PercentRank:
        return builtins::evaluateStatistic(definition.id, arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Dimensions:
        return builtins::evaluateDimensions(arguments);
    case BuiltinId::ArrayRank:
        return builtins::evaluateArrayRank(arguments);
    case BuiltinId::Length:
        return builtins::evaluateLength(arguments);
    case BuiltinId::ArrayGet:
        return builtins::evaluateArrayGet(arguments);
    case BuiltinId::Reshape:
        return builtins::evaluateReshape(arguments);
    case BuiltinId::Identity:
        return builtins::evaluateIdentity(arguments);
    case BuiltinId::Zeros:
        return builtins::evaluateZeros(arguments);
    case BuiltinId::MatrixGet:
        return builtins::evaluateMatrixGet(arguments);
    case BuiltinId::Trace:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateTrace(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateTrace(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Rows:
        return builtins::evaluateRows(arguments);
    case BuiltinId::Cols:
        return builtins::evaluateCols(arguments);
    case BuiltinId::Diag:
        return builtins::evaluateDiag(arguments);
    case BuiltinId::VectorAdd:
        return builtins::evaluateVectorAdd(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorSubtract:
        return builtins::evaluateVectorSubtract(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorScale:
        return builtins::evaluateVectorScale(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorDot:
        return builtins::evaluateDot(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorCross:
        return builtins::evaluateVectorCross(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorNorm:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateNorm(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateNorm(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorManhattan:
        return builtins::evaluateVectorManhattan(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorEuclidean:
        return builtins::evaluateVectorEuclidean(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorNormalize:
        if (const auto* approximation = currentApproximationContext())
            if (const auto result = builtins::evaluateApproximateNormalize(
                arguments, registry_, mathematics_, angleSemantics_, *approximation))
                return *result;
        return builtins::evaluateNormalize(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorProject:
        return builtins::evaluateVectorProject(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorAngle:
        return builtins::evaluateVectorAngle(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorReflect:
        return builtins::evaluateVectorReflect(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorReflectAxis:
        return builtins::evaluateVectorReflectAxis(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::VectorSum:
        return builtins::evaluateVectorSum(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Expm1:
    case BuiltinId::Log1p:
    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc:
    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc:
        return builtins::evaluateStableElementary(
            definition.id, arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Log2:
    case BuiltinId::Log10:
    case BuiltinId::Gamma:
    case BuiltinId::LogGamma:
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
    case BuiltinId::Hypergeometric1F1:
    case BuiltinId::Hypergeometric2F1:
    case BuiltinId::EllipticF:
    case BuiltinId::EllipticE:
    case BuiltinId::EllipticPi:
    case BuiltinId::ExponentialIntegralEi:
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
    case BuiltinId::LogarithmicIntegralLi:
    case BuiltinId::Polylog:
    case BuiltinId::Beta:
    case BuiltinId::BetaLog:
    case BuiltinId::GeneralizedBinomial:
    case BuiltinId::FallingFactorial:
    case BuiltinId::RisingFactorial:
        return builtins::evaluateSpecialFunction(
            definition.id, arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::RandSeed:
        return builtins::evaluateRandSeed(arguments, randomEngine_);
    case BuiltinId::Rand:
        return builtins::evaluateRand(arguments, randomEngine_);
    case BuiltinId::RandInt:
        return builtins::evaluateRandInt(arguments, randomEngine_);
    case BuiltinId::Choice:
        return builtins::evaluateChoice(arguments, randomEngine_);
    case BuiltinId::RandN:
        return builtins::evaluateRandN(
            arguments, randomEngine_, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Sqrt:
        return builtins::evaluateSqrt(arguments, registry_, mathematics_);
    case BuiltinId::Abs:
        return builtins::evaluateAbs(arguments, registry_, mathematics_);
    case BuiltinId::Sign:
        return builtins::evaluateSign(arguments, registry_, mathematics_);
    case BuiltinId::Re:
        return builtins::evaluateRe(arguments, registry_, mathematics_);
    case BuiltinId::Im:
        return builtins::evaluateIm(arguments, registry_, mathematics_);
    case BuiltinId::Conj:
        return builtins::evaluateConj(arguments, registry_, mathematics_);
    case BuiltinId::Sin:
        return builtins::evaluateSin(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Cos:
        return builtins::evaluateCos(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Tan:
        return builtins::evaluateTan(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Cot:
        return builtins::evaluateCot(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Sec:
        return builtins::evaluateSec(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Csc:
        return builtins::evaluateCsc(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Asin:
        return builtins::evaluateAsin(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Acos:
        return builtins::evaluateAcos(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Atan:
        return builtins::evaluateAtan(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Atan2:
        return builtins::evaluateAtan2(arguments, registry_, mathematics_, angleSemantics_);
    case BuiltinId::Sinh:
        return builtins::evaluateSinh(arguments, registry_, mathematics_);
    case BuiltinId::Cosh:
        return builtins::evaluateCosh(arguments, registry_, mathematics_);
    case BuiltinId::Tanh:
        return builtins::evaluateTanh(arguments, registry_, mathematics_);
    case BuiltinId::Asinh:
        return builtins::evaluateAsinh(arguments, registry_, mathematics_);
    case BuiltinId::Acosh:
        return builtins::evaluateAcosh(arguments, registry_, mathematics_);
    case BuiltinId::Atanh:
        return builtins::evaluateAtanh(arguments, registry_, mathematics_);
    case BuiltinId::Csch:
        return builtins::evaluateCsch(arguments, registry_, mathematics_);
    case BuiltinId::Sech:
        return builtins::evaluateSech(arguments, registry_, mathematics_);
    case BuiltinId::Coth:
        return builtins::evaluateCoth(arguments, registry_, mathematics_);
    case BuiltinId::Arg:
        return builtins::evaluateArg(arguments, registry_, mathematics_);
    case BuiltinId::Log:
        return builtins::evaluateLog(arguments, registry_, mathematics_);
    case BuiltinId::Exp:
        return builtins::evaluateExp(arguments, registry_, mathematics_);
    case BuiltinId::NumericalApproximation:
        error::throwCalcError(error::CalcErrorType::Internal,
            "N must be evaluated through the precision-aware evaluation path");
    case BuiltinId::Precision: {
        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
        if (const auto result = builtins::evaluatePrecision(arguments, infinity->symbol))
            return *result;
        emitWarning("precision::unevaluated",
            "precision could not determine the guaranteed precision; the expression remains unevaluated");
        return expression::Expr::call(call.head, {arguments.front()});
    }
    case BuiltinId::Accuracy: {
        const auto* infinity = symbolRegistry_.find("Infinity");
        if (!infinity)
            error::throwCalcError(error::CalcErrorType::Internal, "Infinity symbol is not registered");
        if (const auto result = builtins::evaluateAccuracy(arguments, infinity->symbol))
            return *result;
        emitWarning("accuracy::unevaluated",
            "accuracy could not determine the guaranteed accuracy; the expression remains unevaluated");
        return expression::Expr::call(call.head, {arguments.front()});
    }
    case BuiltinId::Rationalize:
        if (const auto result = builtins::evaluateRationalize(arguments))
            return *result;
        emitWarning("rationalize::unevaluated",
            "rationalize could not convert part of the expression; it remains unevaluated");
        return expression::Expr::call(call.head, std::vector<expression::Expr>{arguments.begin(), arguments.end()});
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify: {
        if (arguments.empty() || arguments.size() > 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "simplify expects an expression and optional assumptions");
        mathematics::AssumptionSet assumptions;
        if (arguments.size() == 2)
            assumptions = mathematics::parseAssumptions(arguments[1], registry_, mathematics_);
        const simplification::SimplificationContext context{
            registry_, mathematics_, angleSemantics_, std::move(assumptions)};
        if (definition.id == BuiltinId::FullSimplify)
            return simplification::fullSimplify(arguments.front(), context);
        return simplification::Simplifier{}.simplify(arguments.front(), context);
    }
    case BuiltinId::Expand:
        return symbolic::expandExpression(
            arguments.front(), registry_, mathematics_, angleSemantics_);
    case BuiltinId::Factor:
        return symbolic::factorExpression(
            arguments.front(), registry_, mathematics_, angleSemantics_);
    case BuiltinId::Collect: {
        if (arguments.size() != 2)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "collect expects an expression and a symbol or symbol array");
        std::vector<expression::Symbol> variables;
        if (arguments[1].isSymbol())
            variables.push_back(arguments[1].asSymbol());
        else if (arguments[1].isArray() && arguments[1].asArray().rank() == 1) {
            const auto& array = arguments[1].asArray();
            for (std::size_t i = 0; i < array.size(); ++i) {
                const expression::Expr item = array.element(i);
                if (!item.isSymbol())
                    error::throwCalcError(
                        error::CalcErrorType::Type,
                        "collect variable array must contain only symbols");
                variables.push_back(item.asSymbol());
            }
        }
        else
            error::throwCalcError(
                error::CalcErrorType::Type,
                "collect expects an expression and a symbol or symbol array");
        return symbolic::collectExpression(
            arguments[0], variables, registry_, mathematics_, angleSemantics_);
    }
    case BuiltinId::Solve: {
        if (arguments.size() < 2 || arguments.size() > 3)
            error::throwCalcError(
                error::CalcErrorType::Type,
                "solve expects equation(s), variable(s), and optional constraints");

        std::vector<expression::Symbol> variables;
        if (arguments[1].isSymbol())
            variables.push_back(arguments[1].asSymbol());
        else if (arguments[1].isArray() && arguments[1].asArray().rank() == 1) {
            const auto& array = arguments[1].asArray();
            for (std::size_t i = 0; i < array.size(); ++i) {
                const expression::Expr item = array.element(i);
                if (!item.isSymbol())
                    error::throwCalcError(
                        error::CalcErrorType::Type,
                        "solve variable array must contain only symbols");
                if (std::find(variables.begin(), variables.end(), item.asSymbol()) != variables.end())
                    error::throwCalcError(
                        error::CalcErrorType::Type,
                        "solve variable array contains a duplicate symbol");
                variables.push_back(item.asSymbol());
            }
        }
        else
            error::throwCalcError(
                error::CalcErrorType::Type,
                "solve expects a symbol or symbol array as its second argument");

        solver::SolveConstraints constraints;
        if (arguments.size() == 3)
            constraints = solver::parseSolveConstraints(
                arguments[2], variables, registry_, mathematics_, angleSemantics_);

        solver::SolutionSet solutions = [&]() {
            if (variables.size() == 1 && !arguments[0].isArray()) {
                const bool realDomain = constraints.domain
                    && mathematics::isSubdomainOf(
                        *constraints.domain, mathematics::NumericDomain::Real);
                if (realDomain) {
                    if (auto transcendental = solver::solveRealInjectiveFunctionRelation(
                            arguments[0], variables.front(), registry_, mathematics_,
                            angleSemantics_, constraints.assumptions))
                        return *transcendental;
                }
                return solver::solveUnivariatePolynomialRelation(
                    arguments[0], variables.front(), registry_, mathematics_, angleSemantics_);
            }

            std::vector<expression::Expr> equations;
            if (arguments[0].isArray() && arguments[0].asArray().rank() == 1)
                equations = arguments[0].asArray().materialize();
            else
                equations.push_back(arguments[0]);

            // 一変数のrelation配列は論理積として扱う。等式があれば先に解いて有限候補を作り、残りをexact constraintとして絞る。
            // 等式がなければ最初の不等式からReal領域branchを作り、残りの不等式を条件として交差させる。
            if (variables.size() == 1 && !equations.empty()) {
                auto first = equations.begin();
                const auto equality = std::find_if(
                    equations.begin(), equations.end(), [&](const expression::Expr& item) {
                        return item.isCall()
                            && item.asCall().head.sameIdentity(registry_.symbol(BuiltinId::Equal));
                    });
                if (equality != equations.end())
                    first = equality;

                solver::SolutionSet result = solver::solveUnivariatePolynomialRelation(
                    *first, variables.front(), registry_, mathematics_, angleSemantics_);
                for (auto iterator = equations.begin(); iterator != equations.end(); ++iterator) {
                    if (iterator == first)
                        continue;
                    const std::array<expression::Symbol, 1> oneVariable{variables.front()};
                    const solver::SolveConstraints relationConstraint =
                        solver::parseSolveConstraints(
                            *iterator, oneVariable, registry_, mathematics_, angleSemantics_);
                    result = solver::applySolveConstraints(
                        std::move(result), relationConstraint,
                        registry_, mathematics_, angleSemantics_);
                }
                return result;
            }

            return solver::solveLinearPolynomialSystem(
                equations, variables, registry_, mathematics_, angleSemantics_);
        }();

        if (constraints.domain
            && *constraints.domain == mathematics::NumericDomain::Complex
            && !solutions.variables().empty()
            && solutions.variables().front().domain == mathematics::NumericDomain::Real) {
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "Ordered inequalities are defined only over Real or a subdomain");
        }

        solutions = solver::applySolveConstraints(
            std::move(solutions), constraints, registry_, mathematics_, angleSemantics_);
        if (containsUnresolvedSolution(solutions))
            emitWarning("solve::unresolved",
                "solve could not determine a complete solution set; unresolved cases remain");
        return expression::Expr::solutionSet(std::move(solutions));
    }
    case BuiltinId::Set:
        return evaluateSet(arguments);
    case BuiltinId::SetDelayed:
        return evaluateSetDelayed(call, arguments);
    case BuiltinId::Less:
    case BuiltinId::LessEqual:
    case BuiltinId::Greater:
    case BuiltinId::GreaterEqual:
    case BuiltinId::Equal:
    case BuiltinId::NotEqual:
        return builtins::evaluateComparison(registry_.symbol(definition.id), arguments);
    case BuiltinId::LogicalAnd:
        return builtins::evaluateLogicalAnd(arguments, registry_);
    case BuiltinId::Element:
        return builtins::evaluateElement(arguments, registry_, mathematics_);
    case BuiltinId::If:
        error::throwCalcError(
            error::CalcErrorType::Internal,
            "If must be handled by the evaluation machine");
    case BuiltinId::History:
        return evaluateHistory(arguments);
    case BuiltinId::InputHistory:
        return evaluateIndexedHistory(arguments, true);
    case BuiltinId::OutputHistory:
        return evaluateIndexedHistory(arguments, false);
    case BuiltinId::Exit:
        if (context_ && context_->exitRequested)
            *context_->exitRequested = true;
        return expression::Expr{true};
    case BuiltinId::Clear:
        if (context_ && context_->clearRequested)
            *context_->clearRequested = true;
        return expression::Expr{true};
    case BuiltinId::Definitions:
        return evaluateDefinitions();
    case BuiltinId::Undefine:
        return evaluateUndefine(arguments);
    case BuiltinId::AngleMode: {
        const auto angleSymbol = [&](mathematics::AngleUnit unit) -> expression::Expr {
            std::string_view name;
            switch (unit) {
            case mathematics::AngleUnit::Degree: name = "Deg"; break;
            case mathematics::AngleUnit::Radian: name = "Rad"; break;
            case mathematics::AngleUnit::Gradian: name = "Grad"; break;
            }

            const auto* predefined = symbolRegistry_.find(name);
            if (!predefined)
                error::throwCalcError(
                    error::CalcErrorType::Internal,
                    "Angle-mode symbol is not registered");
            return expression::Expr{predefined->symbol};
        };

        if (arguments.empty())
            return angleSymbol(angleSemantics_.defaultUnit());

        if (arguments.size() != 1 || !arguments.front().isSymbol())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "angleMode expects no argument or one of Rad, Deg, Grad");

        const auto* predefined = symbolRegistry_.find(arguments.front().asSymbol());
        if (!predefined)
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "angleMode expects one of Rad, Deg, Grad");

        mathematics::AngleUnit unit;
        switch (predefined->id) {
        case symbols::PredefinedSymbolId::DegreeUnit:
            unit = mathematics::AngleUnit::Degree;
            break;
        case symbols::PredefinedSymbolId::RadianUnit:
            unit = mathematics::AngleUnit::Radian;
            break;
        case symbols::PredefinedSymbolId::GradianUnit:
            unit = mathematics::AngleUnit::Gradian;
            break;
        default:
            error::throwCalcError(
                error::CalcErrorType::Domain,
                "angleMode expects one of Rad, Deg, Grad");
        }

        if (!context_ || !context_->angleSemantics)
            error::throwCalcError(
                error::CalcErrorType::Internal,
                "angleMode requires a mutable kernel session");
        context_->angleSemantics->setDefaultUnit(unit);
        return angleSymbol(unit);
    }
    case BuiltinId::UnitApplied: {
        if (arguments.size() != 2 || !arguments[1].isString())
            error::throwCalcError(
                error::CalcErrorType::Type,
                "UnitApplied requires a value and unit name");

        // 角度単位だけは数学層で意味を持つため、綴りを正規化する。
        // 長さ等の単位はまだ演算しないが、将来のunit systemへ渡せるようUnitApplied式として保持し、評価エラーにはしない。
        std::string unit = arguments[1].asString();
        if (const auto angleUnit = mathematics::AngleSemantics::parseUnit(unit))
            unit = std::string{mathematics::AngleSemantics::canonicalName(*angleUnit)};
        return expression::Expr::call(
            registry_.symbol(BuiltinId::UnitApplied),
            {arguments[0], expression::Expr{std::move(unit)}});
    }
    }

    error::throwCalcError(
        error::CalcErrorType::Internal,
        "Builtin dispatch is incomplete");
}

} // namespace mmcal::evaluation
