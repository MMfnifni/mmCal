// Arrayを共通表現とするexact-first線形代数の公開builtin層
#include "linear_algebra.hpp"

#include "approximation/expression_interval.hpp"
#include "builtins/array_helpers.hpp"
#include "builtins/exact_operations.hpp"
#include "error/error_message.hpp"
#include "expression/array_utils.hpp"
#include "linear_algebra/approximate_matrix.hpp"
#include "linear_algebra/decomposition.hpp"
#include "linear_algebra/eigen.hpp"
#include "linear_algebra/exact_matrix.hpp"
#include "linear_algebra/matrix.hpp"
#include "linear_algebra/svd.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/real_number.hpp"
#include "numeric/number.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] bool allNumbers(const ArrayExpr& array) noexcept {
    return array.hasExactNumberStorage();
}

[[nodiscard]] Expr transposeArray(const ArrayExpr& array) {
    return Expr::array(array.transposed());
}

[[nodiscard]] Expr productTerm(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() * rhs.asNumber()};
    return exact::multiply({lhs, rhs}, registry, mathematics, angles);
}

[[nodiscard]] Expr sumTerms(
    std::vector<Expr> terms,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (terms.empty())
        return integer(0);

    bool numeric = true;
    Number sum{BigInt{0}};
    for (const Expr& term : terms) {
        if (!term.isNumber()) {
            numeric = false;
            break;
        }
        sum += term.asNumber();
    }
    if (numeric)
        return Expr{std::move(sum)};
    return exact::add(std::move(terms), registry, mathematics, angles);
}

[[nodiscard]] Expr dotCell(
    const ArrayExpr& lhs,
    const ArrayExpr& rhs,
    std::size_t row,
    std::size_t column,
    std::size_t inner,
    std::size_t lhsColumns,
    std::size_t rhsColumns,
    bool numericInputs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (numericInputs) {
        Number sum{BigInt{0}};
        for (std::size_t k = 0; k < inner; ++k) {
            const Number left = lhs.exactNumber(lhs.rank() == 1
                ? k
                : row * lhsColumns + k);
            const Number right = rhs.exactNumber(rhs.rank() == 1
                ? k
                : k * rhsColumns + column);
            sum += left * right;
        }
        return Expr{std::move(sum)};
    }

    std::vector<Expr> terms;
    terms.reserve(inner);
    for (std::size_t k = 0; k < inner; ++k) {
        const Expr left = lhs.element(lhs.rank() == 1
            ? k
            : row * lhsColumns + k);
        const Expr right = rhs.element(rhs.rank() == 1
            ? k
            : k * rhsColumns + column);
        terms.push_back(productTerm(left, right, registry, mathematics, angles));
    }
    return sumTerms(std::move(terms), registry, mathematics, angles);
}

[[nodiscard]] Expr hermitianNorm(
    const ArrayExpr& vector,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (vector.empty())
        return integer(0);

    if (allNumbers(vector)) {
        Number sum{BigInt{0}};
        for (std::size_t i = 0; i < vector.size(); ++i) {
            const Number element = vector.exactNumber(i);
            sum += element.conjugate() * element;
        }
        return exact::sqrt(Expr{std::move(sum)}, registry, mathematics, angles);
    }

    std::vector<Expr> terms;
    terms.reserve(vector.size());
    for (std::size_t i = 0; i < vector.size(); ++i) {
        const Expr element = vector.element(i);
        const mathematics::ValueFacts facts = mathematics::inferValueFacts(
            element, registry, mathematics);
        Expr conjugate = facts.isProvablyReal()
            ? element
            : exact::call(BuiltinId::Conj, {element}, registry, mathematics, angles);
        terms.push_back(productTerm(conjugate, element, registry, mathematics, angles));
    }
    return exact::sqrt(sumTerms(std::move(terms), registry, mathematics, angles),
        registry, mathematics, angles);
}

[[nodiscard]] bool provablyNonZero(
    const Expr& value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    if (value.isNumber())
        return !value.asNumber().isZero();
    const auto facts = mathematics::inferValueFacts(value, registry, mathematics);
    return facts.sign == mathematics::RealSign::Positive
        || facts.sign == mathematics::RealSign::Negative
        || facts.sign == mathematics::RealSign::NonZero;
}

[[nodiscard]] linear_algebra::ExactMatrixContext exactContext(
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return {registry, mathematics, angles};
}

} // namespace

Expr evaluateTranspose(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry&) {
    const ArrayExpr& array = detail::requireArray(arguments.front(), "transpose");
    if (array.rank() == 1)
        return arguments.front();
    if (array.rank() != 2)
        detail::arrayTypeError("transpose currently supports rank-1 or rank-2 arrays");

    linear_algebra::MatrixView matrix{array};
    if (matrix.size() == 0)
        return Expr::array({matrix.columns(), matrix.rows()}, {});
    return transposeArray(array);
}

Expr evaluateConjugateTranspose(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireArray(arguments.front(), "conjugateTranspose");
    if (array.rank() != 1 && array.rank() != 2)
        detail::arrayTypeError("conjugateTranspose supports rank-1 or rank-2 arrays");

    const auto conjugated = [&](const Expr& value) -> Expr {
        if (value.isNumber())
            return Expr{value.asNumber().conjugate()};
        if (value.isDecimalApproximation())
            return value;
        if (value.isComplexDecimalApproximation()) {
            const auto& complex = value.asComplexDecimalApproximation();
            const auto& imaginary = complex.imaginary();
            numeric::DecimalApproximation conjugateImaginary = [&] {
                if (imaginary.origin() == numeric::ApproximationOrigin::ExactValue)
                    return numeric::DecimalApproximation::fromReal(
                        numeric::RealNumber{-imaginary.displayedValue()},
                        imaginary.requestedFractionalDigits());
                const auto result = numeric::DecimalApproximation::fromCertifiedInterval(
                    -imaginary.certifiedUpper(), -imaginary.certifiedLower(),
                    imaginary.requestedFractionalDigits());
                if (!result)
                    throw std::logic_error(
                        "Conjugating a certified decimal approximation must preserve rounding");
                return *result;
            }();
            return Expr{numeric::ComplexDecimalApproximation::fromComponents(
                complex.real(), std::move(conjugateImaginary),
                complex.realExactlyZero(), complex.imaginaryExactlyZero())};
        }
        const mathematics::ValueFacts facts = mathematics::inferValueFacts(
            value, registry, mathematics);
        if (facts.isProvablyReal())
            return value;
        return exact::call(BuiltinId::Conj, {value}, registry, mathematics, angles);
    };

    if (array.rank() == 1) {
        std::vector<Expr> output;
        output.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i)
            output.push_back(conjugated(array.element(i)));
        return Expr::array(array.shape, std::move(output));
    }

    linear_algebra::MatrixView matrix{array};
    std::vector<Expr> output;
    output.reserve(matrix.size());
    for (std::size_t column = 0; column < matrix.columns(); ++column)
        for (std::size_t row = 0; row < matrix.rows(); ++row)
            output.push_back(conjugated(matrix(row, column)));
    return Expr::array({matrix.columns(), matrix.rows()}, std::move(output));
}

Expr evaluateMatrixAdd(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& first = detail::requireArray(arguments.front(), "madd");

    bool allExact = first.hasExactNumberStorage();
    for (std::size_t a = 1; a < arguments.size(); ++a) {
        const ArrayExpr& next = detail::requireArray(arguments[a], "madd");
        if (next.shape != first.shape)
            error::throwCalcError(error::CalcErrorType::Domain,
                "madd requires identical array shapes");
        allExact = allExact && next.hasExactNumberStorage();
    }
    if (allExact) {
        std::vector<Number> values(first.size(), Number{BigInt{0}});
        for (const Expr& argument : arguments) {
            const ArrayExpr& array = argument.asArray();
            for (std::size_t i = 0; i < values.size(); ++i)
                values[i] += array.exactNumber(i);
        }
        return Expr::numberArray(first.shape, std::move(values));
    }

    std::vector<Expr> output;
    output.reserve(first.size());

    for (std::size_t i = 0; i < first.size(); ++i) {
        std::vector<Expr> terms;
        terms.reserve(arguments.size());
        terms.push_back(first.element(i));
        for (std::size_t a = 1; a < arguments.size(); ++a) {
            const ArrayExpr& next = detail::requireArray(arguments[a], "madd");
            terms.push_back(next.element(i));
        }
        output.push_back(sumTerms(std::move(terms), registry, mathematics, angles));
    }
    return Expr::array(first.shape, std::move(output));
}

Expr evaluateDot(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateDot(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& lhs = detail::requireArray(arguments[0], "dot");
    const ArrayExpr& rhs = detail::requireArray(arguments[1], "dot");
    if (lhs.rank() < 1 || lhs.rank() > 2 || rhs.rank() < 1 || rhs.rank() > 2)
        detail::arrayTypeError("dot currently supports vectors and matrices");

    const std::size_t lhsRows = lhs.rank() == 1 ? 1 : lhs.shape[0];
    const std::size_t lhsColumns = lhs.rank() == 1 ? lhs.shape[0] : lhs.shape[1];
    const std::size_t rhsRows = rhs.shape[0];
    const std::size_t rhsColumns = rhs.rank() == 1 ? 1 : rhs.shape[1];
    if (lhsColumns != rhsRows)
        error::throwCalcError(error::CalcErrorType::Domain,
            "dot inner dimensions do not agree");

    const bool numericInputs = allNumbers(lhs) && allNumbers(rhs);
    const std::size_t outputShape[] = {lhsRows, rhsColumns};
    std::size_t outputSize = 0;
    try {
        outputSize = expression::arrayElementCount(outputShape);
    }
    catch (const std::length_error&) {
        error::throwCalcError(error::CalcErrorType::Overflow,
            "dot result dimensions overflow the addressable element count");
    }

    if (numericInputs) {
        std::vector<Number> output;
        output.reserve(outputSize);
        for (std::size_t row = 0; row < lhsRows; ++row) {
            for (std::size_t column = 0; column < rhsColumns; ++column) {
                Number sum{BigInt{0}};
                for (std::size_t k = 0; k < lhsColumns; ++k) {
                    const Number left = lhs.exactNumber(lhs.rank() == 1
                        ? k
                        : row * lhsColumns + k);
                    const Number right = rhs.exactNumber(rhs.rank() == 1
                        ? k
                        : k * rhsColumns + column);
                    sum += left * right;
                }
                output.push_back(std::move(sum));
            }
        }
        if (lhs.rank() == 1 && rhs.rank() == 1)
            return Expr{std::move(output.front())};
        if (lhs.rank() == 1 || rhs.rank() == 1)
            return Expr::numberArray(
                {lhs.rank() == 1 ? rhsColumns : lhsRows}, std::move(output));
        return Expr::numberArray({lhsRows, rhsColumns}, std::move(output));
    }

    std::vector<Expr> output;
    output.reserve(outputSize);
    for (std::size_t row = 0; row < lhsRows; ++row)
        for (std::size_t column = 0; column < rhsColumns; ++column)
            output.push_back(dotCell(lhs, rhs, row, column, lhsColumns,
                lhsColumns, rhsColumns, false, registry, mathematics, angles));

    if (lhs.rank() == 1 && rhs.rank() == 1)
        return output.front();
    if (lhs.rank() == 1 || rhs.rank() == 1)
        return Expr::array({lhs.rank() == 1 ? rhsColumns : lhsRows}, std::move(output));
    return Expr::array({lhsRows, rhsColumns}, std::move(output));
}

Expr evaluateDeterminant(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateDeterminant(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "det");
    if (array.shape[0] != array.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain, "det requires a square matrix");
    const auto result = linear_algebra::determinant(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (!result)
        return Expr::call(registry.symbol(BuiltinId::Determinant), {arguments.front()});
    return *result;
}

Expr evaluateMatrixInverse(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateInverse(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "inverse");
    if (array.shape[0] != array.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain,
            "inverse requires a square matrix");
    const auto result = linear_algebra::inverse(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (!result)
        return Expr::call(registry.symbol(BuiltinId::Inverse), {arguments.front()});
    return *result;
}

Expr evaluateRref(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateRref(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "rref");
    const auto reduced = linear_algebra::rref(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (!reduced)
        return Expr::call(registry.symbol(BuiltinId::Rref), {arguments.front()});
    return reduced->toExpr();
}

Expr evaluateMatrixRank(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateMatrixRank(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "matrixRank");
    const auto rank = linear_algebra::matrixRank(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (!rank)
        return Expr::call(registry.symbol(BuiltinId::Rank), {arguments.front()});
    return detail::sizeExpr(*rank);
}

Expr evaluateSolveLinear(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateSolveLinear(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& matrix = detail::requireMatrix(arguments[0], "solveLinear");
    const ArrayExpr& rhs = detail::requireVector(arguments[1], "solveLinear");
    if (rhs.shape[0] != matrix.shape[0])
        error::throwCalcError(error::CalcErrorType::Domain,
            "solveLinear right-hand side size must match the matrix row count");

    const auto result = linear_algebra::solveLinear(
        linear_algebra::MatrixView{matrix}, rhs,
        exactContext(registry, mathematics, angles));
    if (!result)
        return Expr::call(registry.symbol(BuiltinId::SolveLinear),
            {arguments[0], arguments[1]});
    return *result;
}

Expr evaluateNullSpace(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "nullSpace");
    const auto result = linear_algebra::nullSpace(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (result)
        return *result;
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto approximate = evaluateApproximateNullSpace(
            arguments, registry, mathematics, angles, *context))
            return *approximate;
    return Expr::call(registry.symbol(BuiltinId::NullSpace), {arguments.front()});
}

Expr evaluateLuDecomposition(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "luDecomposition");
    if (array.shape[0] != array.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain,
            "luDecomposition currently requires a square matrix");
    const auto result = linear_algebra::luDecomposition(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (result)
        return *result;
    return Expr::call(registry.symbol(BuiltinId::LuDecomposition), {arguments.front()});
}

Expr evaluateQrDecomposition(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "qrDecomposition");
    const auto result = linear_algebra::qrDecomposition(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (result)
        return *result;
    return Expr::call(registry.symbol(BuiltinId::QrDecomposition), {arguments.front()});
}

Expr evaluateSingularValueDecomposition(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "svd");
    const auto result = linear_algebra::singularValueDecomposition(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (result)
        return *result;
    return Expr::call(registry.symbol(BuiltinId::SingularValueDecomposition), {arguments.front()});
}

Expr evaluateEigenvalues(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "eigenvalues");
    if (array.shape[0] != array.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain,
            "eigenvalues requires a square matrix");
    const auto result = linear_algebra::eigenvalues(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (result)
        return *result;
    return Expr::call(registry.symbol(BuiltinId::Eigenvalues), {arguments.front()});
}

Expr evaluateEigenvectors(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "eigenvectors");
    if (array.shape[0] != array.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain,
            "eigenvectors requires a square matrix");
    const auto result = linear_algebra::eigenvectors(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (result)
        return *result;
    return Expr::call(registry.symbol(BuiltinId::Eigenvectors), {arguments.front()});
}

Expr evaluateEigensystem(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "eigensystem");
    if (array.shape[0] != array.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain,
            "eigensystem requires a square matrix");
    const auto result = linear_algebra::eigensystem(
        linear_algebra::MatrixView{array}, exactContext(registry, mathematics, angles));
    if (result)
        return *result;
    return Expr::call(registry.symbol(BuiltinId::Eigensystem), {arguments.front()});
}

Expr evaluateNorm(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateNorm(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& vector = detail::requireVector(arguments.front(), "norm");
    return hermitianNorm(vector, registry, mathematics, angles);
}

Expr evaluateNormalize(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateNormalize(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& vector = detail::requireVector(arguments.front(), "normalize");
    const Expr norm = hermitianNorm(vector, registry, mathematics, angles);
    if (norm.isNumber() && norm.asNumber().isZero())
        error::throwCalcError(error::CalcErrorType::Domain,
            "normalize requires a nonzero vector");
    if (!provablyNonZero(norm, registry, mathematics))
        return Expr::call(registry.symbol(BuiltinId::VectorNormalize), {arguments.front()});

    const Expr reciprocal = exact::divide(integer(1), norm, registry, mathematics, angles);
    std::vector<Expr> output;
    output.reserve(vector.size());
    for (std::size_t i = 0; i < vector.size(); ++i) {
        const Expr element = vector.element(i);
        output.push_back(productTerm(element, reciprocal, registry, mathematics, angles));
    }
    return Expr::array(vector.shape, std::move(output));
}

Expr evaluateTrace(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto context = approximation::inferredApproximationContext(arguments))
        if (const auto result = evaluateApproximateTrace(
            arguments, registry, mathematics, angles, *context))
            return *result;

    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "trace");
    if (array.shape[0] != array.shape[1])
        error::throwCalcError(error::CalcErrorType::Domain,
            "trace requires a square matrix");
    linear_algebra::MatrixView matrix{array};
    std::vector<Expr> diagonal;
    diagonal.reserve(matrix.rows());
    for (std::size_t i = 0; i < matrix.rows(); ++i)
        diagonal.push_back(matrix(i, i));
    return sumTerms(std::move(diagonal), registry, mathematics, angles);
}

std::optional<Expr> evaluateApproximateDot(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateDot(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateDeterminant(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateDeterminant(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateInverse(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateInverse(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateRref(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateRref(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateMatrixRank(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateMatrixRank(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateSolveLinear(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateSolveLinear(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateNullSpace(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateNullSpace(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateLuDecomposition(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray())
        return std::nullopt;
    return linear_algebra::approximateLuDecomposition(
        arguments.front().asArray(), registry, mathematics, angles, context);
}

std::optional<Expr> evaluateApproximateQrDecomposition(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray())
        return std::nullopt;
    return linear_algebra::approximateQrDecomposition(
        arguments.front().asArray(), registry, mathematics, angles, context);
}

std::optional<Expr> evaluateApproximateSingularValueDecomposition(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray())
        return std::nullopt;
    return linear_algebra::approximateSingularValueDecomposition(
        arguments.front().asArray(), registry, mathematics, angles, context);
}

std::optional<Expr> evaluateApproximateEigenvalues(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray())
        return std::nullopt;
    return linear_algebra::approximateEigenvalues(
        arguments.front().asArray(), registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateEigenvectors(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray())
        return std::nullopt;
    return linear_algebra::approximateEigenvectors(
        arguments.front().asArray(), registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateEigensystem(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray())
        return std::nullopt;
    return linear_algebra::approximateEigensystem(
        arguments.front().asArray(), registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateNorm(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateNorm(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateNormalize(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateNormalize(
        arguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateTrace(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    return linear_algebra::approximateTrace(
        arguments, registry, mathematics, angles, std::move(context));
}

} // namespace mmcal::builtins
