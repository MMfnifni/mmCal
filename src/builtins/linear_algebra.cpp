// Arrayを共通表現とするexact-first線形代数の公開builtin層
#include "linear_algebra.hpp"

#include "approximation/expression_interval.hpp"
#include "builtins/array_helpers.hpp"
#include "builtins/exact_operations.hpp"
#include "error/error_message.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "expression/array_utils.hpp"
#include "linear_algebra/approximate_matrix.hpp"
#include "linear_algebra/decomposition.hpp"
#include "linear_algebra/eigen.hpp"
#include "linear_algebra/exact_matrix.hpp"
#include "linear_algebra/matrix.hpp"
#include "linear_algebra/svd.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/real_number.hpp"
#include "numeric/number.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
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

[[nodiscard]] bool containsFiniteApproximation(const ArrayExpr& array) {
    for (std::size_t i = 0; i < array.size(); ++i) {
        const Expr value = array.element(i);
        if (value.isDecimalApproximation() || value.isComplexDecimalApproximation())
            return true;
    }
    return false;
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

[[nodiscard]] Expr conjugateValue(
    const Expr& value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (value.isNumber())
        return Expr{value.asNumber().conjugate()};
    if (value.isDecimalApproximation())
        return value;
    if (value.isComplexDecimalApproximation()) {
        const auto& complex = value.asComplexDecimalApproximation();
        return Expr{numeric::ComplexDecimalApproximation::fromComponents(
            complex.real(), complex.imaginary().negated(),
            complex.realExactlyZero(), complex.imaginaryExactlyZero())};
    }
    const mathematics::ValueFacts facts = mathematics::inferValueFacts(
        value, registry, mathematics);
    if (facts.isProvablyReal())
        return value;
    return exact::call(BuiltinId::Conj, {value}, registry, mathematics, angles);
}

[[nodiscard]] Expr hermitianInner(
    const ArrayExpr& lhs,
    const ArrayExpr& rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs.size() != rhs.size())
        error::throwCalcError(error::CalcErrorType::Domain,
            "inner requires vectors with the same length");
    if (lhs.empty())
        return integer(0);

    if (allNumbers(lhs) && allNumbers(rhs)) {
        Number sum{BigInt{0}};
        for (std::size_t i = 0; i < lhs.size(); ++i)
            sum += lhs.exactNumber(i).conjugate() * rhs.exactNumber(i);
        return Expr{std::move(sum)};
    }

    std::vector<Expr> terms;
    terms.reserve(lhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i)
        terms.push_back(productTerm(
            conjugateValue(lhs.element(i), registry, mathematics, angles),
            rhs.element(i), registry, mathematics, angles));
    return sumTerms(std::move(terms), registry, mathematics, angles);
}

[[nodiscard]] Expr hermitianNorm(
    const ArrayExpr& vector,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return exact::sqrt(
        hermitianInner(vector, vector, registry, mathematics, angles),
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

[[nodiscard]] bool provablyZero(
    const Expr& value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    if (value.isNumber())
        return value.asNumber().isZero();
    return mathematics::inferValueFacts(value, registry, mathematics).sign
        == mathematics::RealSign::Zero;
}

[[nodiscard]] Expr matrixRow(const ArrayExpr& matrix, std::size_t row) {
    std::vector<Expr> elements;
    elements.reserve(matrix.shape[1]);
    const std::size_t offset = row * matrix.shape[1];
    for (std::size_t column = 0; column < matrix.shape[1]; ++column)
        elements.push_back(matrix.element(offset + column));
    return Expr::array({matrix.shape[1]}, std::move(elements));
}

[[nodiscard]] Expr subtractVectors(
    const ArrayExpr& lhs,
    const ArrayExpr& rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> output;
    output.reserve(lhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i)
        output.push_back(exact::subtract(
            lhs.element(i), rhs.element(i), registry, mathematics, angles));
    return Expr::array(lhs.shape, std::move(output));
}

[[nodiscard]] Expr scaleVector(
    const ArrayExpr& vector,
    const Expr& factor,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<Expr> output;
    output.reserve(vector.size());
    for (std::size_t i = 0; i < vector.size(); ++i)
        output.push_back(productTerm(
            vector.element(i), factor, registry, mathematics, angles));
    return Expr::array(vector.shape, std::move(output));
}

[[nodiscard]] Expr gramSchmidtCore(
    const ArrayExpr& vectors,
    BuiltinId unresolvedId,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    bool requireFullRank,
    bool returnPredicate) {
    const std::size_t rowCount = vectors.shape[0];
    const std::size_t dimension = vectors.shape[1];
    if (requireFullRank && rowCount > dimension)
        return Expr{false};

    std::vector<Expr> basis;
    std::vector<Expr> squaredNorms;
    basis.reserve(std::min(rowCount, dimension));
    squaredNorms.reserve(std::min(rowCount, dimension));
    for (std::size_t row = 0; row < rowCount; ++row) {
        Expr residual = matrixRow(vectors, row);
        for (std::size_t i = 0; i < basis.size(); ++i) {
            const Expr numerator = hermitianInner(
                basis[i].asArray(), residual.asArray(), registry, mathematics, angles);
            const Expr coefficient = exact::divide(
                numerator, squaredNorms[i], registry, mathematics, angles);
            const Expr component = scaleVector(
                basis[i].asArray(), coefficient, registry, mathematics, angles);
            residual = subtractVectors(
                residual.asArray(), component.asArray(), registry, mathematics, angles);
        }

        const Expr squaredNorm = hermitianInner(
            residual.asArray(), residual.asArray(), registry, mathematics, angles);
        if (provablyZero(squaredNorm, registry, mathematics)) {
            if (requireFullRank)
                return Expr{false};
            continue;
        }
        if (!provablyNonZero(squaredNorm, registry, mathematics))
            return Expr::call(registry.symbol(unresolvedId), {Expr::array(vectors)});

        basis.push_back(std::move(residual));
        squaredNorms.push_back(squaredNorm);
    }

    if (returnPredicate)
        return Expr{basis.size() == rowCount};

    std::vector<Expr> output;
    output.reserve(basis.size() * dimension);
    for (std::size_t row = 0; row < basis.size(); ++row) {
        const Expr norm = exact::sqrt(
            squaredNorms[row], registry, mathematics, angles);
        const Expr reciprocal = exact::divide(
            integer(1), norm, registry, mathematics, angles);
        const Expr unit = scaleVector(
            basis[row].asArray(), reciprocal, registry, mathematics, angles);
        for (std::size_t i = 0; i < dimension; ++i)
            output.push_back(unit.asArray().element(i));
    }
    return Expr::array({basis.size(), dimension}, std::move(output));
}

[[nodiscard]] linear_algebra::ExactMatrixContext exactContext(
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return {registry, mathematics, angles};
}

[[nodiscard]] Expr infinity() {
    return Expr{expression::Symbol{"Infinity"}};
}

[[nodiscard]] Expr dotPair(
    const Expr& lhs,
    const Expr& rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::array<Expr, 2> arguments{lhs, rhs};
    return evaluateDot(arguments, registry, mathematics, angles);
}

[[nodiscard]] Expr conjugateTransposeValue(
    const Expr& value,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::array<Expr, 1> arguments{value};
    return evaluateConjugateTranspose(arguments, registry, mathematics, angles);
}

[[nodiscard]] std::optional<std::vector<std::size_t>> exactPivotColumns(
    const linear_algebra::MatrixBuffer& reduced) {
    std::vector<std::size_t> pivots;
    pivots.reserve(std::min(reduced.rows(), reduced.columns()));
    for (std::size_t row = 0; row < reduced.rows(); ++row) {
        for (std::size_t column = 0; column < reduced.columns(); ++column) {
            const Expr& value = reduced(row, column);
            if (!value.isNumber())
                return std::nullopt;
            if (!value.asNumber().isZero()) {
                pivots.push_back(column);
                break;
            }
        }
    }
    return pivots;
}

// exact Number行列ではrank factorization A=F Gを作り，
// A^+=G^H(GG^H)^-1(F^H F)^-1F^H をそのままexact算術で評価する。
[[nodiscard]] std::optional<Expr> exactPseudoInverse(
    const ArrayExpr& array,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const linear_algebra::MatrixView matrix{array};
    const std::size_t rows = matrix.rows();
    const std::size_t columns = matrix.columns();
    if (rows == 0 || columns == 0)
        return Expr::numberArray({columns, rows}, {});
    if (!linear_algebra::allExactNumbers(matrix))
        return std::nullopt;

    const auto context = exactContext(registry, mathematics, angles);
    const auto rank = linear_algebra::matrixRank(matrix, context);
    if (!rank)
        return std::nullopt;
    if (*rank == 0)
        return Expr::array({columns, rows},
            std::vector<Expr>(columns * rows, integer(0)));

    const auto reduced = linear_algebra::rref(matrix, context);
    if (!reduced)
        return std::nullopt;
    const auto pivots = exactPivotColumns(*reduced);
    if (!pivots || pivots->size() != *rank)
        return std::nullopt;

    std::vector<Expr> fValues;
    fValues.reserve(rows * *rank);
    for (std::size_t row = 0; row < rows; ++row)
        for (const std::size_t column : *pivots)
            fValues.push_back(matrix(row, column));
    Expr f = Expr::array({rows, *rank}, std::move(fValues));
    const linear_algebra::MatrixView fView{f.asArray()};

    std::vector<Expr> gValues(*rank * columns, integer(0));
    for (std::size_t column = 0; column < columns; ++column) {
        std::vector<Expr> rhsValues;
        rhsValues.reserve(rows);
        for (std::size_t row = 0; row < rows; ++row)
            rhsValues.push_back(matrix(row, column));
        const Expr rhs = Expr::array({rows}, std::move(rhsValues));
        const auto coefficients = linear_algebra::solveLinear(
            fView, rhs.asArray(), context);
        if (!coefficients || !coefficients->isArray()
            || !coefficients->asArray().isVector())
            return std::nullopt;
        for (std::size_t row = 0; row < *rank; ++row)
            gValues[row * columns + column] = coefficients->asArray().element(row);
    }
    Expr g = Expr::array({*rank, columns}, std::move(gValues));

    const Expr fh = conjugateTransposeValue(f, registry, mathematics, angles);
    const Expr gh = conjugateTransposeValue(g, registry, mathematics, angles);
    const Expr fGram = dotPair(fh, f, registry, mathematics, angles);
    const Expr gGram = dotPair(g, gh, registry, mathematics, angles);

    const std::array<Expr, 1> fGramArguments{fGram};
    const std::array<Expr, 1> gGramArguments{gGram};
    const Expr fGramInverse = evaluateMatrixInverse(
        fGramArguments, registry, mathematics, angles);
    const Expr gGramInverse = evaluateMatrixInverse(
        gGramArguments, registry, mathematics, angles);

    const Expr left = dotPair(gh, gGramInverse, registry, mathematics, angles);
    const Expr middle = dotPair(left, fGramInverse, registry, mathematics, angles);
    return dotPair(middle, fh, registry, mathematics, angles);
}

[[nodiscard]] std::optional<Expr> exactConditionNumber(
    const ArrayExpr& array,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const linear_algebra::MatrixView matrix{array};
    const std::size_t k = std::min(matrix.rows(), matrix.columns());
    if (k == 0 || !linear_algebra::allExactNumbers(matrix))
        return std::nullopt;

    const auto rank = linear_algebra::matrixRank(
        matrix, exactContext(registry, mathematics, angles));
    if (!rank)
        return std::nullopt;
    if (*rank < k)
        return infinity();
    if (k == 1)
        return integer(1);

    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            const Expr value = matrix(row, column);
            if (!value.isNumber() || !value.asNumber().isReal())
                return std::nullopt;
            if (row != column && !value.asNumber().isZero())
                return std::nullopt;
        }

    numeric::RealNumber minimum = matrix(0, 0).asNumber().asReal().abs();
    numeric::RealNumber maximum = minimum;
    for (std::size_t i = 1; i < k; ++i) {
        const numeric::RealNumber magnitude = matrix(i, i).asNumber().asReal().abs();
        if (magnitude < minimum)
            minimum = magnitude;
        if (magnitude > maximum)
            maximum = magnitude;
    }
    if (minimum.isZero())
        return infinity();
    return Expr{Number{maximum} / Number{minimum}};
}

[[nodiscard]] std::optional<std::array<Expr, 3>> svdFactors(const Expr& value) {
    if (value.isList()) {
        const auto& list = value.asList();
        if (list.size() != 3)
            return std::nullopt;
        return std::array<Expr, 3>{
            list.elements[0], list.elements[1], list.elements[2]};
    }
    if (!value.isArray())
        return std::nullopt;

    const ArrayExpr& array = value.asArray();
    if (array.rank() < 2 || array.shape.front() != 3)
        return std::nullopt;
    std::vector<std::size_t> factorShape(array.shape.begin() + 1, array.shape.end());
    const std::size_t factorSize = expression::arrayElementCount(factorShape);
    return std::array<Expr, 3>{
        Expr::array(array.sliced(factorShape, 0, factorSize)),
        Expr::array(array.sliced(factorShape, factorSize, factorSize)),
        Expr::array(array.sliced(std::move(factorShape), factorSize * 2, factorSize))};
}

[[nodiscard]] std::optional<std::size_t> approximateRank(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    const auto result = linear_algebra::approximateMatrixRank(
        arguments, registry, mathematics, angles, std::move(context));
    if (!result || !result->isNumber() || !result->asNumber().isReal()
        || !result->asNumber().asReal().isInteger())
        return std::nullopt;
    const auto converted = numeric::tryToUint64(
        result->asNumber().asReal().asInteger());
    if (!converted || *converted > std::numeric_limits<std::size_t>::max())
        return std::nullopt;
    return static_cast<std::size_t>(*converted);
}

[[nodiscard]] std::optional<Expr> approximateSvdPseudoInverse(
    const ArrayExpr& source,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    const std::size_t k = std::min(source.shape[0], source.shape[1]);
    if (k == 0)
        return Expr::numberArray({source.shape[1], source.shape[0]}, {});

    const std::array<Expr, 1> rankArguments{Expr::array(source)};
    const auto rank = approximateRank(
        rankArguments, registry, mathematics, angles, context);
    if (!rank)
        return std::nullopt;
    // 有限precisionでrank deficiencyを閾値推測しない。exact入力は後段の
    // rank-factorization経路でminimum-norm解を正確に作れる。
    if (*rank != k)
        return std::nullopt;

    const auto decomposition = linear_algebra::approximateSingularValueDecomposition(
        source, registry, mathematics, angles, context);
    if (!decomposition)
        return std::nullopt;
    const auto factors = svdFactors(*decomposition);
    if (!factors || !(*factors)[0].isArray() || !(*factors)[1].isArray()
        || !(*factors)[2].isArray())
        return std::nullopt;

    const ArrayExpr& sigma = (*factors)[1].asArray();
    if (!sigma.isMatrix() || sigma.shape[0] != k || sigma.shape[1] != k)
        return std::nullopt;
    std::vector<Expr> inverseSigmaValues(k * k, integer(0));
    for (std::size_t i = 0; i < k; ++i) {
        const Expr singular = sigma.element(i * k + i);
        const auto reciprocal = approximation::divideApproximateScalars(
            integer(1), singular);
        if (!reciprocal)
            return std::nullopt;
        inverseSigmaValues[i * k + i] = *reciprocal;
    }
    const Expr inverseSigma = Expr::array({k, k}, std::move(inverseSigmaValues));
    const Expr uh = conjugateTransposeValue(
        (*factors)[0], registry, mathematics, angles);

    const std::array<Expr, 2> firstArguments{(*factors)[2], inverseSigma};
    const auto first = evaluateApproximateDot(
        firstArguments, registry, mathematics, angles, context);
    if (!first)
        return std::nullopt;
    const std::array<Expr, 2> secondArguments{*first, uh};
    return evaluateApproximateDot(
        secondArguments, registry, mathematics, angles, std::move(context));
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

    if (array.rank() == 1) {
        std::vector<Expr> output;
        output.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i)
            output.push_back(conjugateValue(array.element(i), registry, mathematics, angles));
        return Expr::array(array.shape, std::move(output));
    }

    linear_algebra::MatrixView matrix{array};
    std::vector<Expr> output;
    output.reserve(matrix.size());
    for (std::size_t column = 0; column < matrix.columns(); ++column)
        for (std::size_t row = 0; row < matrix.rows(); ++row)
            output.push_back(conjugateValue(matrix(row, column), registry, mathematics, angles));
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
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::DenseArrayElement, outputSize);

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

Expr evaluateInner(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& lhs = detail::requireVector(arguments[0], "inner");
    const ArrayExpr& rhs = detail::requireVector(arguments[1], "inner");
    return hermitianInner(lhs, rhs, registry, mathematics, angles);
}

Expr evaluateOuter(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& lhs = detail::requireVector(arguments[0], "outer");
    const ArrayExpr& rhs = detail::requireVector(arguments[1], "outer");
    const std::size_t shape[] = {lhs.size(), rhs.size()};
    const std::size_t count = expression::arrayElementCount(shape);
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::DenseArrayElement, count);

    if (allNumbers(lhs) && allNumbers(rhs)) {
        std::vector<Number> output;
        output.reserve(count);
        for (std::size_t i = 0; i < lhs.size(); ++i)
            for (std::size_t j = 0; j < rhs.size(); ++j)
                output.push_back(lhs.exactNumber(i) * rhs.exactNumber(j));
        return Expr::numberArray({lhs.size(), rhs.size()}, std::move(output));
    }

    std::vector<Expr> output;
    output.reserve(count);
    for (std::size_t i = 0; i < lhs.size(); ++i)
        for (std::size_t j = 0; j < rhs.size(); ++j)
            output.push_back(productTerm(
                lhs.element(i), rhs.element(j), registry, mathematics, angles));
    return Expr::array({lhs.size(), rhs.size()}, std::move(output));
}

Expr evaluateDistance(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& lhs = detail::requireVector(arguments[0], "distance");
    const ArrayExpr& rhs = detail::requireVector(arguments[1], "distance");
    if (lhs.size() != rhs.size())
        error::throwCalcError(error::CalcErrorType::Domain,
            "distance requires vectors with the same length");

    std::vector<Expr> difference;
    difference.reserve(lhs.size());
    for (std::size_t i = 0; i < lhs.size(); ++i)
        difference.push_back(exact::subtract(
            lhs.element(i), rhs.element(i), registry, mathematics, angles));
    Expr vector = Expr::array({lhs.size()}, std::move(difference));
    return hermitianNorm(vector.asArray(), registry, mathematics, angles);
}

Expr evaluateProjectionForOperation(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    std::string_view operationName,
    BuiltinId unresolvedId) {
    const ArrayExpr& vector = detail::requireVector(arguments[0], operationName);
    const ArrayExpr& onto = detail::requireVector(arguments[1], operationName);
    if (vector.size() != onto.size())
        error::throwCalcError(error::CalcErrorType::Domain,
            std::string{operationName} + " requires vectors with the same length");

    const Expr denominator = hermitianInner(onto, onto, registry, mathematics, angles);
    if (denominator.isNumber() && denominator.asNumber().isZero())
        error::throwCalcError(error::CalcErrorType::Domain,
            std::string{operationName} + " requires a nonzero direction vector");
    if (!provablyNonZero(denominator, registry, mathematics))
        return Expr::call(registry.symbol(unresolvedId),
            {arguments[0], arguments[1]});

    const Expr factor = exact::divide(
        hermitianInner(onto, vector, registry, mathematics, angles),
        denominator, registry, mathematics, angles);
    std::vector<Expr> output;
    output.reserve(onto.size());
    for (std::size_t i = 0; i < onto.size(); ++i)
        output.push_back(productTerm(
            onto.element(i), factor, registry, mathematics, angles));
    return Expr::array(onto.shape, std::move(output));
}

Expr evaluateProjection(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return evaluateProjectionForOperation(arguments, registry, mathematics, angles,
        "projection", BuiltinId::VectorProject);
}

Expr evaluateRejection(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& vector = detail::requireVector(arguments[0], "rejection");
    const ArrayExpr& onto = detail::requireVector(arguments[1], "rejection");
    if (vector.size() != onto.size())
        error::throwCalcError(error::CalcErrorType::Domain,
            "rejection requires vectors with the same length");

    const Expr projected = evaluateProjectionForOperation(arguments, registry, mathematics, angles,
        "rejection", BuiltinId::VectorRejection);
    if (!projected.isArray())
        return Expr::call(registry.symbol(BuiltinId::VectorRejection),
            {arguments[0], arguments[1]});
    return subtractVectors(vector, projected.asArray(), registry, mathematics, angles);
}

Expr evaluateOrthogonalQ(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& vectors = detail::requireMatrix(arguments[0], "orthogonalQ");
    for (std::size_t i = 0; i < vectors.shape[0]; ++i) {
        const Expr lhs = matrixRow(vectors, i);
        for (std::size_t j = i + 1; j < vectors.shape[0]; ++j) {
            const Expr rhs = matrixRow(vectors, j);
            const Expr product = hermitianInner(
                lhs.asArray(), rhs.asArray(), registry, mathematics, angles);
            if (provablyZero(product, registry, mathematics))
                continue;
            if (provablyNonZero(product, registry, mathematics))
                return Expr{false};
            return Expr::call(registry.symbol(BuiltinId::OrthogonalQ), {arguments[0]});
        }
    }
    return Expr{true};
}

Expr evaluateOrthonormalQ(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& vectors = detail::requireMatrix(arguments[0], "orthonormalQ");
    for (std::size_t i = 0; i < vectors.shape[0]; ++i) {
        const Expr lhs = matrixRow(vectors, i);
        const Expr normSquared = hermitianInner(
            lhs.asArray(), lhs.asArray(), registry, mathematics, angles);
        const Expr unitDifference = exact::subtract(
            normSquared, integer(1), registry, mathematics, angles);
        if (!provablyZero(unitDifference, registry, mathematics)) {
            if (provablyNonZero(unitDifference, registry, mathematics))
                return Expr{false};
            return Expr::call(registry.symbol(BuiltinId::OrthonormalQ), {arguments[0]});
        }
        for (std::size_t j = i + 1; j < vectors.shape[0]; ++j) {
            const Expr rhs = matrixRow(vectors, j);
            const Expr product = hermitianInner(
                lhs.asArray(), rhs.asArray(), registry, mathematics, angles);
            if (provablyZero(product, registry, mathematics))
                continue;
            if (provablyNonZero(product, registry, mathematics))
                return Expr{false};
            return Expr::call(registry.symbol(BuiltinId::OrthonormalQ), {arguments[0]});
        }
    }
    return Expr{true};
}

Expr evaluateLinearIndependentQ(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& vectors = detail::requireMatrix(arguments[0], "linearIndependentQ");
    return gramSchmidtCore(vectors, BuiltinId::LinearIndependentQ,
        registry, mathematics, angles, true, true);
}

Expr evaluateGramSchmidt(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& vectors = detail::requireMatrix(arguments[0], "gramSchmidt");
    return gramSchmidtCore(vectors, BuiltinId::GramSchmidt,
        registry, mathematics, angles, false, false);
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
    if (const auto context = approximation::inferredApproximationContext(arguments)) {
        if (const auto result = evaluateApproximateInverse(
            arguments, registry, mathematics, angles, *context))
            return *result;
        // 有限precision入力で特異性を証明できない場合，hidden certified pointを
        // 用いた記号逆行列へ逃げない。入力情報のまま未評価に保つ。
        return Expr::call(registry.symbol(BuiltinId::Inverse), {arguments.front()});
    }

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

Expr evaluateConditionNumber(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "conditionNumber");
    if (array.shape[0] == 0 || array.shape[1] == 0)
        error::throwCalcError(error::CalcErrorType::Domain,
            "conditionNumber requires a non-empty matrix");
    if (const auto result = exactConditionNumber(array, registry, mathematics, angles))
        return *result;
    return Expr::call(registry.symbol(BuiltinId::ConditionNumber), {arguments.front()});
}

Expr evaluatePseudoInverse(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& array = detail::requireMatrix(arguments.front(), "pseudoInverse");
    if (const auto result = exactPseudoInverse(array, registry, mathematics, angles))
        return *result;
    return Expr::call(registry.symbol(BuiltinId::PseudoInverse), {arguments.front()});
}

Expr evaluateLeastSquares(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& matrix = detail::requireMatrix(arguments[0], "leastSquares");
    const ArrayExpr& rhs = detail::requireVector(arguments[1], "leastSquares");
    if (rhs.shape[0] != matrix.shape[0])
        error::throwCalcError(error::CalcErrorType::Domain,
            "leastSquares right-hand side size must match the matrix row count");

    const auto inverse = exactPseudoInverse(matrix, registry, mathematics, angles);
    if (!inverse)
        return Expr::call(registry.symbol(BuiltinId::LeastSquares),
            {arguments[0], arguments[1]});
    return dotPair(*inverse, arguments[1], registry, mathematics, angles);
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
    if (arguments.size() != 1 || !arguments.front().isArray()
        || containsFiniteApproximation(arguments.front().asArray()))
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
    if (arguments.size() != 1 || !arguments.front().isArray()
        || containsFiniteApproximation(arguments.front().asArray()))
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
    if (arguments.size() != 1 || !arguments.front().isArray()
        || containsFiniteApproximation(arguments.front().asArray()))
        return std::nullopt;
    return linear_algebra::approximateSingularValueDecomposition(
        arguments.front().asArray(), registry, mathematics, angles, context);
}

std::optional<Expr> evaluateApproximateConditionNumber(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray()
        || !arguments.front().asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments.front().asArray();
    if (containsFiniteApproximation(source))
        return std::nullopt;
    const std::size_t k = std::min(source.shape[0], source.shape[1]);
    if (k == 0)
        return std::nullopt;

    const auto rank = approximateRank(arguments, registry, mathematics, angles, context);
    if (!rank)
        return std::nullopt;
    if (*rank != k)
        return infinity();

    const auto decomposition = linear_algebra::approximateSingularValueDecomposition(
        source, registry, mathematics, angles, context);
    if (!decomposition)
        return std::nullopt;
    const auto factors = svdFactors(*decomposition);
    if (!factors || !(*factors)[1].isArray())
        return std::nullopt;
    const ArrayExpr& sigma = (*factors)[1].asArray();
    if (!sigma.isMatrix() || sigma.shape[0] != k || sigma.shape[1] != k)
        return std::nullopt;

    const Expr largest = sigma.element(0);
    const Expr smallest = sigma.element((k - 1) * k + (k - 1));
    return approximation::divideApproximateScalars(largest, smallest);
}

std::optional<Expr> evaluateApproximatePseudoInverse(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray()
        || !arguments.front().asArray().isMatrix()
        || containsFiniteApproximation(arguments.front().asArray()))
        return std::nullopt;
    return approximateSvdPseudoInverse(
        arguments.front().asArray(), registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateLeastSquares(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 2 || !arguments[0].isArray() || !arguments[1].isArray()
        || !arguments[0].asArray().isMatrix() || !arguments[1].asArray().isVector())
        return std::nullopt;
    if (arguments[1].asArray().shape[0] != arguments[0].asArray().shape[0])
        return std::nullopt;
    if (containsFiniteApproximation(arguments[0].asArray()))
        return std::nullopt;

    // 擬似逆行列を要求表示桁へ先に丸めると，その丸め幅をdotが再伝播して
    // leastSquaresだけ有効桁を余計に失う。中間値は作業桁まで保持し，
    // 利用者向けの丸めは最後のdotで一度だけ行う。
    approximation::ApproximationContext intermediateContext = context;
    intermediateContext.setDecimalDigits(context.workingDecimalDigits());
    const std::array<Expr, 1> inverseArguments{arguments[0]};
    const auto inverse = evaluateApproximatePseudoInverse(
        inverseArguments, registry, mathematics, angles, intermediateContext);
    if (!inverse)
        return std::nullopt;
    const std::array<Expr, 2> dotArguments{*inverse, arguments[1]};
    return evaluateApproximateDot(
        dotArguments, registry, mathematics, angles, std::move(context));
}

std::optional<Expr> evaluateApproximateEigenvalues(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments.front().isArray()
        || containsFiniteApproximation(arguments.front().asArray()))
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
    if (arguments.size() != 1 || !arguments.front().isArray()
        || containsFiniteApproximation(arguments.front().asArray()))
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
    if (arguments.size() != 1 || !arguments.front().isArray()
        || containsFiniteApproximation(arguments.front().asArray()))
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
