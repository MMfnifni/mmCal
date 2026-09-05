// reduced SVD。exactな自明caseと，BigFloat Householder + one-sided Jacobi numerical backend。
#include "svd.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certified_sqrt.hpp"
#include "approximation/expression_interval.hpp"
#include "approximation/interval_math.hpp"
#include "approximation/precision.hpp"
#include "expression/array_utils.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "linear_algebra/point_arithmetic.hpp"
#include "numeric/big_float.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "svd_complex.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra {
namespace {

using approximation::ComplexInterval;
using approximation::RealInterval;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigFloat;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RoundingMode;

constexpr std::size_t maximumPrecisionRetries = 10;

using detail::absolute;
using detail::add;
using detail::divide;
using detail::midpoint;
using detail::multiply;
using detail::one;
using detail::subtract;
using detail::two;
using detail::zero;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] BigFloat squareRoot(const BigFloat& value, std::size_t bits) {
    if (value.isNegative())
        throw approximation::PrecisionInsufficient("SVD encountered a negative rounded norm");
    const auto enclosure = approximation::encloseSqrt(value.toRational(), bits + 8).interval;
    return midpoint(enclosure, bits);
}

class PointMatrix final {
public:
    PointMatrix(std::size_t rows, std::size_t columns, std::size_t bits)
        : rows_(rows), columns_(columns), values_(rows * columns, zero(bits)) {}

    PointMatrix(std::size_t rows, std::size_t columns, std::vector<BigFloat> values)
        : rows_(rows), columns_(columns), values_(std::move(values)) {
        if (values_.size() != rows_ * columns_)
            throw std::invalid_argument("PointMatrix shape does not match element count");
    }

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] BigFloat& operator()(std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const BigFloat& operator()(std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const std::vector<BigFloat>& values() const noexcept { return values_; }

private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<BigFloat> values_;
};

[[nodiscard]] PointMatrix identity(std::size_t rows, std::size_t columns, std::size_t bits) {
    PointMatrix result{rows, columns, bits};
    for (std::size_t i = 0; i < std::min(rows, columns); ++i)
        result(i, i) = one(bits);
    return result;
}

[[nodiscard]] PointMatrix transpose(const PointMatrix& matrix, std::size_t bits) {
    PointMatrix result{matrix.columns(), matrix.rows(), bits};
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < matrix.columns(); ++column)
            result(column, row) = matrix(row, column);
    return result;
}

[[nodiscard]] PointMatrix multiplyMatrices(
    const PointMatrix& lhs,
    const PointMatrix& rhs,
    std::size_t bits) {
    if (lhs.columns() != rhs.rows())
        throw std::invalid_argument("SVD internal matrix dimensions are incompatible");
    PointMatrix result{lhs.rows(), rhs.columns(), bits};
    for (std::size_t row = 0; row < lhs.rows(); ++row)
        for (std::size_t column = 0; column < rhs.columns(); ++column) {
            BigFloat sum = zero(bits);
            for (std::size_t k = 0; k < lhs.columns(); ++k)
                sum = add(sum, multiply(lhs(row, k), rhs(k, column), bits), bits);
            result(row, column) = std::move(sum);
        }
    return result;
}

struct PointReflector final {
    std::size_t first = 0;
    std::vector<BigFloat> vector;
    BigFloat beta;
};

[[nodiscard]] std::optional<PointReflector> reflectorForValues(
    std::vector<BigFloat> values,
    std::size_t first,
    std::size_t bits) {
    BigFloat normSquared = zero(bits);
    for (const BigFloat& value : values)
        normSquared = add(normSquared, multiply(value, value, bits), bits);
    if (normSquared.isZero())
        return std::nullopt;

    const BigFloat norm = squareRoot(normSquared, bits);
    const BigFloat alpha = values.front().isNegative() ? norm : -norm;
    values.front() = subtract(values.front(), alpha, bits);

    BigFloat vNormSquared = zero(bits);
    for (const BigFloat& value : values)
        vNormSquared = add(vNormSquared, multiply(value, value, bits), bits);
    if (vNormSquared.isZero())
        return std::nullopt;
    return PointReflector{first, std::move(values), divide(two(bits), vNormSquared, bits)};
}

void applyHouseholderLeft(
    PointMatrix& matrix,
    const PointReflector& reflector,
    std::size_t firstColumn,
    std::size_t bits) {
    for (std::size_t column = firstColumn; column < matrix.columns(); ++column) {
        BigFloat dot = zero(bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            dot = add(dot, multiply(reflector.vector[i],
                matrix(reflector.first + i, column), bits), bits);
        const BigFloat scaled = multiply(reflector.beta, dot, bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            matrix(reflector.first + i, column) = subtract(
                matrix(reflector.first + i, column),
                multiply(reflector.vector[i], scaled, bits), bits);
    }
}

void applyHouseholderRight(
    PointMatrix& matrix,
    const PointReflector& reflector,
    std::size_t firstRow,
    std::size_t bits) {
    for (std::size_t row = firstRow; row < matrix.rows(); ++row) {
        BigFloat dot = zero(bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            dot = add(dot, multiply(matrix(row, reflector.first + i),
                reflector.vector[i], bits), bits);
        const BigFloat scaled = multiply(reflector.beta, dot, bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            matrix(row, reflector.first + i) = subtract(
                matrix(row, reflector.first + i),
                multiply(scaled, reflector.vector[i], bits), bits);
    }
}

void bidiagonalize(
    PointMatrix& matrix,
    PointMatrix& left,
    PointMatrix& right,
    std::size_t bits) {
    const std::size_t kMax = std::min(matrix.rows(), matrix.columns());
    for (std::size_t k = 0; k < kMax; ++k) {
        std::vector<BigFloat> column;
        column.reserve(matrix.rows() - k);
        for (std::size_t row = k; row < matrix.rows(); ++row)
            column.push_back(matrix(row, k));
        if (const auto reflector = reflectorForValues(std::move(column), k, bits)) {
            applyHouseholderLeft(matrix, *reflector, k, bits);
            applyHouseholderRight(left, *reflector, 0, bits); // left <- left H
            for (std::size_t row = k + 1; row < matrix.rows(); ++row)
                matrix(row, k) = zero(bits);
        }

        if (k + 1 >= matrix.columns())
            continue;
        std::vector<BigFloat> rowValues;
        rowValues.reserve(matrix.columns() - k - 1);
        for (std::size_t columnIndex = k + 1; columnIndex < matrix.columns(); ++columnIndex)
            rowValues.push_back(matrix(k, columnIndex));
        if (const auto reflector = reflectorForValues(std::move(rowValues), k + 1, bits)) {
            applyHouseholderRight(matrix, *reflector, k, bits);
            applyHouseholderRight(right, *reflector, 0, bits); // right <- right G
            for (std::size_t columnIndex = k + 2; columnIndex < matrix.columns(); ++columnIndex)
                matrix(k, columnIndex) = zero(bits);
        }
    }
}

[[nodiscard]] BigFloat columnDot(
    const PointMatrix& matrix,
    std::size_t lhs,
    std::size_t rhs,
    std::size_t bits) {
    BigFloat result = zero(bits);
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        result = add(result,
            multiply(matrix(row, lhs), matrix(row, rhs), bits), bits);
    return result;
}

[[nodiscard]] bool sufficientlyOrthogonal(
    const BigFloat& alpha,
    const BigFloat& beta,
    const BigFloat& gamma,
    std::size_t thresholdBits,
    std::size_t bits) {
    if (gamma.isZero() || alpha.isZero() || beta.isZero())
        return true;
    const BigFloat epsilon = BigFloat::fromDyadic(
        BigInt{1}, -static_cast<BigFloat::exponent_type>(thresholdBits), bits,
        RoundingMode::NearestEven);
    const BigFloat lhs = multiply(gamma, gamma, bits);
    const BigFloat epsSq = multiply(epsilon, epsilon, bits);
    const BigFloat rhs = multiply(epsSq, multiply(alpha, beta, bits), bits);
    return lhs <= rhs;
}

void rotateColumns(
    PointMatrix& matrix,
    std::size_t p,
    std::size_t q,
    const BigFloat& c,
    const BigFloat& s,
    std::size_t bits) {
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        const BigFloat left = matrix(row, p);
        const BigFloat right = matrix(row, q);
        matrix(row, p) = subtract(multiply(c, left, bits), multiply(s, right, bits), bits);
        matrix(row, q) = add(multiply(s, left, bits), multiply(c, right, bits), bits);
    }
}

void oneSidedJacobi(
    PointMatrix& matrix,
    PointMatrix& right,
    std::size_t targetBits,
    std::size_t bits) {
    const std::size_t columns = matrix.columns();
    if (columns < 2)
        return;
    const std::size_t thresholdBits = std::min<std::size_t>(
        static_cast<std::size_t>(std::numeric_limits<BigFloat::exponent_type>::max() / 2),
        std::max<std::size_t>(16, targetBits));
    const std::size_t maximumSweeps = std::max<std::size_t>(32, 8 * columns + 16);

    for (std::size_t sweep = 0; sweep < maximumSweeps; ++sweep) {
        bool changed = false;
        for (std::size_t p = 0; p + 1 < columns; ++p)
            for (std::size_t q = p + 1; q < columns; ++q) {
                const BigFloat alpha = columnDot(matrix, p, p, bits);
                const BigFloat beta = columnDot(matrix, q, q, bits);
                const BigFloat gamma = columnDot(matrix, p, q, bits);
                if (sufficientlyOrthogonal(alpha, beta, gamma, thresholdBits, bits))
                    continue;
                changed = true;

                const BigFloat denominator = multiply(two(bits), gamma, bits);
                if (denominator.isZero())
                    continue;
                const BigFloat tau = divide(subtract(beta, alpha, bits), denominator, bits);
                const BigFloat root = squareRoot(add(one(bits), multiply(tau, tau, bits), bits), bits);
                const BigFloat tDenominator = add(absolute(tau), root, bits);
                BigFloat t = divide(one(bits), tDenominator, bits);
                if (tau.isNegative())
                    t = -t;
                const BigFloat c = divide(one(bits),
                    squareRoot(add(one(bits), multiply(t, t, bits), bits), bits), bits);
                const BigFloat s = multiply(c, t, bits);
                rotateColumns(matrix, p, q, c, s, bits);
                rotateColumns(right, p, q, c, s, bits);
            }
        if (!changed)
            return;
    }
    throw approximation::PrecisionInsufficient("SVD Jacobi iteration did not converge at the current precision");
}

[[nodiscard]] BigFloat columnNorm(
    const PointMatrix& matrix,
    std::size_t column,
    std::size_t bits) {
    return squareRoot(columnDot(matrix, column, column, bits), bits);
}

void completeOrthonormalColumns(PointMatrix& u, std::size_t bits) {
    // zero singular valueに対応するcolumnだけをdeterministicな標準基底から二回直交化する。
    for (std::size_t column = 0; column < u.columns(); ++column) {
        BigFloat normSq = columnDot(u, column, column, bits);
        if (!normSq.isZero())
            continue;

        bool found = false;
        for (std::size_t candidate = 0; candidate < u.rows() && !found; ++candidate) {
            std::vector<BigFloat> v(u.rows(), zero(bits));
            v[candidate] = one(bits);
            for (int pass = 0; pass < 2; ++pass)
                for (std::size_t previous = 0; previous < column; ++previous) {
                    BigFloat dot = zero(bits);
                    for (std::size_t row = 0; row < u.rows(); ++row)
                        dot = add(dot, multiply(v[row], u(row, previous), bits), bits);
                    for (std::size_t row = 0; row < u.rows(); ++row)
                        v[row] = subtract(v[row], multiply(dot, u(row, previous), bits), bits);
                }
            BigFloat candidateNormSq = zero(bits);
            for (const BigFloat& value : v)
                candidateNormSq = add(candidateNormSq, multiply(value, value, bits), bits);
            if (candidateNormSq.isZero())
                continue;
            const BigFloat norm = squareRoot(candidateNormSq, bits);
            for (std::size_t row = 0; row < u.rows(); ++row)
                u(row, column) = divide(v[row], norm, bits);
            found = true;
        }
        if (!found)
            throw approximation::PrecisionInsufficient("SVD could not complete an orthonormal null-space basis");
    }
}

struct PointSvd final {
    PointMatrix u;
    PointMatrix s;
    PointMatrix v;
};

[[nodiscard]] PointSvd svdTall(PointMatrix source, std::size_t targetBits, std::size_t bits) {
    const std::size_t rows = source.rows();
    const std::size_t columns = source.columns();
    if (rows < columns)
        throw std::invalid_argument("svdTall requires rows >= columns");

    PointMatrix left = identity(rows, rows, bits);
    PointMatrix rightBidiag = identity(columns, columns, bits);
    bidiagonalize(source, left, rightBidiag, bits);

    PointMatrix rightJacobi = identity(columns, columns, bits);
    oneSidedJacobi(source, rightJacobi, targetBits, bits);

    std::vector<BigFloat> singularValues;
    singularValues.reserve(columns);
    for (std::size_t column = 0; column < columns; ++column)
        singularValues.push_back(columnNorm(source, column, bits));

    std::vector<std::size_t> order(columns);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        return singularValues[lhs] > singularValues[rhs];
    });

    PointMatrix uLocal{rows, columns, bits};
    PointMatrix rightSorted{columns, columns, bits};
    PointMatrix sigma{columns, columns, bits};
    for (std::size_t outputColumn = 0; outputColumn < columns; ++outputColumn) {
        const std::size_t sourceColumn = order[outputColumn];
        const BigFloat singular = singularValues[sourceColumn];
        sigma(outputColumn, outputColumn) = singular;
        for (std::size_t row = 0; row < columns; ++row)
            rightSorted(row, outputColumn) = rightJacobi(row, sourceColumn);
        if (!singular.isZero())
            for (std::size_t row = 0; row < rows; ++row)
                uLocal(row, outputColumn) = divide(source(row, sourceColumn), singular, bits);
    }
    completeOrthonormalColumns(uLocal, bits);

    PointMatrix u = multiplyMatrices(left, uLocal, bits);
    PointMatrix v = multiplyMatrices(rightBidiag, rightSorted, bits);
    return PointSvd{std::move(u), std::move(sigma), std::move(v)};
}

[[nodiscard]] PointSvd computePointSvd(
    const PointMatrix& source,
    std::size_t targetBits,
    std::size_t bits) {
    if (source.rows() >= source.columns())
        return svdTall(source, targetBits, bits);

    PointSvd transposed = svdTall(transpose(source, bits), targetBits, bits);
    return PointSvd{
        std::move(transposed.v),
        std::move(transposed.s),
        std::move(transposed.u)
    };
}

[[nodiscard]] std::optional<std::vector<RealInterval>> encloseRealElements(
    const ArrayExpr& array,
    std::size_t bits,
    const approximation::CertifiedEvaluator& certified) {
    std::vector<RealInterval> values;
    values.reserve(array.size());
    for (std::size_t i = 0; i < array.size(); ++i) {
        const Expr element = array.element(i);
        const auto enclosed = approximation::encloseComplexExpression(element, bits, certified);
        if (!enclosed || !enclosed->isProvablyReal())
            return std::nullopt;
        values.push_back(enclosed->real());
    }
    return values;
}

[[nodiscard]] PointMatrix midpointMatrix(
    std::size_t rows,
    std::size_t columns,
    const std::vector<RealInterval>& values,
    std::size_t bits) {
    std::vector<BigFloat> points;
    points.reserve(values.size());
    for (const RealInterval& value : values)
        points.push_back(midpoint(value, bits));
    return PointMatrix{rows, columns, std::move(points)};
}

class IntervalRealMatrix final {
public:
    IntervalRealMatrix(std::size_t rows, std::size_t columns, std::size_t bits)
        : rows_(rows), columns_(columns), values_(rows * columns,
            RealInterval::fromRational(Rational{}, bits)) {}
    IntervalRealMatrix(std::size_t rows, std::size_t columns, std::vector<RealInterval> values)
        : rows_(rows), columns_(columns), values_(std::move(values)) {}
    explicit IntervalRealMatrix(const PointMatrix& source)
        : rows_(source.rows()), columns_(source.columns()) {
        values_.reserve(rows_ * columns_);
        for (const BigFloat& value : source.values())
            values_.push_back(RealInterval::point(value));
    }
    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] RealInterval& operator()(std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const RealInterval& operator()(std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }
private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<RealInterval> values_;
};

[[nodiscard]] IntervalRealMatrix intervalTranspose(
    const IntervalRealMatrix& matrix,
    std::size_t bits) {
    IntervalRealMatrix result{matrix.columns(), matrix.rows(), bits};
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < matrix.columns(); ++column)
            result(column, row) = matrix(row, column);
    return result;
}

[[nodiscard]] IntervalRealMatrix intervalMultiply(
    const IntervalRealMatrix& lhs,
    const IntervalRealMatrix& rhs,
    std::size_t bits) {
    IntervalRealMatrix result{lhs.rows(), rhs.columns(), bits};
    for (std::size_t row = 0; row < lhs.rows(); ++row)
        for (std::size_t column = 0; column < rhs.columns(); ++column) {
            RealInterval sum = RealInterval::fromRational(Rational{}, bits);
            for (std::size_t k = 0; k < lhs.columns(); ++k)
                sum = approximation::add(sum,
                    approximation::multiply(lhs(row, k), rhs(k, column), bits), bits);
            result(row, column) = std::move(sum);
        }
    return result;
}

[[nodiscard]] BigFloat maxResidualUpper(
    const IntervalRealMatrix& lhs,
    const IntervalRealMatrix& rhs,
    std::size_t bits) {
    BigFloat result = zero(bits);
    for (std::size_t row = 0; row < lhs.rows(); ++row)
        for (std::size_t column = 0; column < lhs.columns(); ++column) {
            const RealInterval difference = approximation::subtract(
                lhs(row, column), rhs(row, column), bits);
            const BigFloat upper = approximation::absoluteInterval(difference, bits).upper();
            if (upper > result)
                result = upper;
        }
    return result;
}

[[nodiscard]] BigFloat orthogonalityResidual(
    const PointMatrix& matrix,
    std::size_t bits) {
    const IntervalRealMatrix value{matrix};
    const IntervalRealMatrix gram = intervalMultiply(
        intervalTranspose(value, bits), value, bits);
    BigFloat result = zero(bits);
    for (std::size_t row = 0; row < gram.rows(); ++row)
        for (std::size_t column = 0; column < gram.columns(); ++column) {
            RealInterval target = RealInterval::fromRational(
                Rational{BigInt{row == column ? 1 : 0}}, bits);
            const BigFloat upper = approximation::absoluteInterval(
                approximation::subtract(gram(row, column), target, bits), bits).upper();
            if (upper > result)
                result = upper;
        }
    return result;
}

[[nodiscard]] BigFloat decimalTolerance(std::size_t digits, std::size_t bits) {
    if (digits > std::numeric_limits<std::uint64_t>::max() - 2)
        throw std::overflow_error("SVD requested decimal precision is too large");
    const BigInt denominator = numeric::pow(BigInt{10}, static_cast<std::uint64_t>(digits + 2));
    return BigFloat::fromRational(Rational{BigInt{1}, denominator}, bits, RoundingMode::TowardPositive);
}

[[nodiscard]] bool verifyPointSvd(
    const std::vector<RealInterval>& source,
    std::size_t rows,
    std::size_t columns,
    const PointSvd& result,
    std::size_t digits,
    std::size_t bits) {
    const IntervalRealMatrix a{rows, columns, source};
    const IntervalRealMatrix u{result.u};
    const IntervalRealMatrix s{result.s};
    const IntervalRealMatrix v{result.v};
    const IntervalRealMatrix reconstructed = intervalMultiply(
        intervalMultiply(u, s, bits), intervalTranspose(v, bits), bits);
    const BigFloat residual = maxResidualUpper(a, reconstructed, bits);
    const BigFloat uResidual = orthogonalityResidual(result.u, bits);
    const BigFloat vResidual = orthogonalityResidual(result.v, bits);
    const BigFloat tolerance = decimalTolerance(digits, bits);
    return residual <= tolerance && uResidual <= tolerance && vResidual <= tolerance;
}

[[nodiscard]] std::optional<Expr> decimalPointMatrix(
    const PointMatrix& matrix,
    std::size_t digits) {
    std::vector<Expr> values;
    values.reserve(matrix.values().size());
    for (const BigFloat& value : matrix.values()) {
        const Rational exact = value.toRational();
        const auto decimal = numeric::DecimalApproximation::fromVerifiedValueSignificant(
            numeric::RealNumber{exact}, digits);
        values.emplace_back(std::move(decimal));
    }
    return Expr::array({matrix.rows(), matrix.columns()}, std::move(values));
}

[[nodiscard]] std::optional<Expr> approximateAtPrecision(
    const ArrayExpr& source,
    std::size_t bits,
    std::size_t digits,
    const approximation::CertifiedEvaluator& certified) {
    const auto enclosed = encloseRealElements(source, bits, certified);
    if (!enclosed)
        return detail::approximateComplexSvdAtPrecision(
            source, bits, digits, certified);
    const PointMatrix midpointSource = midpointMatrix(
        source.shape[0], source.shape[1], *enclosed, bits);
    PointSvd result = computePointSvd(midpointSource,
        approximation::decimalDigitsToBinaryBits(digits + 2), bits);
    if (!verifyPointSvd(*enclosed, source.shape[0], source.shape[1], result, digits, bits))
        throw approximation::PrecisionInsufficient(
            "SVD reconstruction or orthogonality is not certified at the current precision");

    const auto u = decimalPointMatrix(result.u, digits);
    const auto s = decimalPointMatrix(result.s, digits);
    const auto v = decimalPointMatrix(result.v, digits);
    if (!u || !s || !v)
        return std::nullopt;
    return expression::braceValue({*u, *s, *v});
}

[[nodiscard]] bool exactRealDiagonal(const MatrixView& matrix) {
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < matrix.columns(); ++column) {
            const Expr value = matrix(row, column);
            if (!value.isNumber() || !value.asNumber().isReal())
                return false;
            if (row != column && !value.asNumber().isZero())
                return false;
        }
    return true;
}

[[nodiscard]] Expr exactDiagonalSvd(const MatrixView& matrix) {
    const std::size_t k = std::min(matrix.rows(), matrix.columns());
    std::vector<Expr> u(matrix.rows() * k, integer(0));
    std::vector<Expr> s(k * k, integer(0));
    std::vector<Expr> v(matrix.columns() * k, integer(0));
    std::vector<std::size_t> order(k);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        const auto lhsAbs = matrix(lhs, lhs).asNumber().asReal().abs();
        const auto rhsAbs = matrix(rhs, rhs).asNumber().asReal().abs();
        return lhsAbs > rhsAbs;
    });

    for (std::size_t outputColumn = 0; outputColumn < k; ++outputColumn) {
        const std::size_t sourceColumn = order[outputColumn];
        const Number diagonal = matrix(sourceColumn, sourceColumn).asNumber();
        const bool negative = diagonal.asReal().isNegative();
        Number singular = negative ? -diagonal : diagonal;
        Number sign = negative ? Number{BigInt{-1}} : Number{BigInt{1}};
        u[sourceColumn * k + outputColumn] = Expr{std::move(sign)};
        s[outputColumn * k + outputColumn] = Expr{std::move(singular)};
        v[sourceColumn * k + outputColumn] = integer(1);
    }
    const Expr uExpr = Expr::array({matrix.rows(), k}, std::move(u));
    const Expr sExpr = Expr::array({k, k}, std::move(s));
    const Expr vExpr = Expr::array({matrix.columns(), k}, std::move(v));
    return expression::braceValue({uExpr, sExpr, vExpr});
}

} // namespace

std::optional<Expr> singularValueDecomposition(
    const MatrixView& matrix,
    const ExactMatrixContext&) {
    if (exactRealDiagonal(matrix))
        return exactDiagonalSvd(matrix);
    return std::nullopt;
}

std::optional<Expr> approximateSingularValueDecomposition(
    const ArrayExpr& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (!matrix.isMatrix())
        return std::nullopt;
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    // PrecisionInsufficientはguard桁を増やせば解消し得るが，BackendUnsupportedは構造的な未対応である。
    // 後者を同じ入力で再試行しても改善しないため，直ちに上位fallbackへ返す。
    for (std::size_t attempt = 0; attempt < maximumPrecisionRetries; ++attempt) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        try {
            if (const auto result = approximateAtPrecision(
                matrix, context.workingBinaryBits(), context.decimalDigits(), certified))
                return result;
        }
        catch (const approximation::PrecisionInsufficient&) {
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return std::nullopt;
        }
        context.setGuardDigits(approximation::nextGuardDigits(context.guardDigits()));
    }
    return std::nullopt;
}

} // namespace mmcal::linear_algebra
