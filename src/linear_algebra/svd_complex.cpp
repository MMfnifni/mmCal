// 複素SVD。A^H Aを形成せず，Householder bidiagonalization + one-sided Jacobiを使う。
#include "svd_complex.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_sqrt.hpp"
#include "approximation/expression_interval.hpp"
#include "approximation/interval_math.hpp"
#include "approximation/precision.hpp"
#include "expression/array_utils.hpp"
#include "linear_algebra/complex_point.hpp"
#include "linear_algebra/point_arithmetic.hpp"
#include "numeric/big_int.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra::detail {
namespace {

using approximation::ComplexInterval;
using approximation::RealInterval;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;
using numeric::RoundingMode;

[[nodiscard]] BigFloat squareRoot(const BigFloat& value, std::size_t bits) {
    if (value.isNegative())
        throw approximation::PrecisionInsufficient("SVD encountered a negative rounded norm");
    const auto enclosure = approximation::encloseSqrt(value.toRational(), bits + 8).interval;
    return midpoint(enclosure, bits);
}

class ComplexPointMatrix final {
public:
    ComplexPointMatrix(std::size_t rows, std::size_t columns, std::size_t bits)
        : rows_(rows), columns_(columns), values_(rows * columns, complexZero(bits)) {}
    ComplexPointMatrix(std::size_t rows, std::size_t columns, std::vector<ComplexPoint> values)
        : rows_(rows), columns_(columns), values_(std::move(values)) {
        if (values_.size() != rows_ * columns_)
            throw std::invalid_argument("ComplexPointMatrix shape does not match element count");
    }
    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] ComplexPoint& operator()(std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const ComplexPoint& operator()(std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const std::vector<ComplexPoint>& values() const noexcept { return values_; }
private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<ComplexPoint> values_;
};

[[nodiscard]] ComplexPointMatrix identity(std::size_t rows, std::size_t columns, std::size_t bits) {
    ComplexPointMatrix result{rows, columns, bits};
    for (std::size_t i = 0; i < std::min(rows, columns); ++i)
        result(i, i) = complexOne(bits);
    return result;
}

[[nodiscard]] ComplexPointMatrix conjugateTranspose(
    const ComplexPointMatrix& matrix,
    std::size_t bits) {
    ComplexPointMatrix result{matrix.columns(), matrix.rows(), bits};
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        for (std::size_t column = 0; column < matrix.columns(); ++column)
            result(column, row) = complexConjugate(matrix(row, column));
    return result;
}

[[nodiscard]] ComplexPointMatrix multiplyMatrices(
    const ComplexPointMatrix& lhs,
    const ComplexPointMatrix& rhs,
    std::size_t bits) {
    if (lhs.columns() != rhs.rows())
        throw std::invalid_argument("SVD internal complex matrix dimensions are incompatible");
    ComplexPointMatrix result{lhs.rows(), rhs.columns(), bits};
    for (std::size_t row = 0; row < lhs.rows(); ++row)
        for (std::size_t column = 0; column < rhs.columns(); ++column) {
            ComplexPoint sum = complexZero(bits);
            for (std::size_t k = 0; k < lhs.columns(); ++k)
                sum = complexAdd(sum,
                    complexMultiply(lhs(row, k), rhs(k, column), bits), bits);
            result(row, column) = std::move(sum);
        }
    return result;
}

struct ComplexReflector final {
    std::size_t first = 0;
    std::vector<ComplexPoint> vector;
    BigFloat beta;
};

[[nodiscard]] std::optional<ComplexReflector> reflectorForValues(
    std::vector<ComplexPoint> values,
    std::size_t first,
    std::size_t bits) {
    BigFloat normSquared = zero(bits);
    for (const ComplexPoint& value : values)
        normSquared = add(normSquared, complexMagnitudeSquared(value, bits), bits);
    if (normSquared.isZero())
        return std::nullopt;

    const BigFloat norm = squareRoot(normSquared, bits);
    ComplexPoint phase = complexOne(bits);
    const BigFloat firstMagnitudeSquared = complexMagnitudeSquared(values.front(), bits);
    if (!firstMagnitudeSquared.isZero())
        phase = complexDivideReal(values.front(), squareRoot(firstMagnitudeSquared, bits), bits);
    const ComplexPoint alpha = complexNegate(complexScale(phase, norm, bits));
    values.front() = complexSubtract(values.front(), alpha, bits);

    BigFloat vNormSquared = zero(bits);
    for (const ComplexPoint& value : values)
        vNormSquared = add(vNormSquared, complexMagnitudeSquared(value, bits), bits);
    if (vNormSquared.isZero())
        return std::nullopt;
    return ComplexReflector{first, std::move(values), divide(two(bits), vNormSquared, bits)};
}

void applyHouseholderLeft(
    ComplexPointMatrix& matrix,
    const ComplexReflector& reflector,
    std::size_t firstColumn,
    std::size_t bits) {
    for (std::size_t column = firstColumn; column < matrix.columns(); ++column) {
        ComplexPoint dot = complexZero(bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            dot = complexAdd(dot,
                complexMultiply(complexConjugate(reflector.vector[i]),
                    matrix(reflector.first + i, column), bits), bits);
        const ComplexPoint scaled = complexScale(dot, reflector.beta, bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            matrix(reflector.first + i, column) = complexSubtract(
                matrix(reflector.first + i, column),
                complexMultiply(reflector.vector[i], scaled, bits), bits);
    }
}

void applyHouseholderRight(
    ComplexPointMatrix& matrix,
    const ComplexReflector& reflector,
    std::size_t firstRow,
    std::size_t bits) {
    for (std::size_t row = firstRow; row < matrix.rows(); ++row) {
        ComplexPoint dot = complexZero(bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            dot = complexAdd(dot,
                complexMultiply(matrix(row, reflector.first + i), reflector.vector[i], bits), bits);
        const ComplexPoint scaled = complexScale(dot, reflector.beta, bits);
        for (std::size_t i = 0; i < reflector.vector.size(); ++i)
            matrix(row, reflector.first + i) = complexSubtract(
                matrix(row, reflector.first + i),
                complexMultiply(scaled, complexConjugate(reflector.vector[i]), bits), bits);
    }
}

void bidiagonalize(
    ComplexPointMatrix& matrix,
    ComplexPointMatrix& left,
    ComplexPointMatrix& right,
    std::size_t bits) {
    const std::size_t kMax = std::min(matrix.rows(), matrix.columns());
    for (std::size_t k = 0; k < kMax; ++k) {
        std::vector<ComplexPoint> column;
        column.reserve(matrix.rows() - k);
        for (std::size_t row = k; row < matrix.rows(); ++row)
            column.push_back(matrix(row, k));
        if (const auto reflector = reflectorForValues(std::move(column), k, bits)) {
            applyHouseholderLeft(matrix, *reflector, k, bits);
            applyHouseholderRight(left, *reflector, 0, bits);
            for (std::size_t row = k + 1; row < matrix.rows(); ++row)
                matrix(row, k) = complexZero(bits);
        }

        if (k + 1 >= matrix.columns())
            continue;
        std::vector<ComplexPoint> rowValues;
        rowValues.reserve(matrix.columns() - k - 1);
        for (std::size_t columnIndex = k + 1; columnIndex < matrix.columns(); ++columnIndex)
            rowValues.push_back(complexConjugate(matrix(k, columnIndex)));
        if (const auto reflector = reflectorForValues(std::move(rowValues), k + 1, bits)) {
            // reflectorForValuesは列vectorに対するHを作るため，rowへ適用する場合は
            // conjugate(row)から作ったreflectorをそのまま右適用する。
            applyHouseholderRight(matrix, *reflector, k, bits);
            applyHouseholderRight(right, *reflector, 0, bits);
            for (std::size_t columnIndex = k + 2; columnIndex < matrix.columns(); ++columnIndex)
                matrix(k, columnIndex) = complexZero(bits);
        }
    }
}

[[nodiscard]] ComplexPoint columnDot(
    const ComplexPointMatrix& matrix,
    std::size_t lhs,
    std::size_t rhs,
    std::size_t bits) {
    ComplexPoint result = complexZero(bits);
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        result = complexAdd(result,
            complexMultiply(complexConjugate(matrix(row, lhs)), matrix(row, rhs), bits), bits);
    return result;
}

[[nodiscard]] BigFloat columnNormSquared(
    const ComplexPointMatrix& matrix,
    std::size_t column,
    std::size_t bits) {
    BigFloat result = zero(bits);
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        result = add(result, complexMagnitudeSquared(matrix(row, column), bits), bits);
    return result;
}

[[nodiscard]] bool sufficientlyOrthogonal(
    const BigFloat& alpha,
    const BigFloat& beta,
    const ComplexPoint& gamma,
    std::size_t thresholdBits,
    std::size_t bits) {
    if (complexIsZero(gamma) || alpha.isZero() || beta.isZero())
        return true;
    const BigFloat epsilon = BigFloat::fromDyadic(
        BigInt{1}, -static_cast<BigFloat::exponent_type>(thresholdBits), bits,
        RoundingMode::NearestEven);
    const BigFloat lhs = complexMagnitudeSquared(gamma, bits);
    const BigFloat epsSq = multiply(epsilon, epsilon, bits);
    const BigFloat rhs = multiply(epsSq, multiply(alpha, beta, bits), bits);
    return lhs <= rhs;
}

void scaleColumnPhase(
    ComplexPointMatrix& matrix,
    std::size_t column,
    const ComplexPoint& factor,
    std::size_t bits) {
    for (std::size_t row = 0; row < matrix.rows(); ++row)
        matrix(row, column) = complexMultiply(matrix(row, column), factor, bits);
}

void rotateColumns(
    ComplexPointMatrix& matrix,
    std::size_t p,
    std::size_t q,
    const BigFloat& c,
    const BigFloat& s,
    std::size_t bits) {
    for (std::size_t row = 0; row < matrix.rows(); ++row) {
        const ComplexPoint left = matrix(row, p);
        const ComplexPoint right = matrix(row, q);
        matrix(row, p) = complexSubtract(
            complexScale(left, c, bits), complexScale(right, s, bits), bits);
        matrix(row, q) = complexAdd(
            complexScale(left, s, bits), complexScale(right, c, bits), bits);
    }
}

void oneSidedJacobi(
    ComplexPointMatrix& matrix,
    ComplexPointMatrix& right,
    std::size_t targetBits,
    std::size_t bits) {
    const std::size_t columns = matrix.columns();
    if (columns < 2)
        return;
    const std::size_t thresholdBits = std::min<std::size_t>(
        static_cast<std::size_t>(std::numeric_limits<BigFloat::exponent_type>::max() / 2),
        std::max<std::size_t>(16, targetBits));
    const std::size_t maximumSweeps = std::max<std::size_t>(32, 10 * columns + 20);

    for (std::size_t sweep = 0; sweep < maximumSweeps; ++sweep) {
        bool changed = false;
        for (std::size_t p = 0; p + 1 < columns; ++p)
            for (std::size_t q = p + 1; q < columns; ++q) {
                const BigFloat alpha = columnNormSquared(matrix, p, bits);
                const BigFloat beta = columnNormSquared(matrix, q, bits);
                ComplexPoint gamma = columnDot(matrix, p, q, bits);
                if (sufficientlyOrthogonal(alpha, beta, gamma, thresholdBits, bits))
                    continue;
                changed = true;

                const BigFloat gammaMagnitude = squareRoot(complexMagnitudeSquared(gamma, bits), bits);
                if (gammaMagnitude.isZero())
                    continue;
                const ComplexPoint phase = complexDivideReal(gamma, gammaMagnitude, bits);
                const ComplexPoint phaseConjugate = complexConjugate(phase);
                // qへphase^-1を掛ければp^H qが正実数になる。Vにも同じunitary変換を積む。
                scaleColumnPhase(matrix, q, phaseConjugate, bits);
                scaleColumnPhase(right, q, phaseConjugate, bits);

                const BigFloat denominator = multiply(two(bits), gammaMagnitude, bits);
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
    throw approximation::PrecisionInsufficient(
        "Complex SVD Jacobi iteration did not converge at the current precision");
}

void completeOrthonormalColumns(ComplexPointMatrix& u, std::size_t bits) {
    for (std::size_t column = 0; column < u.columns(); ++column) {
        if (!columnNormSquared(u, column, bits).isZero())
            continue;
        bool found = false;
        for (std::size_t candidate = 0; candidate < u.rows() && !found; ++candidate) {
            std::vector<ComplexPoint> v(u.rows(), complexZero(bits));
            v[candidate] = complexOne(bits);
            for (int pass = 0; pass < 2; ++pass)
                for (std::size_t previous = 0; previous < column; ++previous) {
                    ComplexPoint coefficient = complexZero(bits);
                    for (std::size_t row = 0; row < u.rows(); ++row)
                        coefficient = complexAdd(coefficient,
                            complexMultiply(complexConjugate(u(row, previous)), v[row], bits), bits);
                    for (std::size_t row = 0; row < u.rows(); ++row)
                        v[row] = complexSubtract(v[row],
                            complexMultiply(u(row, previous), coefficient, bits), bits);
                }
            BigFloat normSquared = zero(bits);
            for (const ComplexPoint& value : v)
                normSquared = add(normSquared, complexMagnitudeSquared(value, bits), bits);
            if (normSquared.isZero())
                continue;
            const BigFloat norm = squareRoot(normSquared, bits);
            for (std::size_t row = 0; row < u.rows(); ++row)
                u(row, column) = complexDivideReal(v[row], norm, bits);
            found = true;
        }
        if (!found)
            throw approximation::PrecisionInsufficient(
                "Complex SVD could not complete an orthonormal null-space basis");
    }
}

struct ComplexPointSvd final {
    ComplexPointMatrix u;
    ComplexPointMatrix s;
    ComplexPointMatrix v;
};

[[nodiscard]] ComplexPointSvd svdTall(
    ComplexPointMatrix source,
    std::size_t targetBits,
    std::size_t bits) {
    const std::size_t rows = source.rows();
    const std::size_t columns = source.columns();
    if (rows < columns)
        throw std::invalid_argument("complex svdTall requires rows >= columns");

    ComplexPointMatrix left = identity(rows, rows, bits);
    ComplexPointMatrix rightBidiag = identity(columns, columns, bits);
    bidiagonalize(source, left, rightBidiag, bits);

    ComplexPointMatrix rightJacobi = identity(columns, columns, bits);
    oneSidedJacobi(source, rightJacobi, targetBits, bits);

    std::vector<BigFloat> singularValues;
    singularValues.reserve(columns);
    for (std::size_t column = 0; column < columns; ++column)
        singularValues.push_back(squareRoot(columnNormSquared(source, column, bits), bits));

    std::vector<std::size_t> order(columns);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        return singularValues[lhs] > singularValues[rhs];
    });

    ComplexPointMatrix uLocal{rows, columns, bits};
    ComplexPointMatrix rightSorted{columns, columns, bits};
    ComplexPointMatrix sigma{columns, columns, bits};
    for (std::size_t outputColumn = 0; outputColumn < columns; ++outputColumn) {
        const std::size_t sourceColumn = order[outputColumn];
        const BigFloat singular = singularValues[sourceColumn];
        sigma(outputColumn, outputColumn) = ComplexPoint{singular, zero(bits)};
        for (std::size_t row = 0; row < columns; ++row)
            rightSorted(row, outputColumn) = rightJacobi(row, sourceColumn);
        if (!singular.isZero())
            for (std::size_t row = 0; row < rows; ++row)
                uLocal(row, outputColumn) = complexDivideReal(
                    source(row, sourceColumn), singular, bits);
    }
    completeOrthonormalColumns(uLocal, bits);

    return ComplexPointSvd{
        multiplyMatrices(left, uLocal, bits),
        std::move(sigma),
        multiplyMatrices(rightBidiag, rightSorted, bits)
    };
}

[[nodiscard]] ComplexPointSvd computePointSvd(
    const ComplexPointMatrix& source,
    std::size_t targetBits,
    std::size_t bits) {
    if (source.rows() >= source.columns())
        return svdTall(source, targetBits, bits);
    ComplexPointSvd transposed = svdTall(conjugateTranspose(source, bits), targetBits, bits);
    return ComplexPointSvd{
        std::move(transposed.v),
        std::move(transposed.s),
        std::move(transposed.u)
    };
}

[[nodiscard]] std::optional<std::vector<ComplexInterval>> encloseElements(
    const ArrayExpr& source,
    std::size_t bits,
    const approximation::CertifiedEvaluator& certified) {
    std::vector<ComplexInterval> result;
    result.reserve(source.size());
    for (std::size_t i = 0; i < source.size(); ++i) {
        const Expr element = source.element(i);
        const auto enclosed = approximation::encloseComplexExpression(element, bits, certified);
        if (!enclosed)
            return std::nullopt;
        result.push_back(*enclosed);
    }
    return result;
}

[[nodiscard]] ComplexPointMatrix midpointMatrix(
    std::size_t rows,
    std::size_t columns,
    const std::vector<ComplexInterval>& values,
    std::size_t bits) {
    std::vector<ComplexPoint> points;
    points.reserve(values.size());
    for (const ComplexInterval& value : values)
        points.push_back(ComplexPoint{
            midpoint(value.real(), bits), midpoint(value.imaginary(), bits)});
    return ComplexPointMatrix{rows, columns, std::move(points)};
}

class IntervalComplexMatrix final {
public:
    IntervalComplexMatrix(std::size_t rows, std::size_t columns, std::vector<ComplexInterval> values)
        : rows_(rows), columns_(columns), values_(std::move(values)) {}
    explicit IntervalComplexMatrix(const ComplexPointMatrix& source)
        : rows_(source.rows()), columns_(source.columns()) {
        values_.reserve(rows_ * columns_);
        for (const ComplexPoint& value : source.values())
            values_.emplace_back(
                RealInterval::point(value.real), RealInterval::point(value.imaginary));
    }
    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] const ComplexInterval& operator()(std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }
private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<ComplexInterval> values_;
};

[[nodiscard]] ComplexInterval conjugate(const ComplexInterval& value) {
    return ComplexInterval{value.real(), approximation::negate(value.imaginary())};
}

[[nodiscard]] IntervalComplexMatrix conjugateTranspose(
    const IntervalComplexMatrix& matrix) {
    std::vector<ComplexInterval> values;
    values.reserve(matrix.rows() * matrix.columns());
    for (std::size_t row = 0; row < matrix.columns(); ++row)
        for (std::size_t column = 0; column < matrix.rows(); ++column)
            values.push_back(conjugate(matrix(column, row)));
    return IntervalComplexMatrix{matrix.columns(), matrix.rows(), std::move(values)};
}

[[nodiscard]] IntervalComplexMatrix intervalMultiply(
    const IntervalComplexMatrix& lhs,
    const IntervalComplexMatrix& rhs,
    std::size_t bits) {
    std::vector<ComplexInterval> values;
    values.reserve(lhs.rows() * rhs.columns());
    const ComplexInterval zeroValue = ComplexInterval::fromReal(
        RealInterval::fromRational(Rational{}, bits));
    for (std::size_t row = 0; row < lhs.rows(); ++row)
        for (std::size_t column = 0; column < rhs.columns(); ++column) {
            ComplexInterval sum = zeroValue;
            for (std::size_t k = 0; k < lhs.columns(); ++k)
                sum = approximation::add(sum,
                    approximation::multiply(lhs(row, k), rhs(k, column), bits), bits);
            values.push_back(std::move(sum));
        }
    return IntervalComplexMatrix{lhs.rows(), rhs.columns(), std::move(values)};
}

[[nodiscard]] BigFloat squaredMagnitudeUpper(
    const ComplexInterval& value,
    std::size_t bits) {
    return approximation::add(
        approximation::squareInterval(value.real(), bits),
        approximation::squareInterval(value.imaginary(), bits), bits).upper();
}

[[nodiscard]] BigFloat maxResidualSquaredUpper(
    const IntervalComplexMatrix& lhs,
    const IntervalComplexMatrix& rhs,
    std::size_t bits) {
    BigFloat result = zero(bits);
    for (std::size_t row = 0; row < lhs.rows(); ++row)
        for (std::size_t column = 0; column < lhs.columns(); ++column) {
            const BigFloat candidate = squaredMagnitudeUpper(
                approximation::subtract(lhs(row, column), rhs(row, column), bits), bits);
            if (candidate > result)
                result = candidate;
        }
    return result;
}

[[nodiscard]] BigFloat orthogonalityResidualSquared(
    const ComplexPointMatrix& matrix,
    std::size_t bits) {
    const IntervalComplexMatrix value{matrix};
    const IntervalComplexMatrix gram = intervalMultiply(
        conjugateTranspose(value), value, bits);
    BigFloat result = zero(bits);
    for (std::size_t row = 0; row < gram.rows(); ++row)
        for (std::size_t column = 0; column < gram.columns(); ++column) {
            const ComplexInterval target = ComplexInterval::fromReal(
                RealInterval::fromRational(Rational{BigInt{row == column ? 1 : 0}}, bits));
            const BigFloat candidate = squaredMagnitudeUpper(
                approximation::subtract(gram(row, column), target, bits), bits);
            if (candidate > result)
                result = candidate;
        }
    return result;
}

[[nodiscard]] BigFloat decimalTolerance(std::size_t digits, std::size_t bits) {
    if (digits > std::numeric_limits<std::uint64_t>::max() - 2)
        throw std::overflow_error("SVD requested decimal precision is too large");
    const BigInt denominator = numeric::pow(BigInt{10}, static_cast<std::uint64_t>(digits + 2));
    return BigFloat::fromRational(
        Rational{BigInt{1}, denominator}, bits, RoundingMode::TowardPositive);
}

[[nodiscard]] bool verifyPointSvd(
    const std::vector<ComplexInterval>& source,
    std::size_t rows,
    std::size_t columns,
    const ComplexPointSvd& result,
    std::size_t digits,
    std::size_t bits) {
    const IntervalComplexMatrix a{rows, columns, source};
    const IntervalComplexMatrix u{result.u};
    const IntervalComplexMatrix s{result.s};
    const IntervalComplexMatrix v{result.v};
    const IntervalComplexMatrix reconstructed = intervalMultiply(
        intervalMultiply(u, s, bits), conjugateTranspose(v), bits);
    const BigFloat residualSquared = maxResidualSquaredUpper(a, reconstructed, bits);
    const BigFloat uResidualSquared = orthogonalityResidualSquared(result.u, bits);
    const BigFloat vResidualSquared = orthogonalityResidualSquared(result.v, bits);
    const BigFloat tolerance = decimalTolerance(digits, bits);
    const BigFloat toleranceSquared = numeric::multiply(
        tolerance, tolerance, bits, RoundingMode::TowardPositive);
    return residualSquared <= toleranceSquared
        && uResidualSquared <= toleranceSquared
        && vResidualSquared <= toleranceSquared;
}

[[nodiscard]] std::optional<Expr> decimalPointMatrix(
    const ComplexPointMatrix& matrix,
    std::size_t digits) {
    std::vector<Expr> values;
    values.reserve(matrix.values().size());
    for (const ComplexPoint& value : matrix.values()) {
        const Rational realExact = value.real.toRational();
        const Rational imaginaryExact = value.imaginary.toRational();
        const auto real = numeric::DecimalApproximation::fromCertifiedIntervalSignificant(
            realExact, realExact, digits);
        const auto imaginary = numeric::DecimalApproximation::fromCertifiedIntervalSignificant(
            imaginaryExact, imaginaryExact, digits);
        if (!real || !imaginary)
            return std::nullopt;
        if (imaginaryExact.isZero())
            values.emplace_back(*real);
        else
            values.emplace_back(numeric::ComplexDecimalApproximation::fromComponents(
                *real, *imaginary, realExact.isZero(), false));
    }
    return Expr::array({matrix.rows(), matrix.columns()}, std::move(values));
}

} // namespace

std::optional<Expr> approximateComplexSvdAtPrecision(
    const ArrayExpr& source,
    std::size_t bits,
    std::size_t digits,
    const approximation::CertifiedEvaluator& certified) {
    const auto enclosed = encloseElements(source, bits, certified);
    if (!enclosed)
        return std::nullopt;
    const ComplexPointMatrix midpointSource = midpointMatrix(
        source.shape[0], source.shape[1], *enclosed, bits);
    ComplexPointSvd result = computePointSvd(
        midpointSource, approximation::decimalDigitsToBinaryBits(digits + 2), bits);
    if (!verifyPointSvd(*enclosed, source.shape[0], source.shape[1], result, digits, bits))
        throw approximation::PrecisionInsufficient(
            "Complex SVD reconstruction or orthogonality is not certified at the current precision");

    const auto u = decimalPointMatrix(result.u, digits);
    const auto s = decimalPointMatrix(result.s, digits);
    const auto v = decimalPointMatrix(result.v, digits);
    if (!u || !s || !v)
        return std::nullopt;
    return expression::braceValue({*u, *s, *v});
}

} // namespace mmcal::linear_algebra::detail
