// precision-aware certified Matrix backend
#include "approximate_matrix.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_evaluator.hpp"
#include "approximation/certified_sqrt.hpp"
#include "approximation/expression_interval.hpp"
#include "approximation/interval_math.hpp"
#include "builtins/array_helpers.hpp"
#include "error/error_message.hpp"
#include "expression/array_utils.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra {
namespace {

using approximation::ComplexInterval;
using approximation::RealInterval;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigInt;
using numeric::Rational;

constexpr std::size_t maximumPrecisionRetries = 12;

[[nodiscard]] ComplexInterval exactComplex(
    std::int64_t real,
    std::size_t precisionBits) {
    return ComplexInterval::fromReal(RealInterval::fromRational(
        Rational{BigInt{real}}, precisionBits));
}

[[nodiscard]] bool exactZero(const RealInterval& value) noexcept {
    return value.isPoint() && value.lower().isZero();
}

[[nodiscard]] bool exactZero(const ComplexInterval& value) noexcept {
    return exactZero(value.real()) && exactZero(value.imaginary());
}

[[nodiscard]] std::size_t augmentedColumns(std::size_t columns) {
    if (columns > std::numeric_limits<std::size_t>::max() / 2)
        throw std::length_error("Augmented matrix column count exceeds the size_t range");
    return columns * 2;
}

class IntervalMatrix final {
public:
    IntervalMatrix(std::size_t rows, std::size_t columns, std::vector<ComplexInterval> values)
        : rows_(rows), columns_(columns), values_(std::move(values)) {
        const std::size_t shape[] = {rows_, columns_};
        if (values_.size() != expression::arrayElementCount(shape))
            throw std::invalid_argument("IntervalMatrix shape does not match the element count");
    }

    [[nodiscard]] std::size_t rows() const noexcept { return rows_; }
    [[nodiscard]] std::size_t columns() const noexcept { return columns_; }
    [[nodiscard]] ComplexInterval& operator()(std::size_t row, std::size_t column) noexcept {
        return values_[row * columns_ + column];
    }
    [[nodiscard]] const ComplexInterval& operator()(std::size_t row, std::size_t column) const noexcept {
        return values_[row * columns_ + column];
    }
    void swapRows(std::size_t lhs, std::size_t rhs) noexcept {
        if (lhs == rhs)
            return;
        for (std::size_t column = 0; column < columns_; ++column)
            std::swap((*this)(lhs, column), (*this)(rhs, column));
    }

private:
    std::size_t rows_ = 0;
    std::size_t columns_ = 0;
    std::vector<ComplexInterval> values_;
};

[[nodiscard]] std::optional<std::vector<ComplexInterval>> encloseElements(
    const ArrayExpr& array,
    std::size_t precisionBits,
    const approximation::CertifiedEvaluator& certified,
    approximation::CertifiedEvaluator::EnclosureKind enclosureKind =
        approximation::CertifiedEvaluator::EnclosureKind::Certified) {
    std::vector<ComplexInterval> values;
    values.reserve(array.size());
    for (std::size_t i = 0; i < array.size(); ++i) {
        const Expr element = array.element(i);
        const auto enclosed = approximation::encloseComplexExpression(
            element, precisionBits, certified, enclosureKind);
        if (!enclosed)
            return std::nullopt;
        values.push_back(*enclosed);
    }
    return values;
}

struct EnclosedElements final {
    std::vector<ComplexInterval> certified;
    std::vector<ComplexInterval> information;
};

[[nodiscard]] std::optional<EnclosedElements> encloseElementsWithInformation(
    const ArrayExpr& array,
    std::size_t precisionBits,
    const approximation::CertifiedEvaluator& certified) {
    auto information = encloseElements(array, precisionBits, certified,
        approximation::CertifiedEvaluator::EnclosureKind::Information);
    if (!information)
        return std::nullopt;
    auto value = encloseElements(array, precisionBits, certified,
        approximation::CertifiedEvaluator::EnclosureKind::Certified);
    if (!value)
        return std::nullopt;
    return EnclosedElements{std::move(*value), std::move(*information)};
}

struct DualIntervalMatrix final {
    IntervalMatrix certified;
    IntervalMatrix information;

    void swapRows(std::size_t lhs, std::size_t rhs) noexcept {
        certified.swapRows(lhs, rhs);
        information.swapRows(lhs, rhs);
    }
};

[[nodiscard]] std::optional<DualIntervalMatrix> encloseMatrix(
    const ArrayExpr& array,
    std::size_t precisionBits,
    const approximation::CertifiedEvaluator& certified) {
    auto values = encloseElementsWithInformation(array, precisionBits, certified);
    if (!values)
        return std::nullopt;
    return DualIntervalMatrix{
        IntervalMatrix{array.shape[0], array.shape[1], std::move(values->certified)},
        IntervalMatrix{array.shape[0], array.shape[1], std::move(values->information)}};
}

[[nodiscard]] std::optional<std::size_t> selectPivot(
    const IntervalMatrix& matrix,
    std::size_t firstRow,
    std::size_t column,
    approximation::PrecisionInsufficientKind kind =
        approximation::PrecisionInsufficientKind::Refinable) {
    bool uncertain = false;
    for (std::size_t row = firstRow; row < matrix.rows(); ++row) {
        const ComplexInterval& value = matrix(row, column);
        if (!value.containsZero())
            return row;
        if (!exactZero(value))
            uncertain = true;
    }
    if (uncertain)
        throw approximation::PrecisionInsufficient(
            "Matrix pivot nonzero status is not certified at the current precision", kind);
    return std::nullopt;
}

struct RrefResult final {
    IntervalMatrix certified;
    IntervalMatrix information;
    std::size_t rank = 0;
    std::vector<std::size_t> pivotColumns;
};

[[nodiscard]] RrefResult intervalRref(
    DualIntervalMatrix matrix,
    std::size_t precisionBits,
    std::size_t pivotColumnLimit) {
    std::size_t pivotRow = 0;
    const std::size_t limit = std::min(pivotColumnLimit, matrix.certified.columns());
    std::vector<std::size_t> pivotColumns;
    pivotColumns.reserve(std::min(matrix.certified.rows(), limit));

    for (std::size_t column = 0;
         column < limit && pivotRow < matrix.certified.rows(); ++column) {
        const auto selected = selectPivot(
            matrix.information, pivotRow, column,
            approximation::PrecisionInsufficientKind::InputInformation);
        if (!selected)
            continue;

        matrix.swapRows(*selected, pivotRow);
        const ComplexInterval certifiedPivot = matrix.certified(pivotRow, column);
        const ComplexInterval informationPivot = matrix.information(pivotRow, column);
        for (std::size_t c = 0; c < matrix.certified.columns(); ++c) {
            matrix.certified(pivotRow, c) = approximation::divide(
                matrix.certified(pivotRow, c), certifiedPivot, precisionBits);
            matrix.information(pivotRow, c) = approximation::divide(
                matrix.information(pivotRow, c), informationPivot, precisionBits);
        }

        for (std::size_t row = 0; row < matrix.certified.rows(); ++row) {
            if (row == pivotRow || exactZero(matrix.information(row, column)))
                continue;
            const ComplexInterval certifiedFactor = matrix.certified(row, column);
            const ComplexInterval informationFactor = matrix.information(row, column);
            for (std::size_t c = 0; c < matrix.certified.columns(); ++c) {
                matrix.certified(row, c) = approximation::subtract(
                    matrix.certified(row, c),
                    approximation::multiply(
                        certifiedFactor, matrix.certified(pivotRow, c), precisionBits),
                    precisionBits);
                matrix.information(row, c) = approximation::subtract(
                    matrix.information(row, c),
                    approximation::multiply(
                        informationFactor, matrix.information(pivotRow, c), precisionBits),
                    precisionBits);
            }
        }
        pivotColumns.push_back(column);
        ++pivotRow;
    }
    return RrefResult{
        std::move(matrix.certified), std::move(matrix.information),
        pivotRow, std::move(pivotColumns)};
}

[[nodiscard]] std::optional<Expr> decimalArray(
    std::vector<std::size_t> shape,
    const std::vector<ComplexInterval>& certified,
    const std::vector<ComplexInterval>& information,
    std::size_t digits) {
    if (certified.size() != information.size())
        throw std::invalid_argument("Certified and information arrays must have the same size");
    std::vector<Expr> output;
    output.reserve(certified.size());
    for (std::size_t i = 0; i < certified.size(); ++i) {
        const auto decimal = approximation::finalizeCertifiedApproximation(
            approximation::CertifiedValue{certified[i]},
            approximation::CertifiedValue{information[i]}, digits);
        if (!decimal)
            return std::nullopt;
        output.push_back(*decimal);
    }
    return Expr::array(std::move(shape), std::move(output));
}

template <class Operation>
[[nodiscard]] std::optional<Expr> retryApproximation(
    approximation::ApproximationContext context,
    Operation&& operation) {
    for (std::size_t attempt = 0; attempt < maximumPrecisionRetries; ++attempt) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        try {
            if (const auto result = operation(context))
                return result;
        }
        catch (const approximation::PrecisionInsufficient& exception) {
            if (!exception.refinable())
                return std::nullopt;
        }
        catch (const approximation::CertifiedBackendUnsupported&) {
            return std::nullopt;
        }
        context.setGuardDigits(approximation::nextGuardDigits(context.guardDigits()));
    }
    return std::nullopt;
}

[[nodiscard]] RealInterval normSquared(
    const std::vector<ComplexInterval>& values,
    std::size_t bits) {
    RealInterval result = RealInterval::fromRational(Rational{}, bits);
    for (const ComplexInterval& value : values) {
        const RealInterval magnitudeSquared = approximation::add(
            approximation::squareInterval(value.real(), bits),
            approximation::squareInterval(value.imaginary(), bits),
            bits);
        result = approximation::add(result, magnitudeSquared, bits);
    }
    return result;
}

} // namespace

std::optional<Expr> approximateDot(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 2 || !arguments[0].isArray() || !arguments[1].isArray())
        return std::nullopt;
    const ArrayExpr& lhs = arguments[0].asArray();
    const ArrayExpr& rhs = arguments[1].asArray();
    if (lhs.rank() < 1 || lhs.rank() > 2 || rhs.rank() < 1 || rhs.rank() > 2)
        return std::nullopt;

    const std::size_t lhsRows = lhs.rank() == 1 ? 1 : lhs.shape[0];
    const std::size_t lhsColumns = lhs.rank() == 1 ? lhs.shape[0] : lhs.shape[1];
    const std::size_t rhsRows = rhs.shape[0];
    const std::size_t rhsColumns = rhs.rank() == 1 ? 1 : rhs.shape[1];
    if (lhsColumns != rhsRows)
        return std::nullopt;

    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto left = encloseElementsWithInformation(lhs, bits, certified);
        const auto right = encloseElementsWithInformation(rhs, bits, certified);
        if (!left || !right)
            return std::nullopt;

        auto lhsAt = [&](const std::vector<ComplexInterval>& values,
                         std::size_t row, std::size_t column) -> const ComplexInterval& {
            return lhs.rank() == 1 ? values[column] : values[row * lhsColumns + column];
        };
        auto rhsAt = [&](const std::vector<ComplexInterval>& values,
                         std::size_t row, std::size_t column) -> const ComplexInterval& {
            return rhs.rank() == 1 ? values[row] : values[row * rhsColumns + column];
        };

        std::vector<ComplexInterval> output;
        std::vector<ComplexInterval> informationOutput;
        const std::size_t outputShape[] = {lhsRows, rhsColumns};
        output.reserve(expression::arrayElementCount(outputShape));
        informationOutput.reserve(expression::arrayElementCount(outputShape));
        for (std::size_t row = 0; row < lhsRows; ++row) {
            for (std::size_t column = 0; column < rhsColumns; ++column) {
                ComplexInterval sum = exactComplex(0, bits);
                ComplexInterval informationSum = exactComplex(0, bits);
                for (std::size_t k = 0; k < lhsColumns; ++k) {
                    sum = approximation::add(sum,
                        approximation::multiply(
                            lhsAt(left->certified, row, k),
                            rhsAt(right->certified, k, column), bits), bits);
                    informationSum = approximation::add(informationSum,
                        approximation::multiply(
                            lhsAt(left->information, row, k),
                            rhsAt(right->information, k, column), bits), bits);
                }
                output.push_back(std::move(sum));
                informationOutput.push_back(std::move(informationSum));
            }
        }

        if (lhs.rank() == 1 && rhs.rank() == 1)
            return approximation::finalizeCertifiedApproximation(
                approximation::CertifiedValue{output.front()},
                approximation::CertifiedValue{informationOutput.front()},
                current.decimalDigits());
        const std::vector<std::size_t> shape = lhs.rank() == 1 || rhs.rank() == 1
            ? std::vector<std::size_t>{lhs.rank() == 1 ? rhsColumns : lhsRows}
            : std::vector<std::size_t>{lhsRows, rhsColumns};
        return decimalArray(shape, output, informationOutput, current.decimalDigits());
    });
}

std::optional<Expr> approximateDeterminant(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& array = arguments[0].asArray();
    if (array.shape[0] != array.shape[1])
        return std::nullopt;

    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        auto matrix = encloseMatrix(array, bits, certified);
        if (!matrix)
            return std::nullopt;

        ComplexInterval determinant = exactComplex(1, bits);
        ComplexInterval informationDeterminant = exactComplex(1, bits);
        bool negative = false;
        for (std::size_t column = 0; column < matrix->certified.columns(); ++column) {
            const auto pivot = selectPivot(
                matrix->information, column, column,
                approximation::PrecisionInsufficientKind::InputInformation);
            if (!pivot) {
                const ComplexInterval zero = exactComplex(0, bits);
                return approximation::finalizeCertifiedApproximation(
                    approximation::CertifiedValue{zero},
                    approximation::CertifiedValue{zero}, current.decimalDigits());
            }
            if (*pivot != column) {
                matrix->swapRows(*pivot, column);
                negative = !negative;
            }

            const ComplexInterval pivotValue = matrix->certified(column, column);
            const ComplexInterval informationPivot = matrix->information(column, column);
            determinant = approximation::multiply(determinant, pivotValue, bits);
            informationDeterminant = approximation::multiply(
                informationDeterminant, informationPivot, bits);
            for (std::size_t row = column + 1; row < matrix->certified.rows(); ++row) {
                if (exactZero(matrix->information(row, column)))
                    continue;
                const ComplexInterval factor = approximation::divide(
                    matrix->certified(row, column), pivotValue, bits);
                const ComplexInterval informationFactor = approximation::divide(
                    matrix->information(row, column), informationPivot, bits);
                matrix->certified(row, column) = exactComplex(0, bits);
                matrix->information(row, column) = exactComplex(0, bits);
                for (std::size_t c = column + 1; c < matrix->certified.columns(); ++c) {
                    matrix->certified(row, c) = approximation::subtract(
                        matrix->certified(row, c),
                        approximation::multiply(
                            factor, matrix->certified(column, c), bits), bits);
                    matrix->information(row, c) = approximation::subtract(
                        matrix->information(row, c),
                        approximation::multiply(
                            informationFactor, matrix->information(column, c), bits), bits);
                }
            }
        }
        if (negative) {
            determinant = approximation::negate(determinant);
            informationDeterminant = approximation::negate(informationDeterminant);
        }
        return approximation::finalizeCertifiedApproximation(
            approximation::CertifiedValue{determinant},
            approximation::CertifiedValue{informationDeterminant}, current.decimalDigits());
    });
}

std::optional<Expr> approximateInverse(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    if (source.shape[0] != source.shape[1])
        return std::nullopt;

    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto input = encloseElementsWithInformation(source, bits, certified);
        if (!input)
            return std::nullopt;
        const std::size_t n = source.shape[0];
        const std::size_t columns = augmentedColumns(n);
        const std::size_t augmentedShape[] = {n, columns};

        std::vector<ComplexInterval> values;
        std::vector<ComplexInterval> informationValues;
        values.reserve(expression::arrayElementCount(augmentedShape));
        informationValues.reserve(expression::arrayElementCount(augmentedShape));
        for (std::size_t row = 0; row < n; ++row) {
            for (std::size_t column = 0; column < n; ++column) {
                values.push_back(input->certified[row * n + column]);
                informationValues.push_back(input->information[row * n + column]);
            }
            for (std::size_t column = 0; column < n; ++column) {
                values.push_back(exactComplex(row == column ? 1 : 0, bits));
                informationValues.push_back(exactComplex(row == column ? 1 : 0, bits));
            }
        }
        DualIntervalMatrix augmented{
            IntervalMatrix{n, columns, std::move(values)},
            IntervalMatrix{n, columns, std::move(informationValues)}};
        const RrefResult reduced = intervalRref(std::move(augmented), bits, n);
        if (reduced.rank != n)
            throw std::domain_error("Matrix is singular");

        std::vector<ComplexInterval> output;
        std::vector<ComplexInterval> informationOutput;
        output.reserve(source.size());
        informationOutput.reserve(source.size());
        for (std::size_t row = 0; row < n; ++row)
            for (std::size_t column = 0; column < n; ++column) {
                output.push_back(reduced.certified(row, n + column));
                informationOutput.push_back(reduced.information(row, n + column));
            }
        return decimalArray({n, n}, output, informationOutput, current.decimalDigits());
    });
}

std::optional<Expr> approximateRref(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        auto matrix = encloseMatrix(source, bits, certified);
        if (!matrix)
            return std::nullopt;
        const RrefResult reduced = intervalRref(std::move(*matrix), bits, source.shape[1]);

        std::vector<ComplexInterval> output;
        std::vector<ComplexInterval> informationOutput;
        output.reserve(source.size());
        informationOutput.reserve(source.size());
        for (std::size_t row = 0; row < reduced.certified.rows(); ++row)
            for (std::size_t column = 0; column < reduced.certified.columns(); ++column) {
                output.push_back(reduced.certified(row, column));
                informationOutput.push_back(reduced.information(row, column));
            }
        return decimalArray(
            source.shape, output, informationOutput, current.decimalDigits());
    });
}

std::optional<Expr> approximateMatrixRank(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        auto matrix = encloseMatrix(source, bits, certified);
        if (!matrix)
            return std::nullopt;
        const RrefResult reduced = intervalRref(std::move(*matrix), bits, source.shape[1]);
        return Expr{numeric::Number{BigInt::parse(std::to_string(reduced.rank))}};
    });
}

std::optional<Expr> approximateSolveLinear(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 2 || !arguments[0].isArray() || !arguments[1].isArray()
        || !arguments[0].asArray().isMatrix() || !arguments[1].asArray().isVector())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    const ArrayExpr& rhs = arguments[1].asArray();
    if (rhs.shape[0] != source.shape[0])
        return std::nullopt;
    if (source.shape[1] == std::numeric_limits<std::size_t>::max())
        throw std::length_error("Linear system augmented column count exceeds the size_t range");

    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};
    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto coefficients = encloseElementsWithInformation(source, bits, certified);
        const auto right = encloseElementsWithInformation(rhs, bits, certified);
        if (!coefficients || !right)
            return std::nullopt;

        const std::size_t rows = source.shape[0];
        const std::size_t variables = source.shape[1];
        const std::size_t augmentedColumns = variables + 1;
        const std::size_t augmentedShape[] = {rows, augmentedColumns};
        std::vector<ComplexInterval> values;
        std::vector<ComplexInterval> informationValues;
        values.reserve(expression::arrayElementCount(augmentedShape));
        informationValues.reserve(expression::arrayElementCount(augmentedShape));
        for (std::size_t row = 0; row < rows; ++row) {
            for (std::size_t column = 0; column < variables; ++column) {
                values.push_back(coefficients->certified[row * variables + column]);
                informationValues.push_back(
                    coefficients->information[row * variables + column]);
            }
            values.push_back(right->certified[row]);
            informationValues.push_back(right->information[row]);
        }

        RrefResult reduced = intervalRref(
            DualIntervalMatrix{
                IntervalMatrix{rows, augmentedColumns, std::move(values)},
                IntervalMatrix{rows, augmentedColumns, std::move(informationValues)}},
            bits, variables);
        for (std::size_t row = reduced.rank; row < rows; ++row) {
            const ComplexInterval& residual = reduced.information(row, variables);
            if (!residual.containsZero())
                throw std::domain_error("Linear system is inconsistent");
            if (!exactZero(residual))
                throw approximation::PrecisionInsufficient(
                    "Linear system consistency is not certified by the input information",
                    approximation::PrecisionInsufficientKind::InputInformation);
        }
        if (reduced.rank != variables)
            throw std::domain_error("Linear system does not have a unique solution");

        std::vector<ComplexInterval> solution;
        std::vector<ComplexInterval> informationSolution;
        solution.reserve(variables);
        informationSolution.reserve(variables);
        for (std::size_t variable = 0; variable < variables; ++variable) {
            solution.push_back(reduced.certified(variable, variables));
            informationSolution.push_back(reduced.information(variable, variables));
        }
        return decimalArray(
            {variables}, solution, informationSolution, current.decimalDigits());
    });
}

std::optional<Expr> approximateNullSpace(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& source = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        auto matrix = encloseMatrix(source, bits, certified);
        if (!matrix)
            return std::nullopt;
        RrefResult reduced = intervalRref(std::move(*matrix), bits, source.shape[1]);

        const std::size_t variables = source.shape[1];
        const std::size_t noPivot = std::numeric_limits<std::size_t>::max();
        std::vector<std::size_t> pivotRowByColumn(variables, noPivot);
        for (std::size_t row = 0; row < reduced.pivotColumns.size(); ++row)
            pivotRowByColumn[reduced.pivotColumns[row]] = row;

        const std::size_t nullity = variables - reduced.pivotColumns.size();
        std::vector<ComplexInterval> basis;
        const std::size_t shape[] = {nullity, variables};
        basis.reserve(expression::arrayElementCount(shape));
        for (std::size_t freeColumn = 0; freeColumn < variables; ++freeColumn) {
            if (pivotRowByColumn[freeColumn] != noPivot)
                continue;
            for (std::size_t variable = 0; variable < variables; ++variable) {
                if (variable == freeColumn) {
                    basis.push_back(exactComplex(1, bits));
                    continue;
                }

                const std::size_t pivotRow = pivotRowByColumn[variable];
                if (pivotRow == noPivot) {
                    basis.push_back(exactComplex(0, bits));
                    continue;
                }
                basis.push_back(approximation::negate(
                    reduced.certified(pivotRow, freeColumn)));
            }
        }
        std::vector<ComplexInterval> informationBasis;
        informationBasis.reserve(basis.size());
        for (std::size_t freeColumn = 0; freeColumn < variables; ++freeColumn) {
            if (pivotRowByColumn[freeColumn] != noPivot)
                continue;
            for (std::size_t variable = 0; variable < variables; ++variable) {
                if (variable == freeColumn) {
                    informationBasis.push_back(exactComplex(1, bits));
                    continue;
                }
                const std::size_t pivotRow = pivotRowByColumn[variable];
                if (pivotRow == noPivot) {
                    informationBasis.push_back(exactComplex(0, bits));
                    continue;
                }
                informationBasis.push_back(approximation::negate(
                    reduced.information(pivotRow, freeColumn)));
            }
        }
        return decimalArray(
            {nullity, variables}, basis, informationBasis, current.decimalDigits());
    });
}

std::optional<Expr> approximateNorm(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isVector())
        return std::nullopt;
    const ArrayExpr& vector = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto values = encloseElementsWithInformation(vector, bits, certified);
        if (!values)
            return std::nullopt;
        const RealInterval squared = normSquared(values->certified, bits);
        const RealInterval informationSquared = normSquared(values->information, bits);
        const RealInterval norm = approximation::encloseSqrt(squared, bits).interval;
        const RealInterval informationNorm = approximation::encloseSqrt(
            informationSquared, bits).interval;
        return approximation::finalizeCertifiedApproximation(
            approximation::CertifiedValue{norm},
            approximation::CertifiedValue{informationNorm}, current.decimalDigits());
    });
}

std::optional<Expr> approximateNormalize(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isVector())
        return std::nullopt;
    const ArrayExpr& vector = arguments[0].asArray();
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto values = encloseElementsWithInformation(vector, bits, certified);
        if (!values)
            return std::nullopt;
        const RealInterval squared = normSquared(values->certified, bits);
        const RealInterval informationSquared = normSquared(values->information, bits);
        const RealInterval norm = approximation::encloseSqrt(squared, bits).interval;
        const RealInterval informationNorm = approximation::encloseSqrt(
            informationSquared, bits).interval;
        if (exactZero(informationNorm))
            throw std::domain_error("normalize requires a nonzero vector");
        if (informationNorm.containsZero())
            throw approximation::PrecisionInsufficient(
                "Vector norm nonzero status is not certified by the input information",
                approximation::PrecisionInsufficientKind::InputInformation);

        const ComplexInterval denominator = ComplexInterval::fromReal(norm);
        const ComplexInterval informationDenominator = ComplexInterval::fromReal(informationNorm);
        std::vector<ComplexInterval> output;
        std::vector<ComplexInterval> informationOutput;
        output.reserve(values->certified.size());
        informationOutput.reserve(values->information.size());
        for (const ComplexInterval& value : values->certified)
            output.push_back(approximation::divide(value, denominator, bits));
        for (const ComplexInterval& value : values->information)
            informationOutput.push_back(
                approximation::divide(value, informationDenominator, bits));
        return decimalArray(
            {vector.shape[0]}, output, informationOutput, current.decimalDigits());
    });
}

std::optional<Expr> approximateTrace(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    approximation::ApproximationContext context) {
    if (arguments.size() != 1 || !arguments[0].isArray()
        || !arguments[0].asArray().isMatrix())
        return std::nullopt;
    const ArrayExpr& matrix = arguments[0].asArray();
    if (matrix.shape[0] != matrix.shape[1])
        return std::nullopt;
    approximation::CertifiedEvaluator certified{builtins, mathematics, angles};

    return retryApproximation(context, [&](const auto& current) -> std::optional<Expr> {
        const std::size_t bits = current.workingBinaryBits();
        const auto values = encloseElementsWithInformation(matrix, bits, certified);
        if (!values)
            return std::nullopt;
        ComplexInterval sum = exactComplex(0, bits);
        ComplexInterval informationSum = exactComplex(0, bits);
        for (std::size_t i = 0; i < matrix.shape[0]; ++i) {
            const std::size_t index = i * matrix.shape[1] + i;
            sum = approximation::add(sum, values->certified[index], bits);
            informationSum = approximation::add(
                informationSum, values->information[index], bits);
        }
        return approximation::finalizeCertifiedApproximation(
            approximation::CertifiedValue{sum},
            approximation::CertifiedValue{informationSum}, current.decimalDigits());
    });
}

} // namespace mmcal::linear_algebra
