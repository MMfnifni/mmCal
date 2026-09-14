#include "algebraic_field_linear.hpp"

#include "evaluation/evaluation_budget.hpp"
#include "numeric/big_int.hpp"

#include <algorithm>
#include <memory>
#include <utility>

namespace mmcal::linear_algebra {
namespace {

using numeric::BigInt;
using numeric::Rational;
using symbolic::AlgebraicElement;
using symbolic::NumberFieldContext;

[[nodiscard]] std::optional<AlgebraicElement> fieldConstant(
    const std::shared_ptr<const NumberFieldContext>& field,
    const Rational& value) {
    if (!field || field->degree() == 0)
        return std::nullopt;
    std::vector<Rational> coefficients(field->degree());
    coefficients[0] = value;
    return AlgebraicElement::create(field, std::move(coefficients));
}

[[nodiscard]] bool sameField(
    std::span<const AlgebraicElement> elements,
    const std::shared_ptr<const NumberFieldContext>& field) noexcept {
    return std::all_of(elements.begin(), elements.end(), [&](const AlgebraicElement& value) {
        return value.field().get() == field.get();
    });
}

} // namespace

std::optional<AlgebraicFieldRref> algebraicFieldRref(
    std::size_t rows,
    std::size_t columns,
    std::span<const AlgebraicElement> input) {
    if (rows == 0 || columns == 0 || input.size() != rows * columns)
        return std::nullopt;
    const auto field = input.front().field();
    if (!field || !sameField(input, field))
        return std::nullopt;

    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::TemporaryMatrixElement,
        rows * columns);

    AlgebraicFieldRref result;
    result.rows = rows;
    result.columns = columns;
    result.elements.assign(input.begin(), input.end());
    result.pivotColumns.reserve(std::min(rows, columns));

    std::size_t pivotRow = 0;
    for (std::size_t column = 0; column < columns && pivotRow < rows; ++column) {
        std::size_t selected = pivotRow;
        while (selected < rows
            && result.elements[selected * columns + column].isZero())
            ++selected;
        if (selected == rows)
            continue;

        if (selected != pivotRow) {
            for (std::size_t j = 0; j < columns; ++j)
                std::swap(
                    result.elements[pivotRow * columns + j],
                    result.elements[selected * columns + j]);
        }

        const AlgebraicElement pivot = result.elements[pivotRow * columns + column];
        for (std::size_t j = column; j < columns; ++j) {
            auto normalized = result.elements[pivotRow * columns + j].divide(pivot);
            if (!normalized)
                return std::nullopt;
            result.elements[pivotRow * columns + j] = std::move(*normalized);
        }

        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow)
                continue;
            const AlgebraicElement factor = result.elements[row * columns + column];
            if (factor.isZero())
                continue;
            for (std::size_t j = column; j < columns; ++j) {
                auto product = factor.multiply(result.elements[pivotRow * columns + j]);
                if (!product)
                    return std::nullopt;
                auto reduced = result.elements[row * columns + j].subtract(*product);
                if (!reduced)
                    return std::nullopt;
                result.elements[row * columns + j] = std::move(*reduced);
            }
        }

        result.pivotColumns.push_back(column);
        ++pivotRow;
    }
    return result;
}

std::optional<std::vector<std::vector<AlgebraicElement>>> algebraicFieldNullSpace(
    std::size_t rows,
    std::size_t columns,
    std::span<const AlgebraicElement> elements) {
    const auto reduced = algebraicFieldRref(rows, columns, elements);
    if (!reduced || elements.empty())
        return std::nullopt;

    const auto field = elements.front().field();
    const auto zero = fieldConstant(field, Rational{BigInt{0}});
    const auto one = fieldConstant(field, Rational{BigInt{1}});
    if (!zero || !one)
        return std::nullopt;

    std::vector<bool> pivot(columns, false);
    for (const std::size_t column : reduced->pivotColumns) {
        if (column >= columns)
            return std::nullopt;
        pivot[column] = true;
    }

    std::vector<std::vector<AlgebraicElement>> basis;
    basis.reserve(columns - reduced->pivotColumns.size());
    for (std::size_t freeColumn = 0; freeColumn < columns; ++freeColumn) {
        if (pivot[freeColumn])
            continue;

        std::vector<AlgebraicElement> vector(columns, *zero);
        vector[freeColumn] = *one;
        for (std::size_t row = 0; row < reduced->pivotColumns.size(); ++row) {
            const std::size_t pivotColumn = reduced->pivotColumns[row];
            auto negative = zero->subtract((*reduced)(row, freeColumn));
            if (!negative)
                return std::nullopt;
            vector[pivotColumn] = std::move(*negative);
        }
        basis.push_back(std::move(vector));
    }
    return basis;
}

} // namespace mmcal::linear_algebra
