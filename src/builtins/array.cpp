// Array shape・index・reshape
#include "array.hpp"

#include "builtins/array_helpers.hpp"
#include "error/error_message.hpp"
#include "expression/array_utils.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using expression::ArrayExpr;
using expression::Expr;

[[nodiscard]] Expr shapeExpr(const std::vector<std::size_t>& shape) {
    std::vector<Expr> result;
    result.reserve(shape.size());
    for (const std::size_t dimension : shape)
        result.push_back(detail::sizeExpr(dimension));
    const std::size_t rank = result.size();
    return Expr::array({rank}, std::move(result));
}

[[nodiscard]] std::vector<std::size_t> requireShape(
    const Expr& expression,
    std::string_view name) {
    const ArrayExpr& shape = detail::requireVector(expression, name);
    if (shape.elements.empty())
        detail::arrayTypeError(std::string{name} + " shape must contain at least one dimension");

    std::vector<std::size_t> dimensions;
    dimensions.reserve(shape.elements.size());
    for (const Expr& item : shape.elements)
        dimensions.push_back(detail::requireSize(item, name));
    return dimensions;
}

} // namespace

Expr evaluateDimensions(std::span<const Expr> arguments) {
    const Expr& value = arguments.front();
    if (!value.isArray() && !value.isList())
        detail::arrayTypeError("dimensions expects a brace value");
    return shapeExpr(expression::commonBraceDimensions(value));
}

Expr evaluateArrayRank(std::span<const Expr> arguments) {
    const Expr& value = arguments.front();
    if (!value.isArray() && !value.isList())
        detail::arrayTypeError("arrayRank expects a brace value");
    return detail::sizeExpr(expression::commonBraceDimensions(value).size());
}

Expr evaluateLength(std::span<const Expr> arguments) {
    const Expr& value = arguments.front();
    if (value.isArray())
        return detail::sizeExpr(value.asArray().shape.front());
    if (value.isList())
        return detail::sizeExpr(value.asList().elements.size());
    detail::arrayTypeError("length expects a brace value");
}

Expr evaluateArrayGet(std::span<const Expr> arguments) {
    if (arguments.front().isList()) {
        const auto& list = arguments.front().asList();
        const std::size_t index = detail::requireSize(arguments[1], "at");
        if (index >= list.elements.size())
            error::throwCalcError(error::CalcErrorType::Domain, "at index is out of range");
        if (arguments.size() == 2)
            return list.elements[index];

        std::vector<Expr> nested;
        nested.reserve(arguments.size() - 1);
        nested.push_back(list.elements[index]);
        for (std::size_t i = 2; i < arguments.size(); ++i)
            nested.push_back(arguments[i]);
        return evaluateArrayGet(nested);
    }

    const ArrayExpr& array = detail::requireArray(arguments.front(), "at");
    if (arguments.size() < 2 || arguments.size() > array.rank() + 1)
        detail::arrayTypeError(
            "at expects between one and arrayRank[A] zero-based indices");

    std::vector<std::size_t> indices;
    indices.reserve(arguments.size() - 1);
    for (std::size_t i = 1; i < arguments.size(); ++i)
        indices.push_back(detail::requireSize(arguments[i], "at"));

    for (std::size_t dimension = 0; dimension < indices.size(); ++dimension)
        if (indices[dimension] >= array.shape[dimension])
            error::throwCalcError(error::CalcErrorType::Domain, "at index is out of range");

    if (indices.size() == array.rank())
        return array.elements[array.flatIndex(indices)];

    // row-major Arrayではprefix indexで選ばれるsubarrayは常に連続領域。
    std::size_t prefix = 0;
    for (std::size_t dimension = 0; dimension < indices.size(); ++dimension)
        prefix = prefix * array.shape[dimension] + indices[dimension];

    std::vector<std::size_t> shape(
        array.shape.begin() + static_cast<std::ptrdiff_t>(indices.size()), array.shape.end());
    const std::size_t sliceSize = expression::arrayElementCount(shape);
    const std::size_t offset = prefix * sliceSize;
    std::vector<Expr> elements(
        array.elements.begin() + static_cast<std::ptrdiff_t>(offset),
        array.elements.begin() + static_cast<std::ptrdiff_t>(offset + sliceSize));
    return Expr::array(std::move(shape), std::move(elements));
}

Expr evaluateReshape(std::span<const Expr> arguments) {
    const ArrayExpr& array = detail::requireArray(arguments[0], "reshape");
    std::vector<std::size_t> shape = requireShape(arguments[1], "reshape");

    std::size_t count = 0;
    try {
        count = expression::arrayElementCount(shape);
    }
    catch (const std::length_error&) {
        error::throwCalcError(error::CalcErrorType::Overflow,
            "reshape dimensions overflow the addressable element count");
    }

    if (count != array.elements.size())
        error::throwCalcError(error::CalcErrorType::Domain,
            "reshape requires the same total element count");

    return Expr::array(std::move(shape), array.elements);
}

} // namespace mmcal::builtins
