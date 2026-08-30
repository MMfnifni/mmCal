// 配列・ベクトル操作
#include "array_vector.hpp"

#include "builtins/array_helpers.hpp"
#include "builtins/exact_operations.hpp"
#include "error/error_message.hpp"
#include "evaluation/evaluation_budget.hpp"
#include "mathematics/value_facts.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <algorithm>
#include <charconv>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::ArrayExpr;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;

[[nodiscard]] Expr integer(std::int64_t value) { return Expr{Number{BigInt{value}}}; }

void requireSameLength(const ArrayExpr& a, const ArrayExpr& b, std::string_view name) {
    if (a.shape[0] != b.shape[0])
        error::throwCalcError(error::CalcErrorType::Domain,
            std::string{name} + " requires vectors with the same length");
}

[[nodiscard]] bool provablyRealVector(
    const ArrayExpr& vector,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics) {
    if (vector.hasExactRealStorage())
        return true;
    for (std::size_t i = 0; i < vector.size(); ++i) {
        const Expr item = vector.element(i);
        if (!mathematics::inferValueFacts(item, registry, mathematics).isProvablyReal())
            return false;
    }
    return true;
}

[[nodiscard]] Expr vectorExpr(std::vector<Expr> elements) {
    const std::size_t size = elements.size();
    return Expr::array({size}, std::move(elements));
}

[[nodiscard]] Expr dot(
    const ArrayExpr& a, const ArrayExpr& b,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireSameLength(a, b, "vdot");
    if (a.hasExactNumberStorage() && b.hasExactNumberStorage()) {
        Number result{BigInt{0}};
        for (std::size_t i = 0; i < a.size(); ++i)
            result += a.exactNumber(i) * b.exactNumber(i);
        return Expr{std::move(result)};
    }
    std::vector<Expr> terms;
    terms.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        terms.push_back(exact::multiply({a.element(i), b.element(i)}, registry, mathematics, angles));
    return exact::add(std::move(terms), registry, mathematics, angles);
}

[[nodiscard]] Expr realNorm(
    const ArrayExpr& vector,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!provablyRealVector(vector, registry, mathematics))
        return Expr::call(registry.symbol(BuiltinId::VectorNorm), {vectorExpr(vector.materialize())});
    std::vector<Expr> squares;
    squares.reserve(vector.size());
    for (std::size_t i = 0; i < vector.size(); ++i) {
        const Expr item = vector.element(i);
        squares.push_back(exact::multiply({item, item}, registry, mathematics, angles));
    }
    return exact::sqrt(exact::add(std::move(squares), registry, mathematics, angles),
        registry, mathematics, angles);
}

[[nodiscard]] Expr vectorDifference(
    const ArrayExpr& a, const ArrayExpr& b,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireSameLength(a, b, "vector distance");
    if (a.hasExactNumberStorage() && b.hasExactNumberStorage()) {
        std::vector<Number> values;
        values.reserve(a.size());
        for (std::size_t i = 0; i < a.size(); ++i)
            values.push_back(a.exactNumber(i) - b.exactNumber(i));
        return Expr::numberArray({a.size()}, std::move(values));
    }
    std::vector<Expr> result;
    result.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        result.push_back(exact::subtract(a.element(i), b.element(i), registry, mathematics, angles));
    return vectorExpr(std::move(result));
}

[[nodiscard]] Expr scaleVector(
    const ArrayExpr& vector, Expr scale,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (vector.hasExactNumberStorage() && scale.isNumber()) {
        std::vector<Number> values;
        values.reserve(vector.size());
        for (std::size_t i = 0; i < vector.size(); ++i)
            values.push_back(vector.exactNumber(i) * scale.asNumber());
        return Expr::numberArray({vector.size()}, std::move(values));
    }
    std::vector<Expr> result;
    result.reserve(vector.size());
    for (std::size_t i = 0; i < vector.size(); ++i) {
        const Expr item = vector.element(i);
        result.push_back(exact::multiply({item, scale}, registry, mathematics, angles));
    }
    return vectorExpr(std::move(result));
}

[[nodiscard]] Expr projection(
    const ArrayExpr& a, const ArrayExpr& b,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireSameLength(a, b, "vproject");
    if (!provablyRealVector(a, registry, mathematics)
        || !provablyRealVector(b, registry, mathematics))
        return Expr::call(registry.symbol(BuiltinId::VectorProject),
            {vectorExpr(a.materialize()), vectorExpr(b.materialize())});
    Expr denominator = dot(b, b, registry, mathematics, angles);
    if (denominator.isNumber() && denominator.asNumber().isZero())
        error::throwCalcError(error::CalcErrorType::Domain,
            "vproject requires a nonzero direction vector");
    Expr factor = exact::divide(dot(a, b, registry, mathematics, angles), denominator,
        registry, mathematics, angles);
    return scaleVector(b, std::move(factor), registry, mathematics, angles);
}

} // namespace

Expr evaluateIdentity(std::span<const Expr> arguments) {
    const std::size_t n = detail::requireSize(arguments[0], "identity");
    const std::size_t shape[] = {n, n};
    const std::size_t count = expression::arrayElementCount(shape);
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::DenseArrayElement, count);
    std::vector<BigInt> elements;
    elements.reserve(count);
    for (std::size_t r = 0; r < n; ++r)
        for (std::size_t c = 0; c < n; ++c)
            elements.emplace_back(r == c ? 1 : 0);
    return Expr::integerArray({n, n}, std::move(elements));
}

Expr evaluateZeros(std::span<const Expr> arguments) {
    const std::size_t rows = detail::requireSize(arguments[0], "zeros");
    const std::size_t cols = detail::requireSize(arguments[1], "zeros");
    const std::size_t shape[] = {rows, cols};
    const std::size_t count = expression::arrayElementCount(shape);
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::DenseArrayElement, count);
    return Expr::integerArray({rows, cols}, std::vector<BigInt>(count, BigInt{0}));
}

Expr evaluateRows(std::span<const Expr> arguments) {
    const ArrayExpr& matrix = detail::requireMatrix(arguments[0], "rows");
    return detail::sizeExpr(matrix.shape[0]);
}

Expr evaluateCols(std::span<const Expr> arguments) {
    const ArrayExpr& matrix = detail::requireMatrix(arguments[0], "cols");
    return detail::sizeExpr(matrix.shape[1]);
}

Expr evaluateDiag(std::span<const Expr> arguments) {
    const ArrayExpr& matrix = detail::requireMatrix(arguments[0], "diag");
    const std::size_t count = std::min(matrix.shape[0], matrix.shape[1]);
    std::vector<Expr> diagonal;
    diagonal.reserve(count);
    for (std::size_t i = 0; i < count; ++i)
        diagonal.push_back(matrix.element(i * matrix.shape[1] + i));
    return vectorExpr(std::move(diagonal));
}

Expr evaluateVectorAdd(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& a = detail::requireVector(arguments[0], "vadd");
    const ArrayExpr& b = detail::requireVector(arguments[1], "vadd");
    requireSameLength(a, b, "vadd");
    if (a.hasExactNumberStorage() && b.hasExactNumberStorage()) {
        std::vector<Number> values;
        values.reserve(a.size());
        for (std::size_t i = 0; i < a.size(); ++i)
            values.push_back(a.exactNumber(i) + b.exactNumber(i));
        return Expr::numberArray({a.size()}, std::move(values));
    }
    std::vector<Expr> result;
    result.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        result.push_back(exact::add({a.element(i), b.element(i)}, registry, mathematics, angles));
    return vectorExpr(std::move(result));
}

Expr evaluateVectorSubtract(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return vectorDifference(detail::requireVector(arguments[0], "vsub"), detail::requireVector(arguments[1], "vsub"),
        registry, mathematics, angles);
}

Expr evaluateVectorScale(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (arguments[1].isArray() || arguments[1].isString() || arguments[1].isBoolean())
        error::throwCalcError(error::CalcErrorType::Type, "vscalar scale must be a scalar expression");
    return scaleVector(detail::requireVector(arguments[0], "vscalar"), arguments[1], registry, mathematics, angles);
}

Expr evaluateVectorCross(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& a = detail::requireVector(arguments[0], "vcross");
    const ArrayExpr& b = detail::requireVector(arguments[1], "vcross");
    if (a.shape[0] != 3 || b.shape[0] != 3)
        error::throwCalcError(error::CalcErrorType::Domain, "vcross requires two 3D vectors");
    return vectorExpr({
        exact::subtract(exact::multiply({a.element(1), b.element(2)}, registry, mathematics, angles),
            exact::multiply({a.element(2), b.element(1)}, registry, mathematics, angles), registry, mathematics, angles),
        exact::subtract(exact::multiply({a.element(2), b.element(0)}, registry, mathematics, angles),
            exact::multiply({a.element(0), b.element(2)}, registry, mathematics, angles), registry, mathematics, angles),
        exact::subtract(exact::multiply({a.element(0), b.element(1)}, registry, mathematics, angles),
            exact::multiply({a.element(1), b.element(0)}, registry, mathematics, angles), registry, mathematics, angles)
    });
}

Expr evaluateVectorNorm(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return realNorm(detail::requireVector(arguments[0], "vnorm"), registry, mathematics, angles);
}

Expr evaluateVectorManhattan(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& a = detail::requireVector(arguments[0], "vmanhattan");
    const ArrayExpr& b = detail::requireVector(arguments[1], "vmanhattan");
    requireSameLength(a, b, "vmanhattan");
    if (!provablyRealVector(a, registry, mathematics) || !provablyRealVector(b, registry, mathematics))
        return Expr::call(registry.symbol(BuiltinId::VectorManhattan), {arguments[0], arguments[1]});
    std::vector<Expr> terms;
    terms.reserve(a.size());
    for (std::size_t i = 0; i < a.size(); ++i)
        terms.push_back(exact::call(BuiltinId::Abs,
            {exact::subtract(a.element(i), b.element(i), registry, mathematics, angles)},
            registry, mathematics, angles));
    return exact::add(std::move(terms), registry, mathematics, angles);
}

Expr evaluateVectorEuclidean(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr difference = vectorDifference(detail::requireVector(arguments[0], "veuclidean"),
        detail::requireVector(arguments[1], "veuclidean"), registry, mathematics, angles);
    return realNorm(difference.asArray(), registry, mathematics, angles);
}

Expr evaluateVectorNormalize(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& vector = detail::requireVector(arguments[0], "vnormalize");
    if (!provablyRealVector(vector, registry, mathematics))
        return Expr::call(registry.symbol(BuiltinId::VectorNormalize), {arguments[0]});
    Expr norm = realNorm(vector, registry, mathematics, angles);
    if (norm.isNumber() && norm.asNumber().isZero())
        error::throwCalcError(error::CalcErrorType::Domain, "vnormalize requires a nonzero vector");
    return scaleVector(vector, exact::divide(integer(1), norm, registry, mathematics, angles),
        registry, mathematics, angles);
}

Expr evaluateVectorProject(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return projection(detail::requireVector(arguments[0], "vproject"), detail::requireVector(arguments[1], "vproject"),
        registry, mathematics, angles);
}

Expr evaluateVectorAngle(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& a = detail::requireVector(arguments[0], "vangle");
    const ArrayExpr& b = detail::requireVector(arguments[1], "vangle");
    requireSameLength(a, b, "vangle");
    if (!provablyRealVector(a, registry, mathematics) || !provablyRealVector(b, registry, mathematics))
        return Expr::call(registry.symbol(BuiltinId::VectorAngle), {arguments[0], arguments[1]});
    Expr na = realNorm(a, registry, mathematics, angles);
    Expr nb = realNorm(b, registry, mathematics, angles);
    if ((na.isNumber() && na.asNumber().isZero()) || (nb.isNumber() && nb.asNumber().isZero()))
        error::throwCalcError(error::CalcErrorType::Domain, "vangle is undefined for the zero vector");
    Expr cosine = exact::divide(dot(a, b, registry, mathematics, angles),
        exact::multiply({na, nb}, registry, mathematics, angles), registry, mathematics, angles);
    return exact::call(BuiltinId::Acos, {std::move(cosine)}, registry, mathematics, angles);
}

Expr evaluateVectorReflect(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& a = detail::requireVector(arguments[0], "vreflect");
    const ArrayExpr& n = detail::requireVector(arguments[1], "vreflect");
    Expr p = projection(a, n, registry, mathematics, angles);
    if (!p.isArray())
        return Expr::call(registry.symbol(BuiltinId::VectorReflect), {arguments[0], arguments[1]});
    Expr twice = scaleVector(p.asArray(), integer(2), registry, mathematics, angles);
    return vectorDifference(a, twice.asArray(), registry, mathematics, angles);
}

Expr evaluateVectorReflectAxis(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& a = detail::requireVector(arguments[0], "vreflect_axis");
    const ArrayExpr& axis = detail::requireVector(arguments[1], "vreflect_axis");
    Expr p = projection(a, axis, registry, mathematics, angles);
    if (!p.isArray())
        return Expr::call(registry.symbol(BuiltinId::VectorReflectAxis), {arguments[0], arguments[1]});
    Expr twice = scaleVector(p.asArray(), integer(2), registry, mathematics, angles);
    return vectorDifference(twice.asArray(), a, registry, mathematics, angles);
}

Expr evaluateVectorSum(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const ArrayExpr& vector = detail::requireVector(arguments[0], "vsum");
    return exact::add(vector.materialize(), registry, mathematics, angles);
}

} // namespace mmcal::builtins
