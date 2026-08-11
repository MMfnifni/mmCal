// DFT・FFT・畳み込み
#include "signal_processing.hpp"

#include "builtins/exact_operations.hpp"
#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

void requireArity(std::span<const Expr> arguments, std::size_t expected, std::string_view name) {
    if (arguments.size() != expected)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " expects " + std::to_string(expected) + " argument(s)");
}

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] BigInt sizeInteger(std::size_t value) {
    return BigInt::parse(std::to_string(value));
}

[[nodiscard]] std::vector<Expr> vectorArgument(const Expr& expression, std::string_view name) {
    if (!expression.isArray() || expression.asArray().rank() != 1)
        error::throwCalcError(
            error::CalcErrorType::Type,
            std::string{name} + " requires a rank-1 array");
    return expression.asArray().elements;
}

[[nodiscard]] Expr zero() {
    return integer(0);
}

[[nodiscard]] Expr multiply(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() * rhs.asNumber()};
    return exact::multiply({std::move(lhs), std::move(rhs)}, registry, mathematics, angles);
}

[[nodiscard]] Expr add(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() + rhs.asNumber()};
    return exact::add({std::move(lhs), std::move(rhs)}, registry, mathematics, angles);
}

[[nodiscard]] Expr subtract(
    Expr lhs,
    Expr rhs,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (lhs.isNumber() && rhs.isNumber())
        return Expr{lhs.asNumber() - rhs.asNumber()};
    return exact::subtract(std::move(lhs), std::move(rhs), registry, mathematics, angles);
}

[[nodiscard]] Expr divideBySize(
    Expr value,
    std::size_t denominator,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (denominator == 0)
        return value;
    return exact::divide(
        std::move(value), Expr{Number{sizeInteger(denominator)}}, registry, mathematics, angles);
}

[[nodiscard]] Expr pi(const mathematics::MathRegistry& mathematics) {
    const auto* definition = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!definition)
        error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
    return Expr{definition->symbol};
}

// Fourier位相はセッションの既定角度に依存しない。常にradを明示する。
[[nodiscard]] Expr twiddle(
    std::size_t numerator,
    std::size_t denominator,
    bool inverse,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (denominator == 0)
        error::throwCalcError(error::CalcErrorType::Internal, "Fourier transform size is zero");

    BigInt signedNumerator = sizeInteger(numerator) * BigInt{2};
    if (!inverse)
        signedNumerator = -signedNumerator;
    const Rational coefficient{std::move(signedNumerator), sizeInteger(denominator)};

    Expr phase = exact::multiply(
        {Expr{Number{coefficient}}, pi(mathematics)}, registry, mathematics, angles);
    Expr radians = Expr::call(
        registry.symbol(BuiltinId::UnitApplied), {std::move(phase), Expr{std::string{"Rad"}}});
    return exact::call(BuiltinId::Cis, {std::move(radians)}, registry, mathematics, angles);
}

[[nodiscard]] std::vector<Expr> directTransform(
    const std::vector<Expr>& input,
    bool inverse,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::size_t n = input.size();
    if (n == 0)
        return {};

    std::vector<Expr> roots;
    roots.reserve(n);
    for (std::size_t m = 0; m < n; ++m)
        roots.push_back(twiddle(m, n, inverse, registry, mathematics, angles));

    std::vector<Expr> output;
    output.reserve(n);
    for (std::size_t k = 0; k < n; ++k) {
        Expr sum = zero();
        std::size_t rootIndex = 0;
        for (std::size_t j = 0; j < n; ++j) {
            sum = add(std::move(sum), multiply(input[j], roots[rootIndex],
                    registry, mathematics, angles), registry, mathematics, angles);

            if (j + 1 < n && k != 0) {
                if (rootIndex >= n - k)
                    rootIndex -= n - k;
                else
                    rootIndex += k;
            }
        }
        if (inverse)
            sum = divideBySize(std::move(sum), n, registry, mathematics, angles);
        output.push_back(std::move(sum));
    }
    return output;
}

[[nodiscard]] bool isPowerOfTwo(std::size_t value) noexcept {
    return value != 0 && (value & (value - 1)) == 0;
}

[[nodiscard]] std::vector<Expr> radix2Transform(
    const std::vector<Expr>& input,
    bool inverse,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::size_t n = input.size();
    if (!isPowerOfTwo(n))
        return directTransform(input, inverse, registry, mathematics, angles);

    std::vector<Expr> data = input;

    // in-place bit reversal permutation.
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; (j & bit) != 0; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            std::swap(data[i], data[j]);
    }

    for (std::size_t length = 2; length <= n; length <<= 1) {
        const std::size_t half = length >> 1;
        std::vector<Expr> roots;
        roots.reserve(half);
        for (std::size_t j = 0; j < half; ++j)
            roots.push_back(twiddle(j, length, inverse, registry, mathematics, angles));

        for (std::size_t block = 0; block < n; block += length) {
            for (std::size_t j = 0; j < half; ++j) {
                Expr even = data[block + j];
                Expr odd = multiply(data[block + j + half], roots[j],
                    registry, mathematics, angles);
                data[block + j] = add(even, odd, registry, mathematics, angles);
                data[block + j + half] = subtract(std::move(even), std::move(odd),
                    registry, mathematics, angles);
            }
        }

        if (length == n)
            break; // size_t overflow guard for the next shift.
    }

    if (inverse)
        for (Expr& value : data)
            value = divideBySize(std::move(value), n, registry, mathematics, angles);
    return data;
}

[[nodiscard]] Expr vectorExpr(std::vector<Expr> elements) {
    const std::size_t count = elements.size();
    return Expr::array({count}, std::move(elements));
}

} // namespace

Expr evaluateDft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::dft);
    return vectorExpr(directTransform(
        vectorArgument(arguments.front(), names::dft), false, registry, mathematics, angles));
}

Expr evaluateFft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::fft);
    return vectorExpr(radix2Transform(
        vectorArgument(arguments.front(), names::fft), false, registry, mathematics, angles));
}

Expr evaluateIfft(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 1, names::ifft);
    return vectorExpr(radix2Transform(
        vectorArgument(arguments.front(), names::ifft), true, registry, mathematics, angles));
}

Expr evaluateConvolution(
    std::span<const Expr> arguments,
    const evaluation::BuiltinRegistry& registry,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    requireArity(arguments, 2, names::convolve);
    const std::vector<Expr> lhs = vectorArgument(arguments[0], names::convolve);
    const std::vector<Expr> rhs = vectorArgument(arguments[1], names::convolve);
    if (lhs.empty() || rhs.empty())
        return vectorExpr({});

    if (lhs.size() > std::numeric_limits<std::size_t>::max() - rhs.size() + 1)
        error::throwCalcError(error::CalcErrorType::Overflow, "convolution result is too large");

    std::vector<Expr> output(lhs.size() + rhs.size() - 1, zero());
    for (std::size_t i = 0; i < lhs.size(); ++i)
        for (std::size_t j = 0; j < rhs.size(); ++j)
            output[i + j] = add(std::move(output[i + j]),
                multiply(lhs[i], rhs[j], registry, mathematics, angles),
                registry, mathematics, angles);
    return vectorExpr(std::move(output));
}

} // namespace mmcal::builtins
