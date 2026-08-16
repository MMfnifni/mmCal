// exact Fourier用cyclotomic quotient field Q[t]/Phi_n(t)
#include "cyclotomic_field.hpp"

#include "numeric/big_int.hpp"

#include <cstddef>
#include <memory>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using numeric::BigInt;
using numeric::Rational;
using Polynomial = std::vector<Rational>;

constexpr std::size_t maximumCyclotomicDegree = 64;

void trim(Polynomial& polynomial) {
    while (!polynomial.empty() && polynomial.back().isZero())
        polynomial.pop_back();
}

[[nodiscard]] std::optional<Polynomial> divideExactMonic(
    Polynomial dividend,
    const Polynomial& divisor) {
    if (divisor.empty() || divisor.back() != Rational{BigInt{1}})
        return std::nullopt;
    trim(dividend);
    if (dividend.size() < divisor.size())
        return std::nullopt;

    Polynomial quotient(dividend.size() - divisor.size() + 1);
    while (dividend.size() >= divisor.size()) {
        const std::size_t shift = dividend.size() - divisor.size();
        const Rational factor = dividend.back();
        quotient[shift] = factor;
        for (std::size_t i = 0; i < divisor.size(); ++i)
            dividend[shift + i] -= factor * divisor[i];
        trim(dividend);
    }
    if (!dividend.empty())
        return std::nullopt;
    trim(quotient);
    return quotient;
}

[[nodiscard]] std::optional<Polynomial> cyclotomicPolynomialImpl(
    std::size_t n,
    std::unordered_map<std::size_t, Polynomial>& cache) {
    if (const auto iterator = cache.find(n); iterator != cache.end())
        return iterator->second;
    if (n == 0)
        return std::nullopt;
    if (n == 1) {
        Polynomial polynomial{Rational{BigInt{-1}}, Rational{BigInt{1}}};
        cache.emplace(n, polynomial);
        return polynomial;
    }

    Polynomial polynomial(n + 1);
    polynomial[0] = Rational{BigInt{-1}};
    polynomial[n] = Rational{BigInt{1}};
    for (std::size_t divisor = 1; divisor < n; ++divisor) {
        if (n % divisor != 0)
            continue;
        const auto factor = cyclotomicPolynomialImpl(divisor, cache);
        if (!factor)
            return std::nullopt;
        auto quotient = divideExactMonic(std::move(polynomial), *factor);
        if (!quotient)
            return std::nullopt;
        polynomial = std::move(*quotient);
    }
    cache.emplace(n, polynomial);
    return polynomial;
}

[[nodiscard]] std::optional<Polynomial> cyclotomicPolynomial(std::size_t n) {
    std::unordered_map<std::size_t, Polynomial> cache;
    return cyclotomicPolynomialImpl(n, cache);
}

[[nodiscard]] std::vector<Rational> multiplyReduced(
    std::span<const Rational> lhs,
    std::span<const Rational> rhs,
    std::span<const Rational> reduction) {
    const std::size_t degree = reduction.size();
    if (lhs.size() != degree || rhs.size() != degree || degree == 0)
        return {};

    std::vector<Rational> product(2 * degree - 1);
    for (std::size_t i = 0; i < degree; ++i)
        for (std::size_t j = 0; j < degree; ++j)
            product[i + j] += lhs[i] * rhs[j];

    for (std::size_t power = product.size(); power-- > degree;) {
        const Rational factor = product[power];
        if (factor.isZero())
            continue;
        const std::size_t shift = power - degree;
        for (std::size_t i = 0; i < degree; ++i)
            product[shift + i] += factor * reduction[i];
    }
    product.resize(degree);
    return product;
}

} // namespace

std::size_t cyclotomicDegree(std::size_t conductor) noexcept {
    if (conductor == 0)
        return 0;
    std::size_t result = conductor;
    std::size_t remaining = conductor;
    for (std::size_t prime = 2; prime <= remaining / prime; ++prime) {
        if (remaining % prime != 0)
            continue;
        result -= result / prime;
        while (remaining % prime == 0)
            remaining /= prime;
    }
    if (remaining > 1)
        result -= result / remaining;
    return result;
}

CyclotomicFieldContext::CyclotomicFieldContext(
    std::size_t conductor,
    std::vector<Rational> polynomial,
    std::vector<Rational> reduction,
    std::vector<std::vector<Rational>> powers)
    : conductor_(conductor),
      polynomial_(std::move(polynomial)),
      reduction_(std::move(reduction)),
      powers_(std::move(powers)) {}

std::shared_ptr<const CyclotomicFieldContext> CyclotomicFieldContext::create(
    std::size_t conductor) {
    const std::size_t expectedDegree = cyclotomicDegree(conductor);
    if (conductor == 0 || expectedDegree == 0 || expectedDegree > maximumCyclotomicDegree)
        return {};
    const auto polynomial = cyclotomicPolynomial(conductor);
    if (!polynomial || polynomial->size() != expectedDegree + 1)
        return {};

    const std::size_t degree = expectedDegree;
    std::vector<Rational> reduction(degree);
    for (std::size_t i = 0; i < degree; ++i)
        reduction[i] = -(*polynomial)[i];

    std::vector<std::vector<Rational>> powers(conductor);
    powers[0].assign(degree, Rational{});
    powers[0][0] = Rational{BigInt{1}};
    if (conductor > 1) {
        powers[1].assign(degree, Rational{});
        if (degree == 1)
            powers[1][0] = reduction[0];
        else
            powers[1][1] = Rational{BigInt{1}};
    }
    for (std::size_t exponent = 2; exponent < conductor; ++exponent)
        powers[exponent] = multiplyReduced(powers[exponent - 1], powers[1], reduction);

    return std::shared_ptr<const CyclotomicFieldContext>{new CyclotomicFieldContext{
        conductor, *polynomial, std::move(reduction), std::move(powers)}};
}

std::size_t CyclotomicFieldContext::conductor() const noexcept { return conductor_; }
std::size_t CyclotomicFieldContext::degree() const noexcept { return polynomial_.size() - 1; }
std::span<const Rational> CyclotomicFieldContext::minimalPolynomial() const noexcept {
    return polynomial_;
}
std::span<const Rational> CyclotomicFieldContext::power(std::size_t exponent) const {
    return powers_[exponent % conductor_];
}

std::vector<Rational> CyclotomicFieldContext::multiply(
    std::span<const Rational> lhs,
    std::span<const Rational> rhs) const {
    return multiplyReduced(lhs, rhs, reduction_);
}

std::optional<std::vector<Rational>> CyclotomicFieldContext::embedGaussianRational(
    const Rational& real,
    const Rational& imaginary) const {
    std::vector<Rational> coordinates(degree());
    coordinates[0] = real;
    if (imaginary.isZero())
        return coordinates;
    if (conductor_ % 4 != 0)
        return std::nullopt;

    // t=exp(-2 Pi I/n) なので I=t^(3n/4)。
    const auto imaginaryUnit = power((3 * conductor_) / 4);
    for (std::size_t i = 0; i < coordinates.size(); ++i)
        coordinates[i] += imaginary * imaginaryUnit[i];
    return coordinates;
}

std::optional<std::pair<Rational, Rational>> CyclotomicFieldContext::exactGaussianRational(
    std::span<const Rational> value) const {
    if (value.size() != degree())
        return std::nullopt;
    if (conductor_ % 4 != 0) {
        for (std::size_t i = 1; i < value.size(); ++i)
            if (!value[i].isZero())
                return std::nullopt;
        return std::pair<Rational, Rational>{value[0], Rational{}};
    }

    const auto imaginaryUnit = power((3 * conductor_) / 4);
    std::optional<Rational> imaginary;
    for (std::size_t i = 1; i < value.size(); ++i) {
        if (imaginaryUnit[i].isZero()) {
            if (!value[i].isZero())
                return std::nullopt;
            continue;
        }
        const Rational candidate = value[i] / imaginaryUnit[i];
        if (!imaginary)
            imaginary = candidate;
        else if (*imaginary != candidate)
            return std::nullopt;
    }
    const Rational imag = imaginary.value_or(Rational{});
    for (std::size_t i = 1; i < value.size(); ++i)
        if (value[i] != imag * imaginaryUnit[i])
            return std::nullopt;
    const Rational real = value[0] - imag * imaginaryUnit[0];
    return std::pair<Rational, Rational>{real, imag};
}

} // namespace mmcal::symbolic
