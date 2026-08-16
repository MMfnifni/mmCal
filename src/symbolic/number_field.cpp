// 代数体Q(theta)とそのpower basis上のexact要素
#include "number_field.hpp"

#include "numeric/big_int.hpp"
#include "rational_linear_basis.hpp"

#include <algorithm>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <utility>

namespace mmcal::symbolic {
namespace {

using numeric::BigInt;
using numeric::Rational;
using Polynomial = std::vector<Rational>;
using RationalVector = std::vector<Rational>;
using RationalMatrix = std::vector<RationalVector>;

constexpr std::size_t maximumSignRefinementBits = 4096;
constexpr std::size_t maximumInternedNumberFields = 256;
constexpr std::size_t maximumCachedReciprocalsPerField = 16;
constexpr std::size_t maximumCachedMinimalPolynomialsPerField = 16;

class NumberFieldInterner final {
public:
    [[nodiscard]] std::shared_ptr<const NumberFieldContext> find(
        const AlgebraicNumber& generator) {
        std::lock_guard lock(mutex_);
        pruneExpired();
        return findLocked(generator);
    }

    [[nodiscard]] std::shared_ptr<const NumberFieldContext> intern(
        std::shared_ptr<const NumberFieldContext> candidate) {
        std::lock_guard lock(mutex_);
        pruneExpired();
        if (auto existing = findLocked(candidate->generator()))
            return existing;

        if (fields_.size() >= maximumInternedNumberFields)
            fields_.erase(fields_.begin());
        fields_.push_back(candidate);
        return candidate;
    }

private:
    std::mutex mutex_;
    std::vector<std::weak_ptr<const NumberFieldContext>> fields_;

    void pruneExpired() {
        fields_.erase(
            std::remove_if(fields_.begin(), fields_.end(),
                [](const std::weak_ptr<const NumberFieldContext>& field) {
                    return field.expired();
                }),
            fields_.end());
    }

    [[nodiscard]] std::shared_ptr<const NumberFieldContext> findLocked(
        const AlgebraicNumber& generator) {
        for (std::size_t i = 0; i < fields_.size(); ++i) {
            auto existing = fields_[i].lock();
            if (!existing || !existing->generator().hasSameRootIdentity(generator))
                continue;

            // bounded LRU。cache policyは性能だけに影響し，数学的identityには影響しない。
            if (i + 1 != fields_.size()) {
                auto hit = std::move(fields_[i]);
                fields_.erase(fields_.begin() + static_cast<std::ptrdiff_t>(i));
                fields_.push_back(std::move(hit));
            }
            return existing;
        }
        return {};
    }
};

[[nodiscard]] NumberFieldInterner& numberFieldInterner() {
    static NumberFieldInterner interner;
    return interner;
}

struct RationalInterval final {
    Rational lower;
    Rational upper;
};

[[nodiscard]] RationalInterval multiplyIntervals(
    const RationalInterval& lhs,
    const RationalInterval& rhs) {
    const Rational p0 = lhs.lower * rhs.lower;
    const Rational p1 = lhs.lower * rhs.upper;
    const Rational p2 = lhs.upper * rhs.lower;
    const Rational p3 = lhs.upper * rhs.upper;
    Rational lower = p0;
    Rational upper = p0;
    for (const Rational* value : {&p1, &p2, &p3}) {
        if (*value < lower) lower = *value;
        if (upper < *value) upper = *value;
    }
    return RationalInterval{std::move(lower), std::move(upper)};
}

[[nodiscard]] RationalInterval evaluateInterval(
    std::span<const Rational> coefficients,
    const RationalRootInterval& theta) {
    if (coefficients.empty())
        return RationalInterval{};

    RationalInterval result{coefficients.back(), coefficients.back()};
    const RationalInterval variable{theta.lower, theta.upper};
    for (std::size_t i = coefficients.size() - 1; i-- > 0;) {
        result = multiplyIntervals(result, variable);
        result.lower += coefficients[i];
        result.upper += coefficients[i];
    }
    return result;
}

void trim(Polynomial& polynomial) {
    while (!polynomial.empty() && polynomial.back().isZero())
        polynomial.pop_back();
}

[[nodiscard]] Polynomial polynomialRemainder(
    Polynomial dividend,
    const Polynomial& divisor) {
    if (divisor.empty())
        throw std::logic_error("Number field polynomial divisor is zero");

    trim(dividend);
    const std::size_t divisorDegree = divisor.size() - 1;
    const Rational divisorLeading = divisor.back();
    while (!dividend.empty() && dividend.size() - 1 >= divisorDegree) {
        const std::size_t shift = dividend.size() - 1 - divisorDegree;
        const Rational factor = dividend.back() / divisorLeading;
        for (std::size_t i = 0; i <= divisorDegree; ++i)
            dividend[shift + i] -= factor * divisor[i];
        trim(dividend);
    }
    return dividend;
}

[[nodiscard]] std::pair<Polynomial, Polynomial> polynomialDivide(
    Polynomial dividend,
    const Polynomial& divisor) {
    if (divisor.empty())
        throw std::logic_error("Number field polynomial divisor is zero");

    trim(dividend);
    const std::size_t divisorDegree = divisor.size() - 1;
    const Rational divisorLeading = divisor.back();
    Polynomial quotient;
    if (!dividend.empty() && dividend.size() - 1 >= divisorDegree)
        quotient.resize(dividend.size() - divisorDegree);

    while (!dividend.empty() && dividend.size() - 1 >= divisorDegree) {
        const std::size_t shift = dividend.size() - 1 - divisorDegree;
        const Rational factor = dividend.back() / divisorLeading;
        quotient[shift] += factor;
        for (std::size_t i = 0; i <= divisorDegree; ++i)
            dividend[shift + i] -= factor * divisor[i];
        trim(dividend);
    }
    trim(quotient);
    return {std::move(quotient), std::move(dividend)};
}

[[nodiscard]] Polynomial polynomialSubtract(
    Polynomial lhs,
    const Polynomial& rhs) {
    if (lhs.size() < rhs.size())
        lhs.resize(rhs.size());
    for (std::size_t i = 0; i < rhs.size(); ++i)
        lhs[i] -= rhs[i];
    trim(lhs);
    return lhs;
}

[[nodiscard]] Polynomial polynomialMultiply(
    const Polynomial& lhs,
    const Polynomial& rhs) {
    if (lhs.empty() || rhs.empty())
        return {};
    Polynomial result(lhs.size() + rhs.size() - 1);
    for (std::size_t i = 0; i < lhs.size(); ++i)
        for (std::size_t j = 0; j < rhs.size(); ++j)
            result[i + j] += lhs[i] * rhs[j];
    trim(result);
    return result;
}

} // namespace

NumberFieldContext::NumberFieldContext(
    AlgebraicNumber generator,
    std::vector<Rational> minimalPolynomial,
    std::vector<Rational> reduction)
    : generator_(std::move(generator)),
      minimalPolynomial_(std::move(minimalPolynomial)),
      reduction_(std::move(reduction)) {}

std::shared_ptr<const NumberFieldContext> NumberFieldContext::create(
    AlgebraicNumber generator) {
    // Contextのgenerator identityに別fieldの算術cacheを抱え込ませない。
    generator = generator.withArithmeticElement({});
    const std::span<const Rational> definingPolynomial = generator.polynomial();
    if (definingPolynomial.size() <= 1
        || definingPolynomial.back() != Rational{BigInt{1}})
        return {};

    if (auto existing = numberFieldInterner().find(generator))
        return existing;

    std::vector<Rational> polynomial(
        definingPolynomial.begin(), definingPolynomial.end());
    std::vector<Rational> reduction(polynomial.size() - 1);
    for (std::size_t i = 0; i < reduction.size(); ++i)
        reduction[i] = -polynomial[i];
    auto candidate = std::shared_ptr<const NumberFieldContext>{new NumberFieldContext{
        std::move(generator), std::move(polynomial), std::move(reduction)}};
    return numberFieldInterner().intern(std::move(candidate));
}

const AlgebraicNumber& NumberFieldContext::generator() const noexcept { return generator_; }
AlgebraicRootDomain NumberFieldContext::domain() const noexcept { return generator_.domain(); }
std::span<const Rational> NumberFieldContext::minimalPolynomial() const noexcept {
    return minimalPolynomial_;
}
std::span<const Rational> NumberFieldContext::reduction() const noexcept { return reduction_; }
std::size_t NumberFieldContext::degree() const noexcept { return minimalPolynomial_.size() - 1; }

std::vector<Rational> NumberFieldContext::multiply(
    std::span<const Rational> lhs,
    std::span<const Rational> rhs) const {
    const std::size_t fieldDegree = degree();
    if (lhs.size() != fieldDegree || rhs.size() != fieldDegree)
        throw std::logic_error("Number field element size mismatch");

    Polynomial product(2 * fieldDegree - 1);
    for (std::size_t i = 0; i < fieldDegree; ++i)
        for (std::size_t j = 0; j < fieldDegree; ++j)
            product[i + j] += lhs[i] * rhs[j];

    for (std::int64_t exponent = static_cast<std::int64_t>(product.size()) - 1;
         exponent >= static_cast<std::int64_t>(fieldDegree); --exponent) {
        const std::size_t index = static_cast<std::size_t>(exponent);
        Rational coefficient = product[index];
        if (coefficient.isZero())
            continue;
        product[index] = Rational{};
        const std::size_t shift = index - fieldDegree;
        for (std::size_t k = 0; k < fieldDegree; ++k)
            product[shift + k] += coefficient * reduction_[k];
    }
    product.resize(fieldDegree);
    return product;
}

std::optional<std::vector<Rational>> NumberFieldContext::findCachedReciprocal(
    std::span<const Rational> value) const {
    std::lock_guard lock(reciprocalCacheMutex_);
    for (std::size_t i = 0; i < reciprocalCache_.size(); ++i) {
        const ReciprocalCacheEntry& entry = reciprocalCache_[i];
        const bool direct = std::equal(
            value.begin(), value.end(), entry.value.begin(), entry.value.end());
        const bool reverse = !direct && std::equal(
            value.begin(), value.end(), entry.reciprocal.begin(), entry.reciprocal.end());
        if (!direct && !reverse)
            continue;

        std::vector<Rational> result = direct ? entry.reciprocal : entry.value;
        if (i + 1 != reciprocalCache_.size()) {
            ReciprocalCacheEntry hit = std::move(reciprocalCache_[i]);
            reciprocalCache_.erase(reciprocalCache_.begin() + static_cast<std::ptrdiff_t>(i));
            reciprocalCache_.push_back(std::move(hit));
        }
        return result;
    }
    return std::nullopt;
}

std::optional<std::vector<Rational>> NumberFieldContext::reciprocalUncached(
    std::span<const Rational> value) const {
    const std::size_t fieldDegree = degree();

    // Qは全てのnumber fieldの部分体なので，定数要素は多項式Euclidを回さない。
    bool constant = true;
    for (std::size_t i = 1; i < value.size(); ++i)
        constant = constant && value[i].isZero();
    if (constant) {
        if (value.front().isZero())
            return std::nullopt;
        std::vector<Rational> result(fieldDegree);
        result.front() = Rational{BigInt{1}} / value.front();
        return result;
    }

    Polynomial oldR(minimalPolynomial_.begin(), minimalPolynomial_.end());
    Polynomial r(value.begin(), value.end());
    trim(r);
    Polynomial oldT;
    Polynomial t{Rational{BigInt{1}}};

    while (!r.empty()) {
        auto division = polynomialDivide(std::move(oldR), r);
        Polynomial quotient = std::move(division.first);
        Polynomial remainder = std::move(division.second);
        oldR = std::move(r);
        r = std::move(remainder);

        Polynomial nextT = polynomialSubtract(
            std::move(oldT), polynomialMultiply(quotient, t));
        oldT = std::move(t);
        t = std::move(nextT);
    }

    if (oldR.size() != 1 || oldR.front().isZero())
        return std::nullopt;
    const Rational inverseGcd = Rational{BigInt{1}} / oldR.front();
    for (Rational& coefficient : oldT)
        coefficient *= inverseGcd;
    oldT = polynomialRemainder(std::move(oldT), minimalPolynomial_);
    oldT.resize(fieldDegree);
    return oldT;
}

void NumberFieldContext::publishReciprocal(
    std::vector<Rational> value,
    std::vector<Rational> reciprocal) const {
    std::lock_guard lock(reciprocalCacheMutex_);

    // reciprocal pairは向きを持たない。競合した計算が先にpublishしていれば既存entryを優先する。
    for (std::size_t i = 0; i < reciprocalCache_.size(); ++i) {
        const ReciprocalCacheEntry& entry = reciprocalCache_[i];
        const bool samePair = (entry.value == value && entry.reciprocal == reciprocal)
            || (entry.value == reciprocal && entry.reciprocal == value);
        if (!samePair)
            continue;
        if (i + 1 != reciprocalCache_.size()) {
            ReciprocalCacheEntry hit = std::move(reciprocalCache_[i]);
            reciprocalCache_.erase(reciprocalCache_.begin() + static_cast<std::ptrdiff_t>(i));
            reciprocalCache_.push_back(std::move(hit));
        }
        return;
    }

    if (reciprocalCache_.size() >= maximumCachedReciprocalsPerField)
        reciprocalCache_.erase(reciprocalCache_.begin());
    reciprocalCache_.push_back(ReciprocalCacheEntry{std::move(value), std::move(reciprocal)});
}

std::optional<std::vector<Rational>> NumberFieldContext::reciprocal(
    std::span<const Rational> value) const {
    const std::size_t fieldDegree = degree();
    if (value.size() != fieldDegree)
        return std::nullopt;
    if (std::all_of(value.begin(), value.end(),
            [](const Rational& coefficient) { return coefficient.isZero(); }))
        return std::nullopt;

    if (auto cached = findCachedReciprocal(value))
        return cached;

    std::vector<Rational> key(value.begin(), value.end());
    auto result = reciprocalUncached(value);
    if (!result)
        return std::nullopt;

    // 1 entryを双方向のreciprocal pairとして保持する。
    publishReciprocal(std::move(key), *result);
    return result;
}

std::optional<std::vector<Rational>> NumberFieldContext::findCachedMinimalPolynomial(
    std::span<const Rational> value) const {
    std::lock_guard lock(minimalPolynomialCacheMutex_);
    for (std::size_t i = 0; i < minimalPolynomialCache_.size(); ++i) {
        const MinimalPolynomialCacheEntry& entry = minimalPolynomialCache_[i];
        if (!std::equal(value.begin(), value.end(), entry.value.begin(), entry.value.end()))
            continue;

        std::vector<Rational> result = entry.polynomial;
        if (i + 1 != minimalPolynomialCache_.size()) {
            MinimalPolynomialCacheEntry hit = std::move(minimalPolynomialCache_[i]);
            minimalPolynomialCache_.erase(
                minimalPolynomialCache_.begin() + static_cast<std::ptrdiff_t>(i));
            minimalPolynomialCache_.push_back(std::move(hit));
        }
        return result;
    }
    return std::nullopt;
}

std::optional<std::vector<Rational>> NumberFieldContext::minimalPolynomialUncached(
    std::span<const Rational> value) const {
    if (value.size() != degree())
        return std::nullopt;

    bool rational = true;
    for (std::size_t i = 1; i < value.size(); ++i)
        rational = rational && value[i].isZero();
    if (rational)
        return Polynomial{-value[0], Rational{BigInt{1}}};

    bool generatorCoordinates = value.size() > 1 && value[0].isZero()
        && value[1] == Rational{BigInt{1}};
    for (std::size_t i = 2; i < value.size(); ++i)
        generatorCoordinates = generatorCoordinates && value[i].isZero();
    if (generatorCoordinates)
        return minimalPolynomial_;

    RationalVector current(degree());
    current.front() = Rational{BigInt{1}};

    // 1,a,a^2,...を逐次消去し，最初の線形従属をminimal polynomialとする。
    detail::RationalLinearBasis krylovBasis(degree());
    if (krylovBasis.append(current))
        return std::nullopt;

    for (std::size_t relationDegree = 1; relationDegree <= degree(); ++relationDegree) {
        current = multiply(current, value);
        if (auto relation = krylovBasis.append(current))
            return relation;
    }
    return std::nullopt;
}

void NumberFieldContext::publishMinimalPolynomial(
    std::vector<Rational> value,
    std::vector<Rational> polynomial) const {
    std::lock_guard lock(minimalPolynomialCacheMutex_);
    for (std::size_t i = 0; i < minimalPolynomialCache_.size(); ++i) {
        if (minimalPolynomialCache_[i].value != value)
            continue;
        if (i + 1 != minimalPolynomialCache_.size()) {
            MinimalPolynomialCacheEntry hit = std::move(minimalPolynomialCache_[i]);
            minimalPolynomialCache_.erase(
                minimalPolynomialCache_.begin() + static_cast<std::ptrdiff_t>(i));
            minimalPolynomialCache_.push_back(std::move(hit));
        }
        return;
    }

    if (minimalPolynomialCache_.size() >= maximumCachedMinimalPolynomialsPerField)
        minimalPolynomialCache_.erase(minimalPolynomialCache_.begin());
    minimalPolynomialCache_.push_back(
        MinimalPolynomialCacheEntry{std::move(value), std::move(polynomial)});
}

std::optional<std::vector<Rational>> NumberFieldContext::minimalPolynomialOf(
    std::span<const Rational> value) const {
    if (value.size() != degree())
        return std::nullopt;
    if (auto cached = findCachedMinimalPolynomial(value))
        return cached;

    auto polynomial = minimalPolynomialUncached(value);
    if (!polynomial)
        return std::nullopt;

    std::vector<Rational> coordinates(value.begin(), value.end());
    std::vector<Rational> result = *polynomial;
    publishMinimalPolynomial(std::move(coordinates), std::move(*polynomial));
    return result;
}


AlgebraicElement::AlgebraicElement(
    std::shared_ptr<const NumberFieldContext> field,
    std::vector<Rational> coefficients)
    : field_(std::move(field)), coefficients_(std::move(coefficients)) {}

std::optional<AlgebraicElement> AlgebraicElement::create(
    std::shared_ptr<const NumberFieldContext> field,
    std::vector<Rational> coefficients) {
    if (!field || coefficients.size() > field->degree())
        return std::nullopt;
    coefficients.resize(field->degree());
    return AlgebraicElement{std::move(field), std::move(coefficients)};
}

std::optional<AlgebraicElement> AlgebraicElement::generator(
    std::shared_ptr<const NumberFieldContext> field) {
    if (!field || field->degree() == 0)
        return std::nullopt;

    std::vector<Rational> coefficients(field->degree());
    if (field->degree() == 1)
        coefficients[0] = field->reduction()[0];
    else
        coefficients[1] = Rational{BigInt{1}};
    return create(std::move(field), std::move(coefficients));
}


const std::shared_ptr<const NumberFieldContext>& AlgebraicElement::field() const noexcept {
    return field_;
}
std::span<const Rational> AlgebraicElement::coefficients() const noexcept { return coefficients_; }
bool AlgebraicElement::isZero() const noexcept {
    return std::all_of(coefficients_.begin(), coefficients_.end(),
        [](const Rational& coefficient) { return coefficient.isZero(); });
}

std::optional<bool> AlgebraicElement::exactEquals(
    const AlgebraicElement& rhs) const noexcept {
    if (!hasSameField(rhs))
        return std::nullopt;
    return coefficients_ == rhs.coefficients_;
}

std::optional<RationalRootInterval> AlgebraicElement::refinedRealInterval(
    std::size_t precisionBits) const {
    if (field_->domain() != AlgebraicRootDomain::Real)
        return std::nullopt;
    const RealAlgebraicNumber* generator = field_->generator().asReal();
    if (!generator)
        return std::nullopt;

    const RationalInterval value = evaluateInterval(
        coefficients_, generator->refined(precisionBits));
    return RationalRootInterval{std::move(value.lower), std::move(value.upper)};
}

std::optional<AlgebraicSign> AlgebraicElement::exactSign() const {
    if (field_->domain() != AlgebraicRootDomain::Real)
        return std::nullopt;
    if (isZero())
        return AlgebraicSign::Zero;

    const RealAlgebraicNumber* generator = field_->generator().asReal();
    if (!generator)
        return std::nullopt;

    for (std::size_t bits = 32; bits <= maximumSignRefinementBits; bits *= 2) {
        const RationalInterval value = evaluateInterval(
            coefficients_, generator->refined(bits));
        if (value.upper < Rational{})
            return AlgebraicSign::Negative;
        if (Rational{} < value.lower)
            return AlgebraicSign::Positive;
    }
    return std::nullopt;
}

bool AlgebraicElement::hasSameField(const AlgebraicElement& rhs) const noexcept {
    return field_.get() == rhs.field_.get();
}

std::optional<AlgebraicElement> AlgebraicElement::add(const AlgebraicElement& rhs) const {
    if (!hasSameField(rhs))
        return std::nullopt;
    RationalVector result(coefficients_.size());
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = coefficients_[i] + rhs.coefficients_[i];
    return AlgebraicElement{field_, std::move(result)};
}

std::optional<AlgebraicElement> AlgebraicElement::subtract(const AlgebraicElement& rhs) const {
    if (!hasSameField(rhs))
        return std::nullopt;
    RationalVector result(coefficients_.size());
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = coefficients_[i] - rhs.coefficients_[i];
    return AlgebraicElement{field_, std::move(result)};
}

std::optional<AlgebraicElement> AlgebraicElement::multiply(const AlgebraicElement& rhs) const {
    if (!hasSameField(rhs))
        return std::nullopt;
    return AlgebraicElement{field_, field_->multiply(coefficients_, rhs.coefficients_)};
}

std::optional<AlgebraicElement> AlgebraicElement::divide(const AlgebraicElement& rhs) const {
    if (!hasSameField(rhs))
        return std::nullopt;
    const auto inverse = field_->reciprocal(rhs.coefficients_);
    if (!inverse)
        return std::nullopt;
    return AlgebraicElement{field_, field_->multiply(coefficients_, *inverse)};
}

std::optional<std::vector<Rational>> AlgebraicElement::minimalPolynomial() const {
    return field_->minimalPolynomialOf(coefficients_);
}

} // namespace mmcal::symbolic
