// 根号の簡約
#include "exact_roots.hpp"

#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/rational.hpp"

#include <cstdint>
#include <stdexcept>

namespace mmcal::mathematics {
namespace {

using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] RealNumber half(const RealNumber& value) {
    return value / RealNumber{BigInt{2}};
}

struct IntegerSquareFactor final {
    BigInt outside{1};
    BigInt inside{1};
};

[[nodiscard]] bool isPrimeTrial(std::uint32_t value) noexcept {
    if (value < 2)
        return false;
    if ((value & 1U) == 0U)
        return value == 2;
    for (std::uint32_t divisor = 3;
         static_cast<std::uint64_t>(divisor) * divisor <= value;
         divisor += 2) {
        if (value % divisor == 0)
            return false;
    }
    return true;
}

[[nodiscard]] IntegerSquareFactor extractKnownSquareFactor(const BigInt& positive) {
    if (positive.isNegative() || positive.isZero())
        throw std::invalid_argument("Square-factor decomposition requires a positive integer");

    BigInt remainder = positive;
    IntegerSquareFactor result;

    // まず全体が完全平方かを調べる。
    // 巨大な perfect square に対して素因数試行を先に何千回も行う必要がなく、かつ integerSqrt はすでに任意精度整数でexactである。
    const auto wholeRoot = numeric::integerSqrt(remainder);
    if (wholeRoot.remainder.isZero())
        return IntegerSquareFactor{wholeRoot.root, BigInt{1}};

    // 一般BigIntの完全因数分解をsqrt(n)まで試すと、巨大な素数入力で計算不能になる。
    // canonicalizerでは、小さい入力ほど十分広く小素数を剥がし、巨大入力では試行範囲を意図的に抑える。
    // 残部全体が完全平方なら最後にもう一度integerSqrtで回収する。
    // したがって難しい巨大合成数で平方因子を取り切れない場合はあっても、誤った因数をradicalの外へ出すことはない。
    const std::size_t bits = positive.bitLength();
    const std::uint32_t trialLimit = bits <= 128 ? 10000U : bits <= 512 ? 1000U : 100U;
    for (std::uint32_t prime = 2; prime <= trialLimit; ++prime) {
        if (!isPrimeTrial(prime))
            continue;

        const BigInt divisor{static_cast<std::int64_t>(prime)};
        std::size_t exponent = 0;
        while ((remainder % divisor).isZero()) {
            remainder /= divisor;
            ++exponent;
        }
        for (std::size_t i = 0; i < exponent / 2; ++i)
            result.outside *= divisor;
        if ((exponent & 1U) != 0U)
            result.inside *= divisor;
        if (remainder == BigInt{1})
            return result;
    }

    const auto residualRoot = numeric::integerSqrt(remainder);
    if (residualRoot.remainder.isZero())
        result.outside *= residualRoot.root;
    else
        result.inside *= remainder;
    return result;
}

} // namespace

RationalSquareRootDecomposition decomposePositiveRationalSquareRoot(
    const Rational& value) {
    if (value <= Rational{BigInt{0}})
        throw std::invalid_argument("Rational square-root decomposition requires a positive value");

    const IntegerSquareFactor numerator = extractKnownSquareFactor(value.numerator());
    const IntegerSquareFactor denominator = extractKnownSquareFactor(value.denominator());

    // sqrt(in/id) = sqrt(in*id) / id。したがって分母の残radicalを有理化し、外側係数 on/od と合わせる。
    const BigInt coefficientDenominator = denominator.outside * denominator.inside;
    return RationalSquareRootDecomposition{
        Rational{numerator.outside, coefficientDenominator},
        numerator.inside * denominator.inside};
}

std::optional<RealNumber> exactSquareRoot(const RealNumber& value) {
    if (value.isNegative())
        return std::nullopt;

    if (value.isInteger()) {
        const auto result = numeric::integerSqrt(value.asInteger());
        if (!result.remainder.isZero())
            return std::nullopt;
        return RealNumber{result.root};
    }

    const Rational& fraction = value.asRational();
    const auto numeratorRoot = numeric::integerSqrt(fraction.numerator());
    const auto denominatorRoot = numeric::integerSqrt(fraction.denominator());

    if (!numeratorRoot.remainder.isZero() || !denominatorRoot.remainder.isZero())
        return std::nullopt;

    return RealNumber{Rational{numeratorRoot.root, denominatorRoot.root}};
}


std::optional<RealNumber> exactRealCubeRoot(const RealNumber& value) {
    if (value.isZero())
        return RealNumber{};

    const bool negative = value.isNegative();
    const Rational magnitude = value.abs().toRational();
    const auto numeratorRoot = numeric::integerCubeRoot(magnitude.numerator());
    const auto denominatorRoot = numeric::integerCubeRoot(magnitude.denominator());
    if (!numeratorRoot.remainder.isZero() || !denominatorRoot.remainder.isZero())
        return std::nullopt;

    Rational root{numeratorRoot.root, denominatorRoot.root};
    if (negative)
        root = -root;
    return RealNumber{std::move(root)};
}

std::optional<Number> exactPrincipalSquareRoot(const Number& value) {
    // 実数は従来の平方根規則をそのまま使う。
    if (value.isReal()) {
        const RealNumber& real = value.asReal();
        if (!real.isNegative()) {
            const auto root = exactSquareRoot(real);
            return root ? std::optional<Number>{Number{*root}} : std::nullopt;
        }

        const auto root = exactSquareRoot(real.abs());
        return root
            ? std::optional<Number>{Number::complex(RealNumber{}, *root)}
            : std::nullopt;
    }

    // z = a + bI のprincipal square rootを u + vI とする。
    //
    //   u^2 - v^2 = a
    //   2uv         = b
    //   u^2 + v^2   = |z| = sqrt(a^2+b^2)
    //
    // よって
    //
    //   u^2 = (|z| + a)/2
    //   v^2 = (|z| - a)/2
    //
    // となる。現在のNumberは実部・虚部がRationalまでなので、|z|, u, |v|が
    // すべてexactなRealNumberとして閉じる場合だけNumberへ畳み込む。
    const auto& complex = value.asComplex();
    const RealNumber& a = complex.real;
    const RealNumber& b = complex.imaginary;

    const RealNumber magnitudeSquared = a * a + b * b;
    const auto magnitude = exactSquareRoot(magnitudeSquared);
    if (!magnitude)
        return std::nullopt;

    const RealNumber uSquared = half(*magnitude + a);
    const RealNumber vSquared = half(*magnitude - a);
    const auto u = exactSquareRoot(uSquared);
    const auto vMagnitude = exactSquareRoot(vSquared);
    if (!u || !vMagnitude)
        return std::nullopt;

    // principal branchは Re(sqrt(z)) >= 0。
    // b != 0 のとき Im(sqrt(z)) は b と同符号。
    // b == 0 かつ a < 0 の負実軸上では +I*sqrt(-a) を採る。
    RealNumber v = *vMagnitude;
    if (b.isNegative())
        v = -v;

    return Number::complex(*u, std::move(v));
}

} // namespace mmcal::mathematics
