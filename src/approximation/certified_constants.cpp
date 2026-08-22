// Piなど定数の保証付き評価
#include "certified_constants.hpp"
#include "certified_precision.hpp"
#include "certified_exponential.hpp"
#include "certified_sqrt.hpp"
#include "evaluation/evaluation_budget.hpp"

#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace mmcal::approximation {
namespace {

using mathematics::ConstantId;
using numeric::BigFloat;
using numeric::BigInt;
using numeric::DecimalApproximation;
using numeric::Rational;

struct AtanEnclosure final {
    RealInterval interval;
    std::size_t termsUsed = 0;
};

[[nodiscard]] RealInterval exactIntegerInterval(
    std::int64_t value,
    std::size_t precisionBits) {
    return RealInterval::fromRational(Rational{BigInt{value}}, precisionBits);
}



[[nodiscard]] AtanEnclosure encloseAtanReciprocal(
    std::uint32_t denominatorValue,
    std::size_t precisionBits) {
    if (denominatorValue <= 1)
        throw std::invalid_argument("Certified atan reciprocal denominator must exceed one");
    if (precisionBits == 0)
        throw std::invalid_argument("Certified atan precision must be at least one bit");

    // Machin公式で使う atan(1/q) は 0 < 1/q < 1 なので、
    //
    //   atan(1/q) = 1/q - 1/(3q^3) + 1/(5q^5) - ...
    //
    // は項の絶対値が単調減少する交代級数になる。
    // したがって真値は、隣り合う2つの部分和 S_k と S_(k+1) の間に必ず存在する。
    // この性質を使えば「項が小さくなったように見える」という経験則ではなく、
    // 厳密な包含区間として収束を判定できる。
    //
    // 各級数項そのものは Rational として定義し、それをBigFloatへ下向き/上向きに
    // 変換してからRealIntervalへ足す。つまり級数の打切り誤差だけでなく、
    // BigFloatの各演算丸めも区間の外向き丸めに吸収される。
    const BigInt q{static_cast<std::int64_t>(denominatorValue)};
    const BigInt qSquared = q * q;
    BigInt qPower = q;  // q^(2k+1), k=0ではq
    BigInt odd{1};      // 2k+1
    bool positive = true;

    const RealInterval zero = exactIntegerInterval(0, precisionBits);
    RealInterval sum = zero;
    std::size_t termsUsed = 0;

    // Pi全体では16倍/4倍するので、作業precisionより十分小さい次項まで進める。
    // 32bit余分に絞っておけば、Machin線形結合と途中の外向き丸めの余裕を大きく取れる。
    // 最終的な正しさはこの固定32bitに依存しない。decimal丸めが確定しなければ
    // approximateConstant側が作業precision自体を増やして再計算する。
    const std::size_t thresholdBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified atan precision is too large");
    const Rational threshold = binaryPrecisionThreshold(thresholdBits);

    for (;;) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        const Rational magnitude{BigInt{1}, odd * qPower};
        const RealInterval term = RealInterval::fromRational(magnitude, precisionBits);
        sum = positive
            ? add(sum, term, precisionBits)
            : subtract(sum, term, precisionBits);
        ++termsUsed;

        // 次の部分和を作る準備。次項を実際に現在値へ確定的に足す前に、
        // S_k と S_(k+1) の両方を包含するhullを作れば、交代級数定理により
        // atan(1/q)そのものも必ずそのhull内にある。
        qPower *= qSquared;
        odd += BigInt{2};
        const Rational nextMagnitude{BigInt{1}, odd * qPower};
        const RealInterval nextTerm = RealInterval::fromRational(
            nextMagnitude, precisionBits);
        const RealInterval nextSum = positive
            ? subtract(sum, nextTerm, precisionBits)
            : add(sum, nextTerm, precisionBits);

        if (nextMagnitude <= threshold)
            return AtanEnclosure{hull(sum, nextSum), termsUsed};

        positive = !positive;
    }
}

[[nodiscard]] std::optional<DecimalApproximation> tryCertifiedDecimal(
    const RealInterval& interval,
    std::size_t fractionalDigits) {
    return DecimalApproximation::fromCertifiedInterval(
        interval.lower().toRational(),
        interval.upper().toRational(),
        fractionalDigits);
}

[[nodiscard]] std::size_t nextGuardDigits(std::size_t guardDigits) {
    // 固定回数で諦めない。Piは有理数の10進丸め境界そのものにはならないので、
    // 区間を狭め続ければ最終的に要求桁の丸め結果は一意に決まる。
    // guardは少なくとも8桁ずつ、十分大きくなった後は約1.5倍で増やす。
    const std::size_t growth = std::max<std::size_t>(8, guardDigits / 2);
    return checkedPrecisionAdd(guardDigits, growth, "Certified constant precision is too large");
}

} // namespace

struct ChudnovskySplit final {
    BigInt p{1};
    BigInt q{1};
    BigInt t{};
};

[[nodiscard]] ChudnovskySplit chudnovskySplit(
    std::uint64_t begin,
    std::uint64_t end) {
    constexpr std::uint64_t a = 13'591'409ULL;
    constexpr std::uint64_t b = 545'140'134ULL;
    constexpr std::uint64_t c3Over24 = 10'939'058'860'032'000ULL;

    if (end != begin + 1) {
        const std::uint64_t middle = begin + (end - begin) / 2;
        ChudnovskySplit left = chudnovskySplit(begin, middle);
        ChudnovskySplit right = chudnovskySplit(middle, end);

        ChudnovskySplit result;
        result.p = left.p * right.p;
        result.q = left.q * right.q;
        result.t = left.t * right.q + left.p * right.t;
        return result;
    }

    if (begin == 0)
        return ChudnovskySplit{BigInt{1}, BigInt{1}, BigInt{static_cast<std::int64_t>(a)}};

    if (begin > std::numeric_limits<std::uint64_t>::max() / 6)
        throw std::overflow_error("Certified Pi Chudnovsky index is too large");

    const std::uint64_t sixK = 6 * begin;
    const BigInt k = BigInt::fromUnsigned(begin);
    ChudnovskySplit result;
    result.p = BigInt::fromUnsigned(sixK - 5)
        * BigInt::fromUnsigned(2 * begin - 1)
        * BigInt::fromUnsigned(sixK - 1);
    result.q = k * k * k * BigInt::fromUnsigned(c3Over24);

    // 旧実装では a + b * begin をuint64_tで先に計算していた。
    // 実用的な桁数では問題にならないが、極端なprecision指定では乗算がwrapし得る。
    // Chudnovsky係数もBigInt上で組み立て、任意精度という層の契約をここだけ64bitへ狭めない。
    const BigInt linear = BigInt::fromUnsigned(a) + BigInt::fromUnsigned(b) * k;
    result.t = result.p * linear;
    if ((begin & 1U) != 0)
        result.t = -result.t;
    return result;
}

[[nodiscard]] CertifiedConstantResult enclosePiChudnovsky(
    std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified Pi precision must be at least one bit");

    // Chudnovsky級数は1項あたり約47bit増える。46bit/項として保守的に多めに取り、
    // 正しさ自体は項数見積りへ依存させず、最初の未採用項による交代級数の上下界で保証する。
    const std::size_t estimatedTerms = precisionBits / 46 + 2;
    if (estimatedTerms >= std::numeric_limits<std::uint64_t>::max())
        throw std::overflow_error("Certified Pi term count is too large");
    const std::uint64_t terms = static_cast<std::uint64_t>(
        std::max<std::size_t>(2, estimatedTerms));
    evaluation::consumeEvaluationBudget(
        evaluation::EvaluationResource::CertifiedRefinement,
        static_cast<std::size_t>(terms));

    const ChudnovskySplit split = chudnovskySplit(0, terms);
    const ChudnovskySplit nextLeaf = chudnovskySplit(terms, terms + 1);

    const Rational partial{split.t, split.q};
    const Rational nextTerm{
        split.p * nextLeaf.t,
        split.q * nextLeaf.q};
    const Rational nextPartial = partial + nextTerm;
    const Rational lowerSum = nextPartial < partial ? nextPartial : partial;
    const Rational upperSum = nextPartial < partial ? partial : nextPartial;

    const std::size_t workBits = checkedPrecisionAdd(
        precisionBits, 32, "Certified Pi precision is too large");
    const RealInterval sum = RealInterval::fromRationalBounds(
        lowerSum, upperSum, workBits);
    const RealInterval root10005 = encloseSqrt(
        Rational{BigInt{10005}}, workBits).interval;
    const RealInterval scale = RealInterval::fromRational(
        Rational{BigInt{426880}}, workBits);
    const RealInterval numerator = multiply(root10005, scale, workBits);

    return CertifiedConstantResult{
        divide(numerator, sum, workBits).roundedOutward(precisionBits),
        static_cast<std::size_t>(terms),
        precisionBits};
}

[[maybe_unused]] CertifiedConstantResult enclosePiMachin(std::size_t precisionBits) {
    if (precisionBits == 0)
        throw std::invalid_argument("Certified Pi precision must be at least one bit");

    // Machinの公式:
    //
    //   pi = 16 atan(1/5) - 4 atan(1/239)
    //
    // 1/5, 1/239はいずれも絶対値が小さいためatan交代級数が速く収束する。
    // 初期Pi実装としてChudnovskyより遅いが、
    //   * 証明が短い
    //   * 交代級数だけで厳密な上下界が作れる
    //   * BigFloat/RealIntervalの検証用referenceとして残せる
    // という利点を優先する。
    const AtanEnclosure atan5 = encloseAtanReciprocal(5, precisionBits);
    const AtanEnclosure atan239 = encloseAtanReciprocal(239, precisionBits);

    const RealInterval sixteen = exactIntegerInterval(16, precisionBits);
    const RealInterval four = exactIntegerInterval(4, precisionBits);
    const RealInterval first = multiply(atan5.interval, sixteen, precisionBits);
    const RealInterval second = multiply(atan239.interval, four, precisionBits);

    return CertifiedConstantResult{
        subtract(first, second, precisionBits),
        atan5.termsUsed + atan239.termsUsed,
        precisionBits
    };
}

CertifiedConstantResult enclosePi(std::size_t precisionBits) {
    /*
    旧実装:
        Pi = 16 atan(1/5) - 4 atan(1/239)
        をRealIntervalの交代Taylor級数で直接評価していた。

    証明が短くreferenceとして有用だが、数千～1万桁では各項の多倍長除算が支配的になる。
    既定経路はbinary-splitting Chudnovskyへ変更し、旧Machin実装は
    MMCAL_USE_MACHIN_PI指定時の比較/referenceとしてこのファイル内に残す。
    */
#ifdef MMCAL_USE_MACHIN_PI
    return enclosePiMachin(precisionBits);
#else
    return enclosePiChudnovsky(precisionBits);
#endif
}


std::optional<CertifiedConstantResult> encloseConstant(
    ConstantId id,
    std::size_t precisionBits) {
    switch (id) {
    case ConstantId::Pi:
        return enclosePi(precisionBits);
    case ConstantId::E: {
        // Eは「保存された10進定数」ではなく exp(1) という数学的定義から生成する。
        // Exp backendが返す区間はTaylor剰余とBigFloat丸めを両方包含している。
        const RealInterval one = RealInterval::fromRational(Rational{BigInt{1}}, precisionBits);
        const CertifiedExponentialResult exponential = encloseExp(one, precisionBits);
        return CertifiedConstantResult{
            exponential.interval,
            exponential.termsUsed,
            precisionBits
        };
    }
    case ConstantId::Phi: {
        // Phi = (1 + sqrt(5)) / 2。sqrt(5)のcertified enclosureから
        // 四則演算を外向き丸めで組み立て、保存小数には依存しない。
        const RealInterval one = exactIntegerInterval(1, precisionBits);
        const RealInterval half = RealInterval::fromRational(
            Rational{BigInt{1}, BigInt{2}}, precisionBits);
        const RealInterval sqrtFive = encloseSqrt(
            Rational{BigInt{5}}, precisionBits).interval;
        return CertifiedConstantResult{
            multiply(add(one, sqrtFive, precisionBits), half, precisionBits),
            0,
            precisionBits
        };
    }
    }

    return std::nullopt;
}

DecimalApproximation approximatePi(std::size_t fractionalDigits) {
    const auto result = approximateConstant(ConstantId::Pi, fractionalDigits);
    if (!result)
        throw std::logic_error("Certified Pi provider is not available");
    return *result;
}

std::optional<DecimalApproximation> approximateConstant(
    ConstantId id,
    std::size_t fractionalDigits) {
    if (fractionalDigits == 0)
        throw std::invalid_argument("Certified constant precision must be greater than zero");

    ApproximationContext context{fractionalDigits};
    for (;;) {
        evaluation::consumeEvaluationBudget(
            evaluation::EvaluationResource::CertifiedRefinement);
        const auto enclosure = encloseConstant(id, context.workingBinaryBits());
        if (!enclosure)
            return std::nullopt;

        if (const auto decimal = tryCertifiedDecimal(enclosure->interval, fractionalDigits))
            return decimal;

        context.setGuardDigits(nextGuardDigits(context.guardDigits()));
    }
}

} // namespace mmcal::approximation
