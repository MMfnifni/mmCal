// BigInt内部の符号なし多倍長整数
#include "big_uint.hpp"

#include <algorithm>
#include <bit>
#include <limits>
#include <span>
#include <stdexcept>
#include <utility>

namespace mmcal::numeric::detail {
namespace {

constexpr std::uint64_t limbMask = std::numeric_limits<std::uint32_t>::max();
constexpr std::uint64_t limbBase = std::uint64_t{1} << 32;
constexpr unsigned limbBits = 32;

// 実測でschoolbookとのcrossoverを決める。benchmark時だけ-Dで上書きできる。
#ifdef MMCAL_KARATSUBA_THRESHOLD_LIMBS
constexpr std::size_t karatsubaThresholdLimbs = MMCAL_KARATSUBA_THRESHOLD_LIMBS;
#else
constexpr std::size_t karatsubaThresholdLimbs = 48;
#endif
#ifdef MMCAL_TOOM3_THRESHOLD_LIMBS
constexpr std::size_t toom3ThresholdLimbs = MMCAL_TOOM3_THRESHOLD_LIMBS;
#else
constexpr std::size_t toom3ThresholdLimbs = 1280;
#endif
#ifdef MMCAL_TOOM3_RECURSIVE_THRESHOLD_LIMBS
constexpr std::size_t toom3RecursiveThresholdLimbs = MMCAL_TOOM3_RECURSIVE_THRESHOLD_LIMBS;
#else
constexpr std::size_t toom3RecursiveThresholdLimbs = 448;
#endif
#ifdef MMCAL_SQUARE_KARATSUBA_THRESHOLD_LIMBS
constexpr std::size_t squareKaratsubaThresholdLimbs = MMCAL_SQUARE_KARATSUBA_THRESHOLD_LIMBS;
#else
constexpr std::size_t squareKaratsubaThresholdLimbs = 48;
#endif

static_assert(karatsubaThresholdLimbs >= 1, "Karatsuba threshold must be at least one limb");
static_assert(squareKaratsubaThresholdLimbs >= 1, "Square Karatsuba threshold must be at least one limb");
static_assert(toom3ThresholdLimbs > karatsubaThresholdLimbs,
    "Toom-3 threshold must be greater than the Karatsuba threshold");
static_assert(toom3RecursiveThresholdLimbs > karatsubaThresholdLimbs,
    "Recursive Toom-3 threshold must be greater than the Karatsuba threshold");

constexpr BigUInt::limb_type decimalChunkBase = 1'000'000'000u;
constexpr unsigned decimalChunkDigits = 9;

// 10進文字列化は巨大値だけdivide-and-conquerへ切り替える。Knuth除算の固定費があるため、
// crossover未満では10^9 chunkの逐次除算を維持する。benchmark時だけ-Dで上書きできる。
#ifdef MMCAL_DECIMAL_DAC_THRESHOLD_LIMBS
constexpr std::size_t decimalDacThresholdLimbs = MMCAL_DECIMAL_DAC_THRESHOLD_LIMBS;
#else
constexpr std::size_t decimalDacThresholdLimbs = 128;
#endif
#ifdef MMCAL_DECIMAL_DAC_LEAF_LIMBS
constexpr std::size_t decimalDacLeafLimbs = MMCAL_DECIMAL_DAC_LEAF_LIMBS;
#else
constexpr std::size_t decimalDacLeafLimbs = 256;
#endif

using Limb = BigUInt::limb_type;
using DoubleLimb = BigUInt::double_limb_type;

void normalizeLimbs(std::vector<Limb>& limbs) noexcept {
    while (!limbs.empty() && limbs.back() == 0)
        limbs.pop_back();
}

[[nodiscard]] std::vector<Limb> addLimbs(
    std::span<const Limb> lhs,
    std::span<const Limb> rhs) {
    const std::size_t size = std::max(lhs.size(), rhs.size());
    std::vector<Limb> result(size, 0);
    DoubleLimb carry = 0;

    for (std::size_t index = 0; index < size; ++index) {
        const DoubleLimb lhsLimb = index < lhs.size() ? lhs[index] : 0;
        const DoubleLimb rhsLimb = index < rhs.size() ? rhs[index] : 0;
        const DoubleLimb sum = lhsLimb + rhsLimb + carry;
        result[index] = static_cast<Limb>(sum & limbMask);
        carry = sum >> limbBits;
    }

    if (carry != 0)
        result.push_back(static_cast<Limb>(carry));
    return result;
}

void subtractLimbsInPlace(std::vector<Limb>& lhs, std::span<const Limb> rhs) {
    DoubleLimb borrow = 0;

    for (std::size_t index = 0; index < lhs.size(); ++index) {
        const DoubleLimb lhsLimb = lhs[index];
        const DoubleLimb rhsLimb = index < rhs.size() ? rhs[index] : 0;
        const DoubleLimb subtrahend = rhsLimb + borrow;

        if (lhsLimb >= subtrahend) {
            lhs[index] = static_cast<Limb>(lhsLimb - subtrahend);
            borrow = 0;
        }
        else {
            lhs[index] = static_cast<Limb>(limbBase + lhsLimb - subtrahend);
            borrow = 1;
        }
    }

    if (borrow != 0)
        throw std::logic_error("BigUInt Karatsuba subtraction underflow");
    normalizeLimbs(lhs);
}

void addShiftedLimbs(
    std::vector<Limb>& destination,
    std::span<const Limb> source,
    std::size_t offset) {
    if (source.empty())
        return;

    if (offset > destination.size() || source.size() > destination.size() - offset)
        throw std::logic_error("BigUInt Karatsuba shifted sum exceeds result size");

    DoubleLimb carry = 0;
    std::size_t index = 0;
    for (; index < source.size(); ++index) {
        const std::size_t target = offset + index;
        const DoubleLimb sum =
            static_cast<DoubleLimb>(destination[target]) + source[index] + carry;
        destination[target] = static_cast<Limb>(sum & limbMask);
        carry = sum >> limbBits;
    }

    std::size_t target = offset + index;
    while (carry != 0) {
        if (target == destination.size())
            throw std::logic_error("BigUInt Karatsuba carry exceeds result size");
        const DoubleLimb sum = static_cast<DoubleLimb>(destination[target]) + carry;
        destination[target] = static_cast<Limb>(sum & limbMask);
        carry = sum >> limbBits;
        ++target;
    }
}

[[nodiscard]] std::vector<Limb> multiplySchoolbook(
    std::span<const Limb> lhs,
    std::span<const Limb> rhs) {
    if (lhs.empty() || rhs.empty())
        return {};

    std::vector<Limb> result(lhs.size() + rhs.size(), 0);

    for (std::size_t lhsIndex = 0; lhsIndex < lhs.size(); ++lhsIndex) {
        DoubleLimb carry = 0;

        for (std::size_t rhsIndex = 0; rhsIndex < rhs.size(); ++rhsIndex) {
            const std::size_t resultIndex = lhsIndex + rhsIndex;
            const DoubleLimb product =
                static_cast<DoubleLimb>(lhs[lhsIndex]) * rhs[rhsIndex]
                + result[resultIndex]
                + carry;

            result[resultIndex] = static_cast<Limb>(product & limbMask);
            carry = product >> limbBits;
        }

        result[lhsIndex + rhs.size()] = static_cast<Limb>(carry);
    }

    normalizeLimbs(result);
    return result;
}

[[nodiscard]] std::vector<Limb> multiplyAdaptive(
    std::span<const Limb> lhs,
    std::span<const Limb> rhs,
    bool insideToom);

void addWideAt(std::vector<Limb>& destination, std::size_t offset, DoubleLimb value) {
    while (value != 0) {
        if (offset >= destination.size())
            throw std::logic_error("BigUInt square carry exceeds result size");

        const DoubleLimb sum =
            static_cast<DoubleLimb>(destination[offset]) + (value & limbMask);
        destination[offset] = static_cast<Limb>(sum & limbMask);
        value = (value >> limbBits) + (sum >> limbBits);
        ++offset;
    }
}

// a^2では非対角項 a[i]a[j] が必ず2回現れる。一般schoolbookのn^2回の
// limb乗算をそのまま行わず、対角n項と上三角n(n-1)/2項だけを計算する。
[[nodiscard]] std::vector<Limb> squareSchoolbook(std::span<const Limb> value) {
    if (value.empty())
        return {};

    std::vector<Limb> result(value.size() * 2, 0);
    for (std::size_t i = 0; i < value.size(); ++i) {
        const DoubleLimb diagonal = static_cast<DoubleLimb>(value[i]) * value[i];
        addWideAt(result, i * 2, diagonal);

        for (std::size_t j = i + 1; j < value.size(); ++j) {
            const DoubleLimb cross = static_cast<DoubleLimb>(value[i]) * value[j];
            // 2*crossは65bitになり得るため、同じ64bit積を2回加えてportableに処理する。
            addWideAt(result, i + j, cross);
            addWideAt(result, i + j, cross);
        }
    }

    normalizeLimbs(result);
    return result;
}

[[nodiscard]] std::vector<Limb> squareAdaptive(
    std::span<const Limb> value,
    bool insideToom);

[[nodiscard]] std::vector<Limb> squareKaratsuba(
    std::span<const Limb> value,
    bool insideToom) {
    if (value.size() <= squareKaratsubaThresholdLimbs)
        return squareSchoolbook(value);

    const std::size_t split = value.size() / 2;
    const auto low = value.first(split);
    const auto high = value.subspan(split);

    auto z0 = squareAdaptive(low, insideToom);
    auto z2 = squareAdaptive(high, insideToom);
    const auto sum = addLimbs(low, high);
    auto z1 = squareAdaptive(sum, insideToom);
    subtractLimbsInPlace(z1, z0);
    subtractLimbsInPlace(z1, z2);

    std::vector<Limb> result(value.size() * 2, 0);
    addShiftedLimbs(result, z0, 0);
    addShiftedLimbs(result, z1, split);
    addShiftedLimbs(result, z2, split * 2);
    normalizeLimbs(result);
    return result;
}

[[nodiscard]] std::vector<Limb> multiplyKaratsuba(
    std::span<const Limb> lhs,
    std::span<const Limb> rhs,
    bool insideToom) {
    if (lhs.empty() || rhs.empty())
        return {};

    const std::size_t smallerSize = std::min(lhs.size(), rhs.size());
    const std::size_t largerSize = std::max(lhs.size(), rhs.size());

    // 小さい積と極端に不均衡な積は、Karatsubaの一時配列・再帰コストが勝る。
    if (smallerSize <= karatsubaThresholdLimbs || largerSize - smallerSize > smallerSize)
        return multiplySchoolbook(lhs, rhs);

    const std::size_t split = largerSize / 2;
    const std::size_t lhsLowSize = std::min(lhs.size(), split);
    const std::size_t rhsLowSize = std::min(rhs.size(), split);

    const auto lhsLow = lhs.first(lhsLowSize);
    const auto lhsHigh = lhs.subspan(lhsLowSize);
    const auto rhsLow = rhs.first(rhsLowSize);
    const auto rhsHigh = rhs.subspan(rhsLowSize);

    auto z0 = multiplyAdaptive(lhsLow, rhsLow, insideToom);
    auto z2 = multiplyAdaptive(lhsHigh, rhsHigh, insideToom);
    const auto lhsSum = addLimbs(lhsLow, lhsHigh);
    const auto rhsSum = addLimbs(rhsLow, rhsHigh);
    auto z1 = multiplyAdaptive(lhsSum, rhsSum, insideToom);
    subtractLimbsInPlace(z1, z0);
    subtractLimbsInPlace(z1, z2);

    std::vector<Limb> result(lhs.size() + rhs.size(), 0);
    addShiftedLimbs(result, z0, 0);
    addShiftedLimbs(result, z1, split);
    addShiftedLimbs(result, z2, split * 2);
    normalizeLimbs(result);
    return result;
}


[[nodiscard]] int compareLimbs(std::span<const Limb> lhs, std::span<const Limb> rhs) noexcept {
    while (!lhs.empty() && lhs.back() == 0)
        lhs = lhs.first(lhs.size() - 1);
    while (!rhs.empty() && rhs.back() == 0)
        rhs = rhs.first(rhs.size() - 1);
    if (lhs.size() != rhs.size())
        return lhs.size() < rhs.size() ? -1 : 1;
    for (std::size_t i = lhs.size(); i-- > 0;) {
        if (lhs[i] != rhs[i])
            return lhs[i] < rhs[i] ? -1 : 1;
    }
    return 0;
}

[[nodiscard]] std::vector<Limb> subtractLimbs(
    std::span<const Limb> lhs,
    std::span<const Limb> rhs) {
    if (compareLimbs(lhs, rhs) < 0)
        throw std::logic_error("BigUInt Toom-3 subtraction underflow");
    std::vector<Limb> result(lhs.begin(), lhs.end());
    subtractLimbsInPlace(result, rhs);
    return result;
}

[[nodiscard]] std::vector<Limb> multiplyLimbsSmall(std::span<const Limb> value, Limb factor) {
    if (value.empty() || factor == 0)
        return {};
    std::vector<Limb> result(value.size(), 0);
    DoubleLimb carry = 0;
    for (std::size_t i = 0; i < value.size(); ++i) {
        const DoubleLimb product = static_cast<DoubleLimb>(value[i]) * factor + carry;
        result[i] = static_cast<Limb>(product & limbMask);
        carry = product >> limbBits;
    }
    if (carry != 0)
        result.push_back(static_cast<Limb>(carry));
    return result;
}

void divideLimbsSmallExact(std::vector<Limb>& value, Limb divisor) {
    DoubleLimb remainder = 0;
    for (std::size_t i = value.size(); i-- > 0;) {
        const DoubleLimb current = (remainder << limbBits) | value[i];
        value[i] = static_cast<Limb>(current / divisor);
        remainder = current % divisor;
    }
    if (remainder != 0)
        throw std::logic_error("BigUInt Toom-3 interpolation division was not exact");
    normalizeLimbs(value);
}

struct SignedLimbs final {
    bool negative = false;
    std::vector<Limb> magnitude;
};

[[nodiscard]] SignedLimbs signedPositive(std::span<const Limb> value) {
    return {false, std::vector<Limb>{value.begin(), value.end()}};
}

[[nodiscard]] SignedLimbs signedNegate(SignedLimbs value) {
    if (!value.magnitude.empty())
        value.negative = !value.negative;
    return value;
}

[[nodiscard]] SignedLimbs signedAdd(const SignedLimbs& lhs, const SignedLimbs& rhs) {
    if (lhs.negative == rhs.negative)
        return {lhs.negative, addLimbs(lhs.magnitude, rhs.magnitude)};
    const int order = compareLimbs(lhs.magnitude, rhs.magnitude);
    if (order == 0)
        return {};
    if (order > 0)
        return {lhs.negative, subtractLimbs(lhs.magnitude, rhs.magnitude)};
    return {rhs.negative, subtractLimbs(rhs.magnitude, lhs.magnitude)};
}

[[nodiscard]] SignedLimbs signedSub(const SignedLimbs& lhs, const SignedLimbs& rhs) {
    return signedAdd(lhs, signedNegate(rhs));
}

[[nodiscard]] SignedLimbs signedDivideExact(SignedLimbs value, Limb divisor) {
    divideLimbsSmallExact(value.magnitude, divisor);
    if (value.magnitude.empty())
        value.negative = false;
    return value;
}

[[nodiscard]] std::span<const Limb> chunkAt(
    std::span<const Limb> value,
    std::size_t offset,
    std::size_t size) noexcept {
    if (offset >= value.size())
        return {};
    return value.subspan(offset, std::min(size, value.size() - offset));
}

[[nodiscard]] std::vector<Limb> evaluateAtOne(
    std::span<const Limb> x0,
    std::span<const Limb> x1,
    std::span<const Limb> x2) {
    auto result = addLimbs(x0, x1);
    return addLimbs(result, x2);
}

[[nodiscard]] SignedLimbs evaluateAtMinusOne(
    std::span<const Limb> x0,
    std::span<const Limb> x1,
    std::span<const Limb> x2) {
    return signedSub({false, addLimbs(x0, x2)}, signedPositive(x1));
}

[[nodiscard]] std::vector<Limb> evaluateAtTwo(
    std::span<const Limb> x0,
    std::span<const Limb> x1,
    std::span<const Limb> x2) {
    auto result = addLimbs(x0, multiplyLimbsSmall(x1, 2));
    return addLimbs(result, multiplyLimbsSmall(x2, 4));
}

[[nodiscard]] SignedLimbs multiplySigned(
    const SignedLimbs& lhs,
    const SignedLimbs& rhs) {
    SignedLimbs result;
    result.negative = lhs.negative != rhs.negative;
    result.magnitude = multiplyAdaptive(lhs.magnitude, rhs.magnitude, true);
    if (result.magnitude.empty())
        result.negative = false;
    return result;
}

void requireNonNegative(const SignedLimbs& value, const char* message) {
    if (value.negative)
        throw std::logic_error(message);
}

// 3分割したoperandを t={0,1,-1,2,∞} で評価し、5回の再帰乗算から補間する。
// Karatsubaより定数項は重いが、十分巨大なbalanced operandでは乗算回数の減少が勝る。
[[nodiscard]] std::vector<Limb> multiplyToom3(
    std::span<const Limb> lhs,
    std::span<const Limb> rhs) {
    const std::size_t split = (std::max(lhs.size(), rhs.size()) + 2) / 3;
    const auto a0 = chunkAt(lhs, 0, split);
    const auto a1 = chunkAt(lhs, split, split);
    const auto a2 = chunkAt(lhs, split * 2, split);
    const auto b0 = chunkAt(rhs, 0, split);
    const auto b1 = chunkAt(rhs, split, split);
    const auto b2 = chunkAt(rhs, split * 2, split);

    const auto v0 = multiplyAdaptive(a0, b0, true);
    const auto v4 = multiplyAdaptive(a2, b2, true);
    const SignedLimbs v1{false, multiplyAdaptive(
        evaluateAtOne(a0, a1, a2), evaluateAtOne(b0, b1, b2), true)};
    const SignedLimbs vm1 = multiplySigned(
        evaluateAtMinusOne(a0, a1, a2), evaluateAtMinusOne(b0, b1, b2));
    const SignedLimbs v2{false, multiplyAdaptive(
        evaluateAtTwo(a0, a1, a2), evaluateAtTwo(b0, b1, b2), true)};

    const SignedLimbs sum13 = signedDivideExact(signedSub(v1, vm1), 2);
    requireNonNegative(sum13, "BigUInt Toom-3 c1+c3 became negative");

    SignedLimbs c2 = signedDivideExact(signedAdd(v1, vm1), 2);
    c2 = signedSub(c2, signedPositive(v0));
    c2 = signedSub(c2, signedPositive(v4));
    requireNonNegative(c2, "BigUInt Toom-3 c2 became negative");

    SignedLimbs t = signedSub(v2, signedPositive(v0));
    t = signedSub(t, {false, multiplyLimbsSmall(c2.magnitude, 4)});
    t = signedSub(t, {false, multiplyLimbsSmall(v4, 16)});
    t = signedDivideExact(std::move(t), 2);
    requireNonNegative(t, "BigUInt Toom-3 interpolation term became negative");

    SignedLimbs c3 = signedDivideExact(signedSub(t, sum13), 3);
    requireNonNegative(c3, "BigUInt Toom-3 c3 became negative");
    SignedLimbs c1 = signedSub(sum13, c3);
    requireNonNegative(c1, "BigUInt Toom-3 c1 became negative");

    std::vector<Limb> result(lhs.size() + rhs.size(), 0);
    addShiftedLimbs(result, v0, 0);
    addShiftedLimbs(result, c1.magnitude, split);
    addShiftedLimbs(result, c2.magnitude, split * 2);
    addShiftedLimbs(result, c3.magnitude, split * 3);
    addShiftedLimbs(result, v4, split * 4);
    normalizeLimbs(result);
    return result;
}

[[nodiscard]] std::vector<Limb> squareAdaptive(
    std::span<const Limb> value,
    bool insideToom) {
    if (value.empty())
        return {};
    if (value.size() <= squareKaratsubaThresholdLimbs)
        return squareSchoolbook(value);

    /*
    初版では一般乗算と同じthresholdで専用Toom-3 squareへ切り替えたが、
    1536～4096 limbsの実測でKaratsuba squareより大幅に遅かった。
    評価点・補間の固定費に対し、squareではKaratsuba自身の対称性が十分強いため、
    現段階では巨大squareもKaratsuba再帰を継続する。Toom squareは再設計時に再評価する。

    if (value.size() >= activeToomThreshold)
        return squareToom3(value);
    */
    return squareKaratsuba(value, insideToom);
}

[[nodiscard]] std::vector<Limb> multiplyAdaptive(
    std::span<const Limb> lhs,
    std::span<const Limb> rhs,
    bool insideToom) {
    if (lhs.empty() || rhs.empty())
        return {};
    const std::size_t smallerSize = std::min(lhs.size(), rhs.size());
    const std::size_t largerSize = std::max(lhs.size(), rhs.size());
    if (smallerSize <= karatsubaThresholdLimbs || largerSize - smallerSize > smallerSize)
        return multiplySchoolbook(lhs, rhs);

    const std::size_t activeToomThreshold = insideToom
        ? toom3RecursiveThresholdLimbs
        : toom3ThresholdLimbs;
    if (smallerSize >= activeToomThreshold)
        return multiplyToom3(lhs, rhs);
    return multiplyKaratsuba(lhs, rhs, insideToom);
}

void validateRadix(unsigned radix) {
    if (radix < 2 || radix > 36)
        throw std::invalid_argument("BigUInt radix must be in the range 2..36");
}

// 被除数の対象limb区間から quotientDigit * divisor を減算する。
// 商の推定値が1大きすぎた場合は true を返す。
bool subtractProduct(
    std::vector<BigUInt::limb_type>& dividend,
    std::size_t offset,
    const std::vector<BigUInt::limb_type>& divisor,
    std::uint64_t quotientDigit) noexcept {
    std::uint64_t carry = 0;
    std::uint64_t borrow = 0;

    for (std::size_t index = 0; index < divisor.size(); ++index) {
        const std::uint64_t product =
            quotientDigit * divisor[index] + carry;
        carry = product >> limbBits;

        const std::uint64_t subtrahend = (product & limbMask) + borrow;
        const std::uint64_t current = dividend[offset + index];

        if (current < subtrahend) {
            dividend[offset + index] = static_cast<BigUInt::limb_type>(
                limbBase + current - subtrahend);
            borrow = 1;
        }
        else {
            dividend[offset + index] = static_cast<BigUInt::limb_type>(
                current - subtrahend);
            borrow = 0;
        }
    }

    const std::uint64_t highSubtrahend = carry + borrow;
    const std::uint64_t high = dividend[offset + divisor.size()];

    if (high < highSubtrahend) {
        dividend[offset + divisor.size()] = static_cast<BigUInt::limb_type>(
            limbBase + high - highSubtrahend);
        return true;
    }

    dividend[offset + divisor.size()] = static_cast<BigUInt::limb_type>(
        high - highSubtrahend);
    return false;
}

void addBack(
    std::vector<BigUInt::limb_type>& dividend,
    std::size_t offset,
    const std::vector<BigUInt::limb_type>& divisor) noexcept {
    std::uint64_t carry = 0;

    for (std::size_t index = 0; index < divisor.size(); ++index) {
        const std::uint64_t sum =
            static_cast<std::uint64_t>(dividend[offset + index])
            + divisor[index]
            + carry;

        dividend[offset + index] = static_cast<BigUInt::limb_type>(sum & limbMask);
        carry = sum >> limbBits;
    }

    dividend[offset + divisor.size()] = static_cast<BigUInt::limb_type>(
        static_cast<std::uint64_t>(dividend[offset + divisor.size()]) + carry);
}

} // namespace

BigUInt::BigUInt(std::uint64_t value) {
    if (value == 0)
        return;

    limbs_.push_back(static_cast<limb_type>(value & limbMask));

    const auto high = static_cast<limb_type>(value >> limbBits);
    if (high != 0)
        limbs_.push_back(high);
}

bool BigUInt::isZero() const noexcept {
    return limbs_.empty();
}

std::size_t BigUInt::limbCount() const noexcept {
    return limbs_.size();
}

std::size_t BigUInt::bitLength() const noexcept {
    if (isZero())
        return 0;

    const auto highest = limbs_.back();
    return (limbs_.size() - 1) * limbBits + std::bit_width(highest);
}

std::size_t BigUInt::trailingZeroBits() const noexcept {
    if (isZero())
        return 0;

    std::size_t bits = 0;
    for (const auto limb : limbs_) {
        if (limb == 0) {
            bits += limbBits;
            continue;
        }

        bits += std::countr_zero(limb);
        break;
    }
    return bits;
}

BigUInt BigUInt::parse(std::string_view text, unsigned radix) {
    validateRadix(radix);

    if (text.empty())
        throw std::invalid_argument("BigUInt cannot parse an empty string");

    std::size_t position = 0;
    if (text.front() == '+')
        position = 1;
    else if (text.front() == '-')
        throw std::invalid_argument("BigUInt cannot parse a negative value");

    if (position == text.size())
        throw std::invalid_argument("BigUInt requires at least one digit");

    BigUInt result;

    if (radix == 10) {
        /*
        旧実装は10進入力を1桁ずつ multiplySmall(10) していた。巨大整数では全limb走査が
        桁数回発生するため、9桁chunkで multiplySmall(10^9) し、走査回数を最大1/9へ減らす。

        for (; position < text.size(); ++position) {
            const unsigned digit = digitValue(text[position]);
            if (digit >= radix)
                throw std::invalid_argument(
                    "BigUInt contains a digit invalid for the selected radix");
            result.multiplySmall(static_cast<limb_type>(radix));
            result.addSmall(static_cast<limb_type>(digit));
        }
        */
        const std::size_t digits = text.size() - position;
        std::size_t chunkDigits = digits % decimalChunkDigits;
        if (chunkDigits == 0)
            chunkDigits = decimalChunkDigits;
        while (position < text.size()) {
            limb_type chunk = 0;
            limb_type power = 1;
            for (std::size_t i = 0; i < chunkDigits; ++i) {
                const unsigned digit = digitValue(text[position++]);
                if (digit >= 10)
                    throw std::invalid_argument(
                        "BigUInt contains a digit invalid for the selected radix");
                chunk = static_cast<limb_type>(chunk * 10 + digit);
                power *= 10;
            }
            result.multiplySmall(chunkDigits == decimalChunkDigits ? decimalChunkBase : power);
            result.addSmall(chunk);
            chunkDigits = decimalChunkDigits;
        }
        return result;
    }

    for (; position < text.size(); ++position) {
        const unsigned digit = digitValue(text[position]);
        if (digit >= radix)
            throw std::invalid_argument(
                "BigUInt contains a digit invalid for the selected radix");
        result.multiplySmall(static_cast<limb_type>(radix));
        result.addSmall(static_cast<limb_type>(digit));
    }
    return result;
}

std::string BigUInt::toString(unsigned radix) const {
    validateRadix(radix);

    if (isZero())
        return "0";

    BigUInt remaining = *this;
    if (radix == 10) {
        /*
        旧実装は10で1桁ずつdivideSmallしていたため、Karatsuba導入後の巨大factorialでは
        計算本体より10進変換が圧倒的に重かった。10^9で9桁ずつ取り出して除算回数を減らす。

        std::string result;
        result.reserve((bitLength() + 2) / 3);
        while (!remaining.isZero()) {
            const auto remainder = remaining.divideSmall(static_cast<limb_type>(radix));
            result.push_back(digitCharacter(remainder));
        }
        std::reverse(result.begin(), result.end());
        return result;
        */
        const auto chunkedConversion = [&](BigUInt value) {
            if (value.isZero())
                return std::string{"0"};

            std::vector<limb_type> chunks;
            chunks.reserve((value.bitLength() + 28) / 29);
            while (!value.isZero())
                chunks.push_back(value.divideSmall(decimalChunkBase));

            std::string result = std::to_string(chunks.back());
            result.reserve(chunks.size() * decimalChunkDigits);
            for (std::size_t i = chunks.size() - 1; i-- > 0;) {
                const std::string chunk = std::to_string(chunks[i]);
                result.append(decimalChunkDigits - chunk.size(), '0');
                result += chunk;
            }
            return result;
        };

        if (limbs_.size() < decimalDacThresholdLimbs)
            return chunkedConversion(std::move(remaining));

        /*
        10^9 chunk化後も巨大値では「全limbを10^9で割る」処理をchunk数だけ繰り返すため
        O(n^2)的な走査が残っていた。初版divide-and-conquerは「値以下で最大の10^(9*2^k)」
        を常にsplitに使ったため、10進chunk数が2^kを少し超えた値で上側が極端に小さくなり、
        サイズ境界ごとに性能の谷ができた。

        現在は概算10進chunk数の半分以下で最大の2^k chunkをsplit幅に選び、上下をほぼ
        balancedにする。また小さくなった再帰葉ではKnuth除算を続けず10^9逐次変換へ戻す。

        std::vector<limb_type> chunks;
        chunks.reserve((bitLength() + 28) / 29);
        while (!remaining.isZero())
            chunks.push_back(remaining.divideSmall(decimalChunkBase));
        ...
        */
        const auto decimalChunkEstimate = [](const BigUInt& value) -> std::size_t {
            const std::size_t bits = value.bitLength();
            // log10(2) < 30103/100000。積のoverflowを避けて10進桁数の安全な上界を作る。
            const std::size_t digits =
                (bits / 100000) * 30103
                + ((bits % 100000) * 30103) / 100000
                + 1;
            return (digits + decimalChunkDigits - 1) / decimalChunkDigits;
        };

        const std::size_t estimatedChunks = decimalChunkEstimate(*this);
        const std::size_t maxLevel = std::bit_width(std::max<std::size_t>(estimatedChunks / 2, 1)) - 1;
        std::vector<BigUInt> powers;
        powers.reserve(maxLevel + 1);
        powers.emplace_back(decimalChunkBase);
        for (std::size_t level = 1; level <= maxLevel; ++level) {
            BigUInt next = powers.back();
            next *= next;
            powers.push_back(std::move(next));
        }

        const auto fixedWidth = [](std::size_t level) -> std::size_t {
            if (level >= std::numeric_limits<std::size_t>::digits)
                throw std::length_error("BigUInt decimal conversion width is too large");
            const std::size_t chunks = std::size_t{1} << level;
            if (chunks > std::numeric_limits<std::size_t>::max() / decimalChunkDigits)
                throw std::length_error("BigUInt decimal conversion width is too large");
            return chunks * decimalChunkDigits;
        };

        const auto convertFixed = [&](const auto& self, const BigUInt& value, std::size_t level) -> std::string {
            const std::size_t width = fixedWidth(level);
            if (level == 0 || value.limbCount() < decimalDacLeafLimbs) {
                std::string text = chunkedConversion(value);
                if (text.size() > width)
                    throw std::logic_error("BigUInt decimal fixed-width conversion overflow");
                text.insert(text.begin(), width - text.size(), '0');
                return text;
            }

            auto parts = divmod(value, powers[level - 1]);
            std::string high = self(self, parts.quotient, level - 1);
            std::string low = self(self, parts.remainder, level - 1);
            high += low;
            return high;
        };

        const auto convertVariable = [&](const auto& self, const BigUInt& value) -> std::string {
            if (value.isZero() || value.limbCount() < decimalDacThresholdLimbs)
                return chunkedConversion(value);

            const std::size_t chunks = decimalChunkEstimate(value);
            if (chunks <= 1)
                return chunkedConversion(value);
            const std::size_t halfChunks = std::max<std::size_t>(chunks / 2, 1);
            std::size_t level = std::bit_width(halfChunks) - 1;
            level = std::min(level, powers.size() - 1);

            auto parts = divmod(value, powers[level]);
            while (parts.quotient.isZero() && level != 0) {
                --level;
                parts = divmod(value, powers[level]);
            }

            std::string high = self(self, parts.quotient);
            std::string low = convertFixed(convertFixed, parts.remainder, level);
            high.reserve(high.size() + fixedWidth(level));
            high += low;
            return high;
        };

        return convertVariable(convertVariable, *this);
    }

    std::string result;
    result.reserve((bitLength() + 2) / 3);
    while (!remaining.isZero()) {
        const auto remainder = remaining.divideSmall(static_cast<limb_type>(radix));
        result.push_back(digitCharacter(remainder));
    }
    std::reverse(result.begin(), result.end());
    return result;
}

BigUInt& BigUInt::operator+=(const BigUInt& rhs) {
    // 自己加算ではresizeによりrhs参照が無効化され得るため、先に退避する。
    if (this == &rhs) {
        const BigUInt copy = rhs;
        return *this += copy;
    }

    const std::size_t requiredSize =
        std::max(limbs_.size(), rhs.limbs_.size());
    limbs_.resize(requiredSize, 0);

    double_limb_type carry = 0;

    for (std::size_t index = 0; index < requiredSize; ++index) {
        const double_limb_type rhsLimb =
            index < rhs.limbs_.size() ? rhs.limbs_[index] : 0;

        const double_limb_type sum =
            static_cast<double_limb_type>(limbs_[index]) + rhsLimb + carry;

        limbs_[index] = static_cast<limb_type>(sum & limbMask);
        carry = sum >> limbBits;
    }

    if (carry != 0)
        limbs_.push_back(static_cast<limb_type>(carry));

    return *this;
}

BigUInt& BigUInt::operator-=(const BigUInt& rhs) {
    if (*this < rhs)
        throw std::underflow_error("BigUInt subtraction would produce a negative value");

    if (this == &rhs) {
        limbs_.clear();
        return *this;
    }

    double_limb_type borrow = 0;

    for (std::size_t index = 0; index < limbs_.size(); ++index) {
        const double_limb_type lhsLimb = limbs_[index];
        const double_limb_type rhsLimb =
            index < rhs.limbs_.size() ? rhs.limbs_[index] : 0;
        const double_limb_type subtrahend = rhsLimb + borrow;

        if (lhsLimb >= subtrahend) {
            limbs_[index] = static_cast<limb_type>(lhsLimb - subtrahend);
            borrow = 0;
        }
        else {
            limbs_[index] = static_cast<limb_type>(
                limbBase + lhsLimb - subtrahend);
            borrow = 1;
        }
    }

    normalize();
    return *this;
}

BigUInt& BigUInt::operator*=(const BigUInt& rhs) {
    if (isZero() || rhs.isZero()) {
        limbs_.clear();
        return *this;
    }

    if (limbs_.size() > limbs_.max_size() - rhs.limbs_.size())
        throw std::length_error("BigUInt multiplication result is too large");

    /*
    旧実装では1 limbだけのBigInt×BigIntも、通常の多倍長乗算と同じ一時resultを確保していた。
    factorialのproduct tree末端ではこの形が非常に多いため、既存multiplySmallを直接使って
    allocationと二重loopのsetupを避ける。

    auto product = multiplyKaratsuba(limbs_, rhs.limbs_);
    limbs_ = std::move(product);
    */
    /*
    旧実装では x*x も一般のmultiplyAdaptiveへそのまま流し、非対角項を左右から2回計算していた。
    squareは対称性を使えばlimb乗算をほぼ半減でき、Karatsubaでも平方専用の再帰式を使えるため、
    値が等しいoperandは専用squareへ送る。operator*(lhs, rhs)ではlhsがcopyされるので、
    alias(this == &rhs)だけでは検出できず値の一致で判定する。

    auto product = multiplyAdaptive(limbs_, rhs.limbs_, false);
    limbs_ = std::move(product);
    */
    if (limbs_ == rhs.limbs_) {
        auto squared = squareAdaptive(limbs_, false);
        limbs_ = std::move(squared);
        return *this;
    }

    if (rhs.limbs_.size() == 1) {
        const limb_type factor = rhs.limbs_.front();
        multiplySmall(factor);
        return *this;
    }

    if (limbs_.size() == 1) {
        const limb_type factor = limbs_.front();
        *this = rhs;
        multiplySmall(factor);
        return *this;
    }

    /*
    旧実装はすべてのBigInt×BigIntをO(n^2)のschoolbook法で計算していた。
    巨大階乗のproduct treeでは同程度の大きさの巨大整数同士を繰り返し掛けるため、
    limb数が増えるほどこの二重loopが支配的になる。

    BigUInt result;
    result.limbs_.assign(limbs_.size() + rhs.limbs_.size(), 0);
    for (std::size_t lhsIndex = 0; lhsIndex < limbs_.size(); ++lhsIndex) {
        double_limb_type carry = 0;
        for (std::size_t rhsIndex = 0; rhsIndex < rhs.limbs_.size(); ++rhsIndex) {
            const std::size_t resultIndex = lhsIndex + rhsIndex;
            const double_limb_type product =
                static_cast<double_limb_type>(limbs_[lhsIndex]) * rhs.limbs_[rhsIndex]
                + result.limbs_[resultIndex] + carry;
            result.limbs_[resultIndex] = static_cast<limb_type>(product & limbMask);
            carry = product >> limbBits;
        }
        result.limbs_[lhsIndex + rhs.limbs_.size()] = static_cast<limb_type>(carry);
    }
    */

    /*
    Karatsuba導入直後は、1-limb fast path以外をすべて multiplyKaratsuba へ渡していた。
    さらに巨大なbalanced operandではToom-3の方が実測で速いため、schoolbook / Karatsuba /
    Toom-3をoperand sizeで選ぶadaptive dispatcherへ置き換える。

    auto product = multiplyKaratsuba(limbs_, rhs.limbs_);
    */
    auto product = multiplyAdaptive(limbs_, rhs.limbs_, false);
    limbs_ = std::move(product);
    return *this;
}

BigUInt& BigUInt::operator/=(const BigUInt& rhs) {
    auto result = divmod(*this, rhs);
    *this = std::move(result.quotient);
    return *this;
}

BigUInt& BigUInt::operator%=(const BigUInt& rhs) {
    auto result = divmod(*this, rhs);
    *this = std::move(result.remainder);
    return *this;
}

BigUInt& BigUInt::operator<<=(std::size_t bits) {
    if (isZero() || bits == 0)
        return *this;

    const std::size_t limbShift = bits / limbBits;
    const unsigned bitShift = static_cast<unsigned>(bits % limbBits);
    const std::size_t extraLimb = bitShift == 0 ? 0 : 1;

    // limb単位の移動量と最終carryを含む領域を確保する。
    if (limbShift > limbs_.max_size() - limbs_.size()
        || extraLimb > limbs_.max_size() - limbs_.size() - limbShift)
        throw std::length_error("BigUInt left shift result is too large");

    std::vector<limb_type> shifted(
        limbs_.size() + limbShift + extraLimb,
        0);

    if (bitShift == 0) {
        for (std::size_t index = 0; index < limbs_.size(); ++index)
            shifted[index + limbShift] = limbs_[index];
    }
    else {
        double_limb_type carry = 0;

        for (std::size_t index = 0; index < limbs_.size(); ++index) {
            const double_limb_type current =
                (static_cast<double_limb_type>(limbs_[index]) << bitShift) | carry;

            shifted[index + limbShift] =
                static_cast<limb_type>(current & limbMask);
            carry = current >> limbBits;
        }

        shifted[limbs_.size() + limbShift] =
            static_cast<limb_type>(carry);
    }

    limbs_ = std::move(shifted);
    normalize();
    return *this;
}

BigUInt& BigUInt::operator>>=(std::size_t bits) {
    if (isZero() || bits == 0)
        return *this;

    const std::size_t limbShift = bits / limbBits;
    const unsigned bitShift = static_cast<unsigned>(bits % limbBits);

    if (limbShift >= limbs_.size()) {
        limbs_.clear();
        return *this;
    }

    const std::size_t resultSize = limbs_.size() - limbShift;
    std::vector<limb_type> shifted(resultSize, 0);

    if (bitShift == 0) {
        for (std::size_t index = 0; index < resultSize; ++index)
            shifted[index] = limbs_[index + limbShift];
    }
    else {
        // 各出力limbは、現在の入力limbと一つ上位から流れ込むビットを合成する。
        for (std::size_t index = 0; index < resultSize; ++index) {
            const std::size_t sourceIndex = index + limbShift;
            const double_limb_type current = limbs_[sourceIndex];
            const double_limb_type high =
                sourceIndex + 1 < limbs_.size() ? limbs_[sourceIndex + 1] : 0;

            shifted[index] = static_cast<limb_type>(
                (current >> bitShift)
                | ((high << (limbBits - bitShift)) & limbMask));
        }
    }

    limbs_ = std::move(shifted);
    normalize();
    return *this;
}

std::strong_ordering BigUInt::operator<=>(const BigUInt& rhs) const noexcept {
    if (limbs_.size() < rhs.limbs_.size())
        return std::strong_ordering::less;
    if (limbs_.size() > rhs.limbs_.size())
        return std::strong_ordering::greater;

    for (std::size_t index = limbs_.size(); index-- > 0;) {
        if (limbs_[index] < rhs.limbs_[index])
            return std::strong_ordering::less;
        if (limbs_[index] > rhs.limbs_[index])
            return std::strong_ordering::greater;
    }

    return std::strong_ordering::equal;
}

bool BigUInt::operator==(const BigUInt& rhs) const noexcept {
    return limbs_ == rhs.limbs_;
}

void BigUInt::normalize() noexcept {
    while (!limbs_.empty() && limbs_.back() == 0)
        limbs_.pop_back();
}

void BigUInt::addSmall(limb_type value) {
    if (value == 0)
        return;

    double_limb_type carry = value;
    std::size_t index = 0;

    while (carry != 0 && index < limbs_.size()) {
        const double_limb_type sum =
            static_cast<double_limb_type>(limbs_[index]) + carry;
        limbs_[index] = static_cast<limb_type>(sum & limbMask);
        carry = sum >> limbBits;
        ++index;
    }

    if (carry != 0)
        limbs_.push_back(static_cast<limb_type>(carry));
}

void BigUInt::multiplySmall(limb_type value) {
    if (isZero() || value == 1)
        return;

    if (value == 0) {
        limbs_.clear();
        return;
    }

    double_limb_type carry = 0;

    for (auto& limb : limbs_) {
        const double_limb_type product =
            static_cast<double_limb_type>(limb) * value + carry;

        limb = static_cast<limb_type>(product & limbMask);
        carry = product >> limbBits;
    }

    if (carry != 0)
        limbs_.push_back(static_cast<limb_type>(carry));
}

BigUInt::limb_type BigUInt::divideSmall(limb_type divisor) {
    if (divisor == 0)
        throw std::domain_error("BigUInt division by zero");

    double_limb_type remainder = 0;

    for (std::size_t index = limbs_.size(); index-- > 0;) {
        const double_limb_type current =
            (remainder << limbBits) | limbs_[index];

        limbs_[index] = static_cast<limb_type>(current / divisor);
        remainder = current % divisor;
    }

    normalize();
    return static_cast<limb_type>(remainder);
}

unsigned BigUInt::digitValue(char ch) noexcept {
    if (ch >= '0' && ch <= '9')
        return static_cast<unsigned>(ch - '0');
    if (ch >= 'A' && ch <= 'Z')
        return static_cast<unsigned>(ch - 'A') + 10;
    if (ch >= 'a' && ch <= 'z')
        return static_cast<unsigned>(ch - 'a') + 10;

    return 36;
}

char BigUInt::digitCharacter(unsigned value) noexcept {
    if (value < 10)
        return static_cast<char>('0' + value);

    return static_cast<char>('A' + (value - 10));
}

BigUIntDivModResult divmod(const BigUInt& dividend, const BigUInt& divisor) {
    if (divisor.isZero())
        throw std::domain_error("BigUInt division by zero");

    if (dividend < divisor)
        return {BigUInt{}, dividend};

    if (dividend == divisor)
        return {BigUInt{1}, BigUInt{}};

    /*
    旧実装では巨大な2^k divisorも、この下のKnuth長除算へそのまま流していた。
    2^kでの商は右shift、余りは下位k bitだけなので、長除算を使う必要がない。
    O(n)のbit操作へ落とすことで、巨大なpower-of-two除算と整数算法の補助経路を軽くする。

    // 旧経路: 特別扱いせず divisor.limbs_.size()==1 またはKnuth長除算へ続行
    */
    const std::size_t divisorBitLength = divisor.bitLength();
    const std::size_t divisorTrailingZeros = divisor.trailingZeroBits();
    if (divisorBitLength == divisorTrailingZeros + 1) {
        BigUInt quotient = dividend >> divisorTrailingZeros;
        BigUInt remainder = dividend;
        const std::size_t wholeLimbs = divisorTrailingZeros / limbBits;
        const unsigned remainingBits = static_cast<unsigned>(divisorTrailingZeros % limbBits);

        if (remainingBits == 0) {
            remainder.limbs_.resize(std::min(wholeLimbs, remainder.limbs_.size()));
        }
        else {
            const std::size_t keep = std::min(wholeLimbs + 1, remainder.limbs_.size());
            remainder.limbs_.resize(keep);
            if (wholeLimbs < remainder.limbs_.size()) {
                const BigUInt::limb_type mask = static_cast<BigUInt::limb_type>((std::uint64_t{1} << remainingBits) - 1);
                remainder.limbs_[wholeLimbs] &= mask;
            }
        }
        remainder.normalize();
        return {std::move(quotient), std::move(remainder)};
    }

    if (divisor.limbs_.size() == 1) {
        BigUInt quotient = dividend;
        const auto remainder = quotient.divideSmall(divisor.limbs_.front());
        return {std::move(quotient), BigUInt{remainder}};
    }

    // Knuth式の長除算。正規化により、商の推定補正を原則1回以内に収める。
    const unsigned normalizationShift = static_cast<unsigned>(
        std::countl_zero(divisor.limbs_.back()));

    BigUInt normalizedDividend = dividend << normalizationShift;
    BigUInt normalizedDivisor = divisor << normalizationShift;

    auto dividendLimbs = std::move(normalizedDividend.limbs_);
    const auto& divisorLimbs = normalizedDivisor.limbs_;
    dividendLimbs.push_back(0);

    const std::size_t divisorSize = divisorLimbs.size();
    const std::size_t quotientSize = dividendLimbs.size() - divisorSize;

    BigUInt quotient;
    quotient.limbs_.assign(quotientSize, 0);

    const std::uint64_t divisorHigh = divisorLimbs[divisorSize - 1];
    const std::uint64_t divisorNext = divisorLimbs[divisorSize - 2];

    for (std::size_t offset = quotientSize; offset-- > 0;) {
        const std::uint64_t numerator =
            (static_cast<std::uint64_t>(dividendLimbs[offset + divisorSize])
                << limbBits)
            | dividendLimbs[offset + divisorSize - 1];

        std::uint64_t quotientDigit = numerator / divisorHigh;
        std::uint64_t remainderEstimate = numerator % divisorHigh;

        while (quotientDigit >= limbBase
            || quotientDigit * divisorNext
                > limbBase * remainderEstimate
                    + dividendLimbs[offset + divisorSize - 2]) {
            --quotientDigit;
            remainderEstimate += divisorHigh;

            if (remainderEstimate >= limbBase)
                break;
        }

        if (subtractProduct(
                dividendLimbs, offset, divisorLimbs, quotientDigit)) {
            --quotientDigit;
            addBack(dividendLimbs, offset, divisorLimbs);
        }

        quotient.limbs_[offset] = static_cast<BigUInt::limb_type>(quotientDigit);
    }

    BigUInt remainder;
    remainder.limbs_.assign(
        dividendLimbs.begin(),
        dividendLimbs.begin() + static_cast<std::ptrdiff_t>(divisorSize));

    if (normalizationShift != 0)
        remainder >>= normalizationShift;

    quotient.normalize();
    remainder.normalize();
    return {std::move(quotient), std::move(remainder)};
}

BigUInt operator+(BigUInt lhs, const BigUInt& rhs) {
    lhs += rhs;
    return lhs;
}

BigUInt operator-(BigUInt lhs, const BigUInt& rhs) {
    lhs -= rhs;
    return lhs;
}

BigUInt operator*(BigUInt lhs, const BigUInt& rhs) {
    lhs *= rhs;
    return lhs;
}

BigUInt operator/(BigUInt lhs, const BigUInt& rhs) {
    lhs /= rhs;
    return lhs;
}

BigUInt operator%(BigUInt lhs, const BigUInt& rhs) {
    lhs %= rhs;
    return lhs;
}

BigUInt operator<<(BigUInt value, std::size_t bits) {
    value <<= bits;
    return value;
}

BigUInt operator>>(BigUInt value, std::size_t bits) {
    value >>= bits;
    return value;
}

} // namespace mmcal::numeric::detail
