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
static_assert(karatsubaThresholdLimbs >= 1, "Karatsuba threshold must be at least one limb");

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

[[nodiscard]] std::vector<Limb> multiplyKaratsuba(
    std::span<const Limb> lhs,
    std::span<const Limb> rhs) {
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

    auto z0 = multiplyKaratsuba(lhsLow, rhsLow);
    auto z2 = multiplyKaratsuba(lhsHigh, rhsHigh);
    const auto lhsSum = addLimbs(lhsLow, lhsHigh);
    const auto rhsSum = addLimbs(rhsLow, rhsHigh);
    auto z1 = multiplyKaratsuba(lhsSum, rhsSum);
    subtractLimbsInPlace(z1, z0);
    subtractLimbsInPlace(z1, z2);

    std::vector<Limb> result(lhs.size() + rhs.size(), 0);
    addShiftedLimbs(result, z0, 0);
    addShiftedLimbs(result, z1, split);
    addShiftedLimbs(result, z2, split * 2);
    normalizeLimbs(result);
    return result;
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
    std::string result;
    result.reserve((bitLength() + 2) / 3);

    while (!remaining.isZero()) {
        const auto remainder =
            remaining.divideSmall(static_cast<limb_type>(radix));
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

    // 小さい積は従来のschoolbook、大きく釣り合った積はKaratsubaへ自動分岐する。
    // thresholdはtools/BigInt_benchmarkで実測し、balanced乗算のcrossover付近である48 limbsに置く。
    auto product = multiplyKaratsuba(limbs_, rhs.limbs_);
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
