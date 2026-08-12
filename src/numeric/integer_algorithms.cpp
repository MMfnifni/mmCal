// 整数平方根・階乗などの整数算法
#include "integer_algorithms.hpp"

#include <charconv>
#include <stdexcept>
#include <system_error>
#include <utility>

namespace mmcal::numeric {
namespace {

[[nodiscard]] BigInt productRange(std::uint64_t first, std::uint64_t last) {
    if (first > last)
        return BigInt{1};

    // 小区間は再帰よりstraight loopの方が軽い。大区間だけproduct treeにして、極端に大きさの違うBigIntを順次掛け続ける形を避ける。
    if (last - first <= 15) {
        BigInt result{1};
        for (std::uint64_t value = first; value <= last; ++value) {
            /*
            旧実装:
            result *= BigInt::parse(std::to_string(value));

            product treeの葉ごとにuint64_t→decimal文字列→BigUIntと往復していた。
            算術上不要な変換なので、uint64_tからmagnitudeを直接構築する。
            */
            result *= BigInt::fromUnsigned(value);
            if (value == last)
                break;
        }
        return result;
    }

    const std::uint64_t middle = first + (last - first) / 2;
    return productRange(first, middle) * productRange(middle + 1, last);
}

} // namespace

BigInt gcd(BigInt lhs, BigInt rhs) {
    lhs = lhs.abs();
    rhs = rhs.abs();

    while (!rhs.isZero()) {
        lhs %= rhs;
        std::swap(lhs, rhs);
    }

    return lhs;
}

BigInt lcm(const BigInt& lhs, const BigInt& rhs) {
    if (lhs.isZero() || rhs.isZero())
        return BigInt{};

    return ((lhs / gcd(lhs, rhs)) * rhs).abs();
}

BigInt pow(BigInt base, std::uint64_t exponent) {
    BigInt result{1};

    // 二乗しながら累乗し、乗算回数を指数の対数程度に抑える。
    while (exponent != 0) {
        if ((exponent & 1u) != 0)
            result *= base;

        exponent >>= 1;
        if (exponent != 0)
            base *= base;
    }

    return result;
}

BigInt factorial(std::uint64_t n) {
    if (n < 2)
        return BigInt{1};
    return productRange(2, n);
}

std::optional<std::uint64_t> tryToUint64(const BigInt& value) {
    if (value.isNegative())
        return std::nullopt;

    const std::string text = value.toString();
    std::uint64_t result = 0;
    const auto conversion = std::from_chars(text.data(), text.data() + text.size(), result);
    if (conversion.ec != std::errc{} || conversion.ptr != text.data() + text.size())
        return std::nullopt;
    return result;
}

IntegerSqrtResult integerSqrt(const BigInt& value) {
    if (value.isNegative())
        throw std::domain_error("integerSqrt requires a nonnegative value");

    if (value.isZero())
        return {BigInt{}, BigInt{}};
    if (value == BigInt{1})
        return {BigInt{1}, BigInt{}};

    const auto initialBit = (value.magnitude_.bitLength() + 1) / 2;
    detail::BigUInt guessMagnitude{1};
    guessMagnitude <<= initialBit;
    BigInt guess{std::move(guessMagnitude), false};

    // Newton法で上側から収束させ、floor(sqrt(value)) を求める。
    while (true) {
        BigInt next = (guess + value / guess) / BigInt{2};
        if (next >= guess)
            break;
        guess = std::move(next);
    }

    BigInt remainder = value - guess * guess;
    return {std::move(guess), std::move(remainder)};
}


IntegerCubeRootResult integerCubeRoot(const BigInt& value) {
    if (value.isNegative())
        throw std::domain_error("integerCubeRoot requires a nonnegative value");
    if (value.isZero())
        return {BigInt{}, BigInt{}};
    if (value == BigInt{1})
        return {BigInt{1}, BigInt{}};

    // sqrtと同様にbit lengthから上側初期値を作り、整数Newton法で収束させる。
    // 二分探索のO(bitLength)回の巨大乗算を避け、巨大BigIntでも反復回数を抑える。
    const std::size_t initialBit = (value.bitLength() + 2) / 3;
    BigInt guess{1};
    guess <<= initialBit;

    while (true) {
        const BigInt square = guess * guess;
        BigInt next = (BigInt{2} * guess + value / square) / BigInt{3};
        if (next >= guess)
            break;
        guess = std::move(next);
    }

    // 整数Newtonの停止点をfloor(cuberoot(value))へ厳密に補正する。
    BigInt cube = guess * guess * guess;
    while (cube > value) {
        guess -= BigInt{1};
        cube = guess * guess * guess;
    }
    while (true) {
        const BigInt next = guess + BigInt{1};
        const BigInt nextCube = next * next * next;
        if (nextCube > value)
            break;
        guess = next;
        cube = nextCube;
    }

    return {guess, value - cube};
}

bool isPerfectSquare(const BigInt& value) {
    if (value.isNegative())
        return false;
    return integerSqrt(value).remainder.isZero();
}

} // namespace mmcal::numeric
