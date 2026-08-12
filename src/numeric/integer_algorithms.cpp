// 整数平方根・階乗などの整数算法
#include "integer_algorithms.hpp"

#include <bit>
#include <charconv>
#include <limits>
#include <stdexcept>
#include <vector>
#include <system_error>
#include <utility>

namespace mmcal::numeric {
namespace {

[[nodiscard]] BigInt productRange(std::uint64_t first, std::uint64_t last) {
    if (first > last)
        return BigInt{1};

    if (last - first <= 15) {
        BigInt result{1};
        for (std::uint64_t value = first; value <= last; ++value) {
            result *= BigInt::fromUnsigned(value);
            if (value == last)
                break;
        }
        return result;
    }

    const std::uint64_t middle = first + (last - first) / 2;
    return productRange(first, middle) * productRange(middle + 1, last);
}

#ifdef MMCAL_USE_PRIME_SWING_FACTORIAL
[[nodiscard]] std::vector<std::uint64_t> primesUpTo(std::uint64_t n) {
    if (n < 2)
        return {};
    if (n > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max() - 1))
        throw std::length_error("Factorial prime sieve is too large");

    std::vector<bool> composite(static_cast<std::size_t>(n) + 1, false);
    for (std::uint64_t candidate = 2; candidate <= n / candidate; ++candidate) {
        if (composite[static_cast<std::size_t>(candidate)])
            continue;
        for (std::uint64_t multiple = candidate * candidate; multiple <= n; multiple += candidate)
            composite[static_cast<std::size_t>(multiple)] = true;
    }

    std::vector<std::uint64_t> primes;
    for (std::uint64_t candidate = 2; candidate <= n; ++candidate) {
        if (!composite[static_cast<std::size_t>(candidate)])
            primes.push_back(candidate);
    }
    return primes;
}


[[nodiscard]] BigInt productFactors(
    const std::vector<std::uint64_t>& factors,
    std::size_t first,
    std::size_t last) {
    if (last - first <= 15) {
        BigInt result{1};
        for (std::size_t i = first; i <= last; ++i)
            result *= BigInt::fromUnsigned(factors[i]);
        return result;
    }
    const std::size_t middle = first + (last - first) / 2;
    return productFactors(factors, first, middle)
        * productFactors(factors, middle + 1, last);
}

[[nodiscard]] BigInt swingFactor(
    std::uint64_t n,
    const std::vector<std::uint64_t>& primes) {
    if (n < 2)
        return BigInt{1};

    std::vector<std::uint64_t> factors;
    for (const std::uint64_t prime : primes) {
        if (prime > n)
            break;
        if (prime == 2)
            continue;

        // v_p(n!) - 2v_p((n/2)!) は各 p^k について floor(n/p^k) が奇数なら1、
        // 偶数なら0になる。Legendre和を二重計算せず、奇数になった段だけprimeを掛ける。
        std::uint64_t quotient = n;
        std::uint64_t primePower = 1;
        while (quotient != 0) {
            quotient /= prime;
            if ((quotient & 1u) != 0)
                primePower *= prime;
        }
        if (primePower != 1)
            factors.push_back(primePower);
    }

    if (factors.empty())
        return BigInt{1};
    return productFactors(factors, 0, factors.size() - 1);
}

[[nodiscard]] BigInt oddFactorial(
    std::uint64_t n,
    const std::vector<std::uint64_t>& primes) {
    if (n < 2)
        return BigInt{1};

    BigInt half = oddFactorial(n / 2, primes);
    return half * half * swingFactor(n, primes);
}

#endif

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

    /*
    旧実装はbalanced product treeで 2..n をすべて掛けていた。逐次乗算より十分速いが、
    2の因子まで巨大中間値へ抱えたまま計算する。Prime-Swingでは2の冪を最後のshiftへ
    分離し、奇数部だけを再帰的に構築することで巨大乗算の総仕事量を減らす。

    return productRange(2, n);
    */
#ifdef MMCAL_USE_PRIME_SWING_FACTORIAL
    // Prime-Swingも実装・比較できるよう残すが、現行BigIntではproduct treeより遅かったため既定にはしない。
    const auto primes = primesUpTo(n);
    BigInt result = oddFactorial(n, primes);
    const std::uint64_t powerOfTwo = n - std::popcount(n);
    if (powerOfTwo > static_cast<std::uint64_t>(std::numeric_limits<std::size_t>::max()))
        throw std::length_error("Factorial shift is too large");
    result <<= static_cast<std::size_t>(powerOfTwo);
    return result;
#else
    return productRange(2, n);
#endif
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
