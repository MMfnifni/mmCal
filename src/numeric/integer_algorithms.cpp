// 整数平方根・階乗などの整数算法
#include "integer_algorithms.hpp"

#include <bit>
#include <charconv>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>
#include <system_error>
#include <utility>

namespace mmcal::numeric {
namespace {

[[nodiscard]] std::uint64_t addMod64(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::uint64_t modulus) noexcept {
    return lhs >= modulus - rhs ? lhs - (modulus - rhs) : lhs + rhs;
}

[[nodiscard]] std::uint64_t multiplyMod64(
    std::uint64_t lhs,
    std::uint64_t rhs,
    std::uint64_t modulus) noexcept {
    std::uint64_t result = 0;
    lhs %= modulus;
    while (rhs != 0) {
        if ((rhs & 1U) != 0)
            result = addMod64(result, lhs, modulus);
        rhs >>= 1U;
        if (rhs != 0)
            lhs = addMod64(lhs, lhs, modulus);
    }
    return result;
}

[[nodiscard]] std::uint64_t powerMod64(
    std::uint64_t base,
    std::uint64_t exponent,
    std::uint64_t modulus) noexcept {
    std::uint64_t result = 1 % modulus;
    base %= modulus;
    while (exponent != 0) {
        if ((exponent & 1U) != 0)
            result = multiplyMod64(result, base, modulus);
        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyMod64(base, base, modulus);
    }
    return result;
}

[[nodiscard]] std::uint64_t rhoPolynomial(
    std::uint64_t x,
    std::uint64_t c,
    std::uint64_t modulus) noexcept {
    return addMod64(multiplyMod64(x, x, modulus), c % modulus, modulus);
}

[[nodiscard]] std::optional<std::uint64_t> pollardRho64(std::uint64_t value) {
    if (value % 2 == 0)
        return 2;
    if (value % 3 == 0)
        return 3;

    for (std::uint64_t c = 1; c < 128; ++c) {
        std::uint64_t x = (2 + c) % value;
        std::uint64_t y = x;
        std::uint64_t divisor = 1;
        for (std::size_t iteration = 0; iteration < 2'000'000 && divisor == 1; ++iteration) {
            x = rhoPolynomial(x, c, value);
            y = rhoPolynomial(rhoPolynomial(y, c, value), c, value);
            const std::uint64_t difference = x > y ? x - y : y - x;
            divisor = std::gcd(difference, value);
        }
        if (divisor > 1 && divisor < value)
            return divisor;
    }
    return std::nullopt;
}

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


bool isPrimeUint64(std::uint64_t value) noexcept {
    if (value < 2)
        return false;
    constexpr std::uint64_t smallPrimes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    };
    for (const std::uint64_t prime : smallPrimes) {
        if (value == prime)
            return true;
        if (value % prime == 0)
            return false;
    }

    std::uint64_t d = value - 1;
    unsigned s = 0;
    while ((d & 1U) == 0) {
        d >>= 1U;
        ++s;
    }
    for (const std::uint64_t base : smallPrimes) {
        if (base >= value)
            continue;
        std::uint64_t x = powerMod64(base, d, value);
        if (x == 1 || x == value - 1)
            continue;
        bool witness = true;
        for (unsigned r = 1; r < s; ++r) {
            x = multiplyMod64(x, x, value);
            if (x == value - 1) {
                witness = false;
                break;
            }
        }
        if (witness)
            return false;
    }
    return true;
}

bool factorUint64(std::uint64_t value, std::vector<std::uint64_t>& factors) {
    if (value == 1)
        return true;
    if (isPrimeUint64(value)) {
        factors.push_back(value);
        return true;
    }
    const auto divisor = pollardRho64(value);
    if (!divisor)
        return false;
    return factorUint64(*divisor, factors)
        && factorUint64(value / *divisor, factors);
}

BigInt gcd(BigInt lhs, BigInt rhs) {
    lhs = lhs.abs();
    rhs = rhs.abs();

    // binary GCDも比較したが、現BigIntでは巨大random operandでEuclid法より大幅に遅かった。
    // Knuth除算の改善前に置き換えるのは逆効果なので、既存Euclid法を維持する。
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

    /*
    旧実装は範囲外の巨大整数まで10進文字列へ変換してからfrom_charsしていた。
    100000bit級では「uint64_tに入らない」と判定するだけのために巨大な10進変換が走る。
    uint64_tは最大64bitなので、それを超える値は表示変換なしで即座に棄却する。

    const std::string text = value.toString();
    std::uint64_t result = 0;
    const auto conversion = std::from_chars(text.data(), text.data() + text.size(), result);
    */
    if (value.bitLength() > 64)
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
