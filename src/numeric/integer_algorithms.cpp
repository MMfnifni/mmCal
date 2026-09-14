// 整数平方根・階乗などの整数算法
#include "integer_algorithms.hpp"

#include <algorithm>
#include <bit>
#include <charconv>
#include <limits>
#include <map>
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


using FactorMap = std::map<BigInt, std::uint64_t>;

struct BigFactorContext final {
    std::size_t rhoStepsRemaining = 4'000'000;
    unsigned proofDepthLimit = 24;
};

[[nodiscard]] const std::vector<std::uint32_t>& smallTrialPrimes() {
    static const std::vector<std::uint32_t> primes = [] {
        constexpr std::uint32_t limit = 10'000;
        std::vector<bool> composite(static_cast<std::size_t>(limit) + 1, false);
        std::vector<std::uint32_t> result;
        for (std::uint32_t candidate = 2; candidate <= limit; ++candidate) {
            if (composite[candidate])
                continue;
            result.push_back(candidate);
            if (candidate > limit / candidate)
                continue;
            for (std::uint32_t multiple = candidate * candidate; multiple <= limit; multiple += candidate)
                composite[multiple] = true;
        }
        return result;
    }();
    return primes;
}

[[nodiscard]] BigInt absDifference(const BigInt& lhs, const BigInt& rhs) {
    return lhs >= rhs ? lhs - rhs : rhs - lhs;
}

[[nodiscard]] BigInt powerModBigInt(BigInt base, BigInt exponent, const BigInt& modulus) {
    BigInt result{1};
    base %= modulus;
    while (!exponent.isZero()) {
        if (exponent.modulo(2) != 0)
            result = (result * base) % modulus;
        exponent >>= 1;
        if (!exponent.isZero())
            base = (base * base) % modulus;
    }
    return result;
}

[[nodiscard]] bool isProbablePrimeBigInt(const BigInt& value) {
    if (value < BigInt{2})
        return false;
    if (const auto small = tryToUint64(value))
        return isPrimeUint64(*small);

    constexpr std::uint32_t trialPrimes[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    };
    for (const std::uint32_t prime : trialPrimes) {
        if (value.modulo(prime) == 0)
            return false;
    }

    BigInt d = value - BigInt{1};
    const std::size_t s = d.trailingZeroBits();
    d >>= s;

    // ここでは素数証明ではなくcomposite witness探索にだけ使う。
    // 全baseを通過しても後段のPocklington証明なしにはprimeとして確定しない。
    constexpr std::uint32_t bases[] = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37
    };
    const BigInt minusOne = value - BigInt{1};
    for (const std::uint32_t baseValue : bases) {
        BigInt x = powerModBigInt(BigInt::fromUnsigned(baseValue), d, value);
        if (x == BigInt{1} || x == minusOne)
            continue;

        bool strongProbablePrime = false;
        for (std::size_t r = 1; r < s; ++r) {
            x = (x * x) % value;
            if (x == minusOne) {
                strongProbablePrime = true;
                break;
            }
            if (x == BigInt{1})
                return false;
        }
        if (!strongProbablePrime)
            return false;
    }
    return true;
}

[[nodiscard]] BigInt rhoStep(const BigInt& value, const BigInt& c, const BigInt& modulus) {
    return (value * value + c) % modulus;
}

[[nodiscard]] std::optional<BigInt> pollardRhoBigInt(
    const BigInt& value,
    BigFactorContext& context) {
    if (value.modulo(2) == 0)
        return BigInt{2};
    if (value.modulo(3) == 0)
        return BigInt{3};

    // Brent法でGCD回数を抑える。全探索量はcontextで上限を共有し、
    // 巨大な難しい入力が無制限に走り続けないようにする。
    constexpr std::size_t batchSize = 64;
    for (std::uint64_t attempt = 1; attempt <= 32 && context.rhoStepsRemaining != 0; ++attempt) {
        const BigInt c = BigInt::fromUnsigned(2 * attempt + 1);
        BigInt y = BigInt::fromUnsigned(2 + attempt);
        BigInt g{1};
        BigInt x;
        BigInt ys;
        std::size_t r = 1;

        while (g == BigInt{1} && context.rhoStepsRemaining != 0) {
            x = y;
            for (std::size_t i = 0; i < r && context.rhoStepsRemaining != 0; ++i) {
                y = rhoStep(y, c, value);
                --context.rhoStepsRemaining;
            }
            if (context.rhoStepsRemaining == 0)
                break;

            std::size_t k = 0;
            BigInt q{1};
            while (k < r && g == BigInt{1} && context.rhoStepsRemaining != 0) {
                ys = y;
                const std::size_t count = std::min(batchSize, r - k);
                for (std::size_t i = 0; i < count && context.rhoStepsRemaining != 0; ++i) {
                    y = rhoStep(y, c, value);
                    --context.rhoStepsRemaining;
                    const BigInt difference = absDifference(x, y);
                    q = (q * difference) % value;
                }
                g = gcd(q, value);
                k += count;
            }

            if (g == value) {
                do {
                    if (context.rhoStepsRemaining == 0)
                        break;
                    ys = rhoStep(ys, c, value);
                    --context.rhoStepsRemaining;
                    g = gcd(absDifference(x, ys), value);
                } while (g == BigInt{1});
            }

            if (g > BigInt{1} && g < value)
                return g;
            if (g == value)
                break;
            if (r > (std::numeric_limits<std::size_t>::max() >> 1))
                break;
            r <<= 1;
        }
    }
    return std::nullopt;
}

[[nodiscard]] bool addFactor(
    FactorMap& factors,
    const BigInt& prime,
    std::uint64_t exponent = 1) {
    auto& current = factors[prime];
    if (current > std::numeric_limits<std::uint64_t>::max() - exponent)
        return false;
    current += exponent;
    return true;
}

[[nodiscard]] bool factorBigIntInternal(
    BigInt value,
    FactorMap& factors,
    BigFactorContext& context,
    unsigned proofDepth);

[[nodiscard]] bool provePrimeBigInt(
    const BigInt& value,
    BigFactorContext& context,
    unsigned proofDepth) {
    if (const auto small = tryToUint64(value))
        return isPrimeUint64(*small);
    if (proofDepth > context.proofDepthLimit || value.bitLength() > 512)
        return false;
    if (!isProbablePrimeBigInt(value))
        return false;

    // Pocklingtonを完全因数分解したn-1に適用する。probable-prime判定は
    // compositeを早く落とすためだけで、prime確定には必ずこの証明を通す。
    FactorMap predecessorFactors;
    if (!factorBigIntInternal(value - BigInt{1}, predecessorFactors, context, proofDepth + 1))
        return false;

    const BigInt exponent = value - BigInt{1};
    for (const auto& [prime, ignoredExponent] : predecessorFactors) {
        static_cast<void>(ignoredExponent);
        const BigInt reducedExponent = exponent / prime;
        bool foundWitness = false;
        for (std::uint64_t witnessValue = 2; witnessValue <= 64; ++witnessValue) {
            const BigInt witness = BigInt::fromUnsigned(witnessValue);
            if (powerModBigInt(witness, exponent, value) != BigInt{1})
                continue;
            const BigInt residue = powerModBigInt(witness, reducedExponent, value);
            if (gcd(residue - BigInt{1}, value) == BigInt{1}) {
                foundWitness = true;
                break;
            }
        }
        if (!foundWitness)
            return false;
    }
    return true;
}

[[nodiscard]] bool stripSmallPrimeFactors(BigInt& value, FactorMap& factors) {
    for (const std::uint32_t prime : smallTrialPrimes()) {
        if (value == BigInt{1})
            return true;
        std::uint64_t exponent = 0;
        while (value.modulo(prime) == 0) {
            value /= BigInt::fromUnsigned(prime);
            ++exponent;
        }
        if (exponent != 0 && !addFactor(factors, BigInt::fromUnsigned(prime), exponent))
            return false;
        if (const auto small = tryToUint64(value)) {
            if (*small < static_cast<std::uint64_t>(prime) * prime)
                break;
        }
    }
    return true;
}

[[nodiscard]] bool factorBigIntInternal(
    BigInt value,
    FactorMap& factors,
    BigFactorContext& context,
    unsigned proofDepth) {
    if (value == BigInt{1})
        return true;
    if (!stripSmallPrimeFactors(value, factors))
        return false;
    if (value == BigInt{1})
        return true;

    if (const auto small = tryToUint64(value)) {
        std::vector<std::uint64_t> smallFactors;
        if (!factorUint64(*small, smallFactors))
            return false;
        for (const std::uint64_t factor : smallFactors) {
            if (!addFactor(factors, BigInt::fromUnsigned(factor)))
                return false;
        }
        return true;
    }

    // 大きな完全平方はrhoより先に割る。巨大なprime powerにも効く。
    const auto square = integerSqrt(value);
    if (square.remainder.isZero()) {
        FactorMap rootFactors;
        if (!factorBigIntInternal(square.root, rootFactors, context, proofDepth))
            return false;
        for (const auto& [prime, exponent] : rootFactors) {
            if (exponent > std::numeric_limits<std::uint64_t>::max() / 2
                || !addFactor(factors, prime, exponent * 2))
                return false;
        }
        return true;
    }

    if (isProbablePrimeBigInt(value)) {
        if (!provePrimeBigInt(value, context, proofDepth + 1))
            return false;
        return addFactor(factors, value);
    }

    const auto divisor = pollardRhoBigInt(value, context);
    if (!divisor)
        return false;
    return factorBigIntInternal(*divisor, factors, context, proofDepth)
        && factorBigIntInternal(value / *divisor, factors, context, proofDepth);
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


bool factorBigInt(const BigInt& value, std::vector<PrimePowerFactor>& factors) {
    factors.clear();
    if (value.isZero() || value.isNegative())
        return false;
    if (value == BigInt{1})
        return true;

    FactorMap grouped;
    BigFactorContext context;
    if (!factorBigIntInternal(value, grouped, context, 0))
        return false;

    factors.reserve(grouped.size());
    for (auto& [prime, exponent] : grouped)
        factors.push_back(PrimePowerFactor{std::move(prime), exponent});
    return true;
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
