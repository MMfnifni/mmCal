#include "random_engine.hpp"

#include <bit>
#include <chrono>
#include <random>
#include <stdexcept>
#include <string>

namespace mmcal::random {
namespace {

[[nodiscard]] numeric::BigInt fromUint64(std::uint64_t value) {
    return numeric::BigInt::parse(std::to_string(value));
}

[[nodiscard]] std::uint64_t stableSeedHash(const numeric::BigInt& value) {
    constexpr std::uint64_t offset = 14695981039346656037ULL;
    constexpr std::uint64_t prime = 1099511628211ULL;

    std::uint64_t hash = offset;
    for (const unsigned char character : value.toString()) {
        hash ^= character;
        hash *= prime;
    }
    return hash;
}

[[nodiscard]] std::uint64_t splitMix64(std::uint64_t& state) noexcept {
    std::uint64_t value = (state += 0x9E3779B97F4A7C15ULL);
    value = (value ^ (value >> 30U)) * 0xBF58476D1CE4E5B9ULL;
    value = (value ^ (value >> 27U)) * 0x94D049BB133111EBULL;
    return value ^ (value >> 31U);
}

[[nodiscard]] std::uint64_t entropyWord() {
    std::random_device device;
    const std::uint64_t high = static_cast<std::uint64_t>(device()) << 32U;
    const std::uint64_t low = static_cast<std::uint64_t>(device());
    const std::uint64_t clock = static_cast<std::uint64_t>(
        std::chrono::high_resolution_clock::now().time_since_epoch().count());
    return high ^ low ^ std::rotl(clock, 17);
}

} // namespace

RandomEngine::RandomEngine() {
    static_cast<void>(reseedFromEntropy());
}

void RandomEngine::seed(const numeric::BigInt& seedValue) {
    std::uint64_t state = stableSeedHash(seedValue);
    for (std::uint64_t& word : state_)
        word = splitMix64(state);

    // xoshiro256**は全zero stateだけを禁止する。splitmixでは実質起こらないが防御する。
    if (state_[0] == 0 && state_[1] == 0 && state_[2] == 0 && state_[3] == 0)
        state_[0] = 1;
}

numeric::BigInt RandomEngine::reseedFromEntropy() {
    const numeric::BigInt seedValue = fromUint64(entropyWord());
    seed(seedValue);
    return seedValue;
}

std::uint64_t RandomEngine::nextUint64() noexcept {
    const std::uint64_t result = std::rotl(state_[1] * 5ULL, 7) * 9ULL;
    const std::uint64_t temporary = state_[1] << 17U;

    state_[2] ^= state_[0];
    state_[3] ^= state_[1];
    state_[1] ^= state_[2];
    state_[0] ^= state_[3];
    state_[2] ^= temporary;
    state_[3] = std::rotl(state_[3], 45);
    return result;
}

numeric::BigInt RandomEngine::uniformBelow(const numeric::BigInt& exclusiveUpper) {
    if (exclusiveUpper <= numeric::BigInt{})
        throw std::invalid_argument("Random upper bound must be positive");
    if (exclusiveUpper == numeric::BigInt{1})
        return numeric::BigInt{};

    const std::size_t bits = exclusiveUpper.bitLength();
    for (;;) {
        numeric::BigInt candidate{};
        std::size_t shift = 0;
        std::size_t remaining = bits;

        while (remaining != 0) {
            const std::size_t take = remaining < 64 ? remaining : 64;
            std::uint64_t word = nextUint64();
            if (take < 64)
                word &= (std::uint64_t{1} << take) - 1ULL;

            numeric::BigInt part = fromUint64(word);
            part <<= shift;
            candidate += part;
            shift += take;
            remaining -= take;
        }

        if (candidate < exclusiveUpper)
            return candidate;
    }
}

numeric::Rational RandomEngine::unit53() {
    constexpr std::uint64_t mask53 = (std::uint64_t{1} << 53U) - 1ULL;
    const std::uint64_t numerator = (nextUint64() >> 11U) & mask53;
    numeric::BigInt denominator{1};
    denominator <<= 53;
    return numeric::Rational{fromUint64(numerator), std::move(denominator)};
}

numeric::Rational RandomEngine::positiveUnit53() {
    constexpr std::uint64_t mask53 = (std::uint64_t{1} << 53U) - 1ULL;
    const std::uint64_t numerator = ((nextUint64() >> 11U) & mask53) + 1ULL;
    numeric::BigInt denominator{1};
    denominator <<= 53;
    return numeric::Rational{fromUint64(numerator), std::move(denominator)};
}

} // namespace mmcal::random
