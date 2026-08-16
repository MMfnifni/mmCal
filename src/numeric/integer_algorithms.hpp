#pragma once

#include "big_int.hpp"

#include <cstdint>
#include <optional>
#include <vector>

namespace mmcal::numeric {

struct IntegerSqrtResult final {
    BigInt root;
    BigInt remainder;
};

struct IntegerCubeRootResult final {
    BigInt root;
    BigInt remainder;
};

[[nodiscard]] BigInt gcd(BigInt lhs, BigInt rhs);
[[nodiscard]] BigInt lcm(const BigInt& lhs, const BigInt& rhs);
[[nodiscard]] BigInt pow(BigInt base, std::uint64_t exponent);
[[nodiscard]] BigInt factorial(std::uint64_t n);
[[nodiscard]] std::optional<std::uint64_t> tryToUint64(const BigInt& value);
[[nodiscard]] IntegerSqrtResult integerSqrt(const BigInt& value);
[[nodiscard]] IntegerCubeRootResult integerCubeRoot(const BigInt& value);
[[nodiscard]] bool isPerfectSquare(const BigInt& value);
[[nodiscard]] bool isPrimeUint64(std::uint64_t value) noexcept;
[[nodiscard]] bool factorUint64(std::uint64_t value, std::vector<std::uint64_t>& factors);

} // namespace mmcal::numeric
