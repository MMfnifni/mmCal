#pragma once

#include "../numeric/rounding_mode.hpp"

#include <cstddef>

namespace mmcal::approximation {

using RoundingMode = numeric::RoundingMode;

// 近似値そのものとは分離し、要求精度と内部ガード桁だけを保持する。
class ApproximationContext final {
public:
    static constexpr std::size_t defaultDecimalDigits = 16;
    static constexpr std::size_t defaultGuardDigits = 8;

    ApproximationContext() = default;
    explicit ApproximationContext(
        std::size_t decimalDigits,
        std::size_t guardDigits = defaultGuardDigits,
        RoundingMode roundingMode = RoundingMode::NearestEven);

    [[nodiscard]] std::size_t decimalDigits() const noexcept;
    [[nodiscard]] std::size_t guardDigits() const noexcept;
    [[nodiscard]] std::size_t workingDecimalDigits() const;
    [[nodiscard]] std::size_t workingBinaryBits() const;
    [[nodiscard]] RoundingMode roundingMode() const noexcept;

    void setDecimalDigits(std::size_t decimalDigits);
    void setGuardDigits(std::size_t guardDigits) noexcept;
    void setRoundingMode(RoundingMode roundingMode) noexcept;

private:
    std::size_t decimalDigits_ = defaultDecimalDigits;
    std::size_t guardDigits_ = defaultGuardDigits;
    RoundingMode roundingMode_ = RoundingMode::NearestEven;
};

} // namespace mmcal::approximation
