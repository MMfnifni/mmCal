#pragma once

#include "numeric/big_float.hpp"
#include "numeric/rational.hpp"

#include <cstddef>

namespace mmcal::approximation {

// 真の実数値を必ず含む閉区間 [lower, upper]。
//
// 端点BigFloatはそれぞれ厳密なdyadic rationalであり、RealIntervalの責務は「未知の真値がこの範囲から絶対に外れない」という包含保証を維持すること。
// したがって四則演算では下端を -infinity 方向、上端を +infinity 方向へ丸める。
class RealInterval final {
public:
    RealInterval(numeric::BigFloat lower, numeric::BigFloat upper);

    [[nodiscard]] static RealInterval point(numeric::BigFloat value);
    [[nodiscard]] static RealInterval fromRational(
        const numeric::Rational& value,
        std::size_t precisionBits);

    // exactなRational上下界を指定precisionのdyadic端点へ外向きに変換する。
    [[nodiscard]] static RealInterval fromRationalBounds(
        const numeric::Rational& lower,
        const numeric::Rational& upper,
        std::size_t precisionBits);

    [[nodiscard]] const numeric::BigFloat& lower() const noexcept;
    [[nodiscard]] const numeric::BigFloat& upper() const noexcept;

    [[nodiscard]] bool isPoint() const noexcept;
    [[nodiscard]] bool containsZero() const noexcept;
    [[nodiscard]] bool contains(const numeric::Rational& value) const;

    // 同じ真値を含むまま指定precisionへ外向きに丸め直す。
    [[nodiscard]] RealInterval roundedOutward(std::size_t precisionBits) const;

private:
    numeric::BigFloat lower_;
    numeric::BigFloat upper_;
};

[[nodiscard]] RealInterval hull(const RealInterval& lhs, const RealInterval& rhs);
[[nodiscard]] RealInterval negate(const RealInterval& value);

[[nodiscard]] RealInterval add(
    const RealInterval& lhs,
    const RealInterval& rhs,
    std::size_t precisionBits);

[[nodiscard]] RealInterval subtract(
    const RealInterval& lhs,
    const RealInterval& rhs,
    std::size_t precisionBits);

[[nodiscard]] RealInterval multiply(
    const RealInterval& lhs,
    const RealInterval& rhs,
    std::size_t precisionBits);

[[nodiscard]] RealInterval divide(
    const RealInterval& lhs,
    const RealInterval& rhs,
    std::size_t precisionBits);

} // namespace mmcal::approximation
