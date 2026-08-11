#pragma once

namespace mmcal::mathematics {

// 数体系の包含関係を表す。
//     Integer ⊂ Rational ⊂ Real ⊂ Complex
// Unknownは「何も分からない」であり、Complexの別名ではない。
enum class NumericDomain {
    Unknown,
    Integer,
    Rational,
    Real,
    Complex
};

// 実数であることが証明できた式について、その符号までどこまで分かるか。
// 複素数やdomain不明の式ではUnknownを使う。
enum class RealSign {
    Unknown,
    Negative,
    Zero,
    Positive,
    NonPositive,
    NonNegative,
    NonZero
};

// lhsがrhsの部分集合として含まれる数体系ならtrue。
// Unknownは数学的なdomainではないので、Unknown同士でもfalseとする。
[[nodiscard]] constexpr bool isSubdomainOf(
    NumericDomain lhs,
    NumericDomain rhs) noexcept {
    if (lhs == NumericDomain::Unknown || rhs == NumericDomain::Unknown)
        return false;

    return static_cast<int>(lhs) <= static_cast<int>(rhs);
}

} // namespace mmcal::mathematics
