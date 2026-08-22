#pragma once

#include "complex_interval.hpp"

#include <variant>

namespace mmcal::approximation {

// certified interval backendで共通に扱う実/複素数値domain。
// Realは必要時だけComplexへ昇格し，ComplexからRealへの縮約は虚部区間が
// 厳密なpoint zeroと証明できる場合だけ許す。
class CertifiedValue final {
public:
    CertifiedValue(RealInterval real);
    CertifiedValue(ComplexInterval complex);

    [[nodiscard]] bool isReal() const noexcept;
    [[nodiscard]] bool isComplex() const noexcept;
    [[nodiscard]] const RealInterval& asReal() const;
    [[nodiscard]] const ComplexInterval& asComplex() const;
    [[nodiscard]] ComplexInterval toComplex() const;

private:
    std::variant<RealInterval, ComplexInterval> value_;
};

} // namespace mmcal::approximation
