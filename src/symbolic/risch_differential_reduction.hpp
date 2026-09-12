#pragma once

#include "symbolic/differential_tower.hpp"
#include "symbolic/risch_core.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace mmcal::symbolic::risch {

// K[t]，K=Q(x)。coefficients[i]はt^iの係数である。
class DifferentialPolynomial final {
public:
    DifferentialPolynomial();
    explicit DifferentialPolynomial(std::vector<RationalFunction> coefficients);

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] std::size_t degree() const noexcept;
    [[nodiscard]] const RationalFunction& coefficient(
        std::size_t exponent) const noexcept;
    [[nodiscard]] const std::vector<RationalFunction>& coefficients() const noexcept;

private:
    std::vector<RationalFunction> coefficients_;
    void normalize();
};

struct DifferentialDerivation final {
    DifferentialExtensionKind kind = DifferentialExtensionKind::Primitive;
    // PrimitiveではDt，ExponentialではDt/t。
    RationalFunction differentialCoefficient;
};

enum class DifferentialPolynomialClass {
    Constant,
    Normal,
    Special,
    Mixed
};

struct DifferentialPolynomialClassification final {
    DifferentialPolynomialClass classification =
        DifferentialPolynomialClass::Constant;
    DifferentialPolynomial gcdWithDerivative;
    bool exactVerified = false;
};

struct DifferentialNormalReduction final {
    // numerator/f^originalPower = sum D(B_i/f^power_i)
    //                               + squareFreeNumerator/f。
    std::vector<std::pair<DifferentialPolynomial, std::size_t>> rationalTerms;
    DifferentialPolynomial squareFreeNumerator;
    std::size_t steps = 0;
    bool exactVerified = false;
};

[[nodiscard]] RischStageResult<DifferentialPolynomial>
differentiateDifferentialPolynomial(
    const DifferentialPolynomial& polynomial,
    const DifferentialDerivation& derivation,
    const RischOptions& options = {});

[[nodiscard]] RischStageResult<DifferentialPolynomialClassification>
classifyDifferentialPolynomial(
    const DifferentialPolynomial& polynomial,
    const DifferentialDerivation& derivation,
    const RischOptions& options = {});

// normalなfactor fについてA/f^kをsquare-free denominatorまで下げる。
// special/mixed factorはNonNormalDifferentialFactorで停止し，非初等とは断定しない。
[[nodiscard]] RischStageResult<DifferentialNormalReduction>
hermiteReduceNormalDifferentialPower(
    DifferentialPolynomial numerator,
    const DifferentialPolynomial& factor,
    std::size_t denominatorPower,
    const DifferentialDerivation& derivation,
    const RischOptions& options = {});

[[nodiscard]] bool verifyNormalDifferentialReduction(
    const DifferentialPolynomial& originalNumerator,
    const DifferentialPolynomial& factor,
    std::size_t originalPower,
    const DifferentialDerivation& derivation,
    const DifferentialNormalReduction& reduction,
    const RischOptions& options = {});

[[nodiscard]] bool equivalentDifferentialPolynomials(
    const DifferentialPolynomial& lhs,
    const DifferentialPolynomial& rhs);

} // namespace mmcal::symbolic::risch
