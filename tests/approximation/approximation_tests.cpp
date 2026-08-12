// approximationの回帰テスト
#include "approximation_tests.hpp"

#include "approximation/approximation_context.hpp"
#include "approximation/certified_trigonometry.hpp"
#include "approximation/certified_evaluator.hpp"
#include "evaluation/builtin_registry.hpp"
#include "mathematics/math_registry.hpp"
#include "mathematics/angle.hpp"
#include "symbols/symbol_table.hpp"
#include "numeric/number.hpp"
#include "approximation/complex_interval.hpp"
#include "approximation/precision.hpp"
#include "numeric/big_int.hpp"
#include "numeric/real_number.hpp"
#include "test_framework.hpp"

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace mmcal::tests {

void runApproximationTests(TestRunner& tests) {
    using approximation::ApproximationContext;
    using approximation::RoundingMode;
    using numeric::BigInt;
    using numeric::RealNumber;

    const ApproximationContext defaults;
    tests.expectEqual(defaults.decimalDigits(), std::size_t{16},
        "ApproximationContext: default requested digits");
    tests.expectEqual(defaults.guardDigits(), std::size_t{8},
        "ApproximationContext: default guard digits");
    tests.expectEqual(defaults.workingDecimalDigits(), std::size_t{24},
        "ApproximationContext: working digits include guard digits");
    tests.expectEqual(defaults.workingBinaryBits(), std::size_t{80},
        "ApproximationContext: converts working decimal digits to safe binary bits");
    tests.expectEqual(approximation::decimalDigitsToBinaryBits(100), std::size_t{333},
        "Approximation precision: 100 decimal digits require a safe 333-bit bound");

    ApproximationContext context{50, 12, RoundingMode::TowardZero};
    tests.expectEqual(context.decimalDigits(), std::size_t{50},
        "ApproximationContext: stores requested digits");
    tests.expectEqual(context.workingDecimalDigits(), std::size_t{62},
        "ApproximationContext: stores explicit guard digits");
    tests.expect(context.roundingMode() == RoundingMode::TowardZero,
        "ApproximationContext: stores rounding mode");

    tests.expectThrows<std::invalid_argument>([] {
        static_cast<void>(ApproximationContext{0});
    }, "ApproximationContext: rejects zero precision");

    context.setDecimalDigits(std::numeric_limits<std::size_t>::max());
    context.setGuardDigits(1);
    tests.expectThrows<std::overflow_error>([&] {
        static_cast<void>(context.workingDecimalDigits());
    }, "ApproximationContext: detects working precision overflow");

    // ここから先はfloat/double/<cmath>を一切使わないcertified trigonometryの試験。
    const auto sin16 = approximation::approximateSin(RealNumber{BigInt{1}}, 16);
    tests.expectEqual(std::string{sin16.value.text()},
        std::string{"0.8414709848078965"},
        "CertifiedTrig: sin(1) at 16 fractional digits");

    const auto sin100 = approximation::approximateSin(RealNumber{BigInt{1}}, 100);
    tests.expectEqual(std::string{sin100.value.text()},
        std::string{"0.8414709848078965066525023216302989996225630607983710656727517099919104043912396689486397435430526959"},
        "CertifiedTrig: sin(1) at 100 fractional digits");

    // 100桁専用の固定反復ではないことを確認するため、半端な137桁でも検証する。
    const auto sin137 = approximation::approximateSin(RealNumber{BigInt{1}}, 137);
    tests.expectEqual(std::string{sin137.value.text()},
        std::string{"0.84147098480789650665250232163029899962256306079837106567275170999191040439123966894863974354305269585434903790792067429325911892099189888"},
        "CertifiedTrig: sin(1) adapts to an arbitrary 137-digit request");
    tests.expect(sin137.termsUsed > sin100.termsUsed && sin100.termsUsed > sin16.termsUsed,
        "CertifiedTrig: requested precision increases the proven convergence work");

    const auto cos100 = approximation::approximateCos(RealNumber{BigInt{1}}, 100);
    tests.expectEqual(std::string{cos100.value.text()},
        std::string{"0.5403023058681397174009366074429766037323104206179222276700972553811003947744717645179518560871830893"},
        "CertifiedTrig: cos(1) at 100 fractional digits");

    const auto tan100 = approximation::approximateTan(RealNumber{BigInt{1}}, 100);
    tests.expectEqual(std::string{tan100.value.text()},
        std::string{"1.5574077246549022305069748074583601730872507723815200383839466056988613971517272895550999652022429838"},
        "CertifiedTrig: tan(1) is certified from sin/cos enclosures");

    const auto sinOneDegree = approximation::approximateSinTurns(
        numeric::Rational{BigInt{1}, BigInt{360}}, 100);
    tests.expectEqual(std::string{sinOneDegree.value.text()},
        std::string{"0.0174524064372835128194189785163161924722527203071396426836124276405973842039280700420019267910213469"},
        "CertifiedTrig: one degree uses certified Pi and turn reduction");

    const auto cosOneDegree = approximation::approximateCosTurns(
        numeric::Rational{BigInt{1}, BigInt{360}}, 100);
    tests.expectEqual(std::string{cosOneDegree.value.text()},
        std::string{"0.9998476951563912391570115588139148516927403105831859396583207145115391811033372153972993952881103455"},
        "CertifiedTrig: cosine of one degree is certified through the same turn path");

    const auto tanOneDegree = approximation::approximateTanTurns(
        numeric::Rational{BigInt{1}, BigInt{360}}, 100);
    tests.expectEqual(std::string{tanOneDegree.value.text()},
        std::string{"0.0174550649282175857651288952197278243141015888398752769047114271021048548564623676228896891582992038"},
        "CertifiedTrig: tangent of one degree uses exact turn reduction and certified Pi");

    const auto sin361Degrees = approximation::approximateSinTurns(
        numeric::Rational{BigInt{361}, BigInt{360}}, 50);
    tests.expectEqual(std::string{sin361Degrees.value.text()},
        std::string{"0.01745240643728351281941897851631619247225272030714"},
        "CertifiedTrig: turn reduction removes complete periods exactly");


    const auto hugeRadianSin = approximation::approximateSin(
        RealNumber{BigInt{1'000'000}}, 20);
    tests.expectEqual(std::string{hugeRadianSin.value.text()},
        std::string{"-0.34999350217129295212"},
        "CertifiedTrig: huge raw radian uses certified Pi/2 argument reduction");

    const auto hugeRadianCos = approximation::approximateCos(
        RealNumber{BigInt{1'000'000}}, 20);
    tests.expectEqual(std::string{hugeRadianCos.value.text()},
        std::string{"0.93675212753314478694"},
        "CertifiedTrig: huge raw radian preserves quadrant mapping");

    // ComplexIntervalはRealIntervalを直交座標へ拡張したcertified complex domain。
    // 実部・虚部それぞれが真値を含むことを、exact Rationalで直接検査する。
    const auto realOne = approximation::RealInterval::fromRational(
        numeric::Rational{BigInt{1}}, 32);
    const auto realTwo = approximation::RealInterval::fromRational(
        numeric::Rational{BigInt{2}}, 32);
    const auto realThree = approximation::RealInterval::fromRational(
        numeric::Rational{BigInt{3}}, 32);
    const auto realFour = approximation::RealInterval::fromRational(
        numeric::Rational{BigInt{4}}, 32);

    const approximation::ComplexInterval z1{realOne, realTwo};
    const approximation::ComplexInterval z2{realThree, realFour};
    const auto product = approximation::multiply(z1, z2, 32);
    tests.expect(product.real().contains(numeric::Rational{BigInt{-5}}),
        "ComplexInterval: multiplication encloses exact real component");
    tests.expect(product.imaginary().contains(numeric::Rational{BigInt{10}}),
        "ComplexInterval: multiplication encloses exact imaginary component");

    const auto quotient = approximation::divide(z1, z2, 48);
    tests.expect(quotient.real().contains(numeric::Rational{BigInt{11}, BigInt{25}}),
        "ComplexInterval: division encloses exact real component");
    tests.expect(quotient.imaginary().contains(numeric::Rational{BigInt{2}, BigInt{25}}),
        "ComplexInterval: division encloses exact imaginary component");
    tests.expect(approximation::ComplexInterval::fromReal(realOne).isProvablyReal(),
        "ComplexInterval: only exact zero imaginary interval proves a real value");

    // CertifiedEvaluatorは再帰評価器なので、病的に深いASTをOSのstack overflowへ
    // 到達させずunsupportedとして返す。Windows Debugの小さいstackでも安全に失敗する。
    symbols::SymbolTable symbols;
    const auto builtins = evaluation::BuiltinRegistry::defaults(symbols);
    const auto mathematics = mathematics::MathRegistry::defaults(symbols, builtins);
    const mathematics::AngleSemantics angles;
    approximation::CertifiedEvaluator evaluator{builtins, mathematics, angles};
    expression::Expr deep{numeric::Number{BigInt{1}}};
    for (std::size_t i = 0; i < 128; ++i)
        deep = expression::Expr::call(
            builtins.symbol(evaluation::BuiltinId::Sin), {std::move(deep)});
    tests.expect(!evaluator.enclose(deep, 80),
        "CertifiedEvaluator: rejects pathological AST depth before recursive enclosure");
}

} // namespace mmcal::tests
