// 実数・複素数を統合するNumberの回帰テスト
#include "number_tests.hpp"

#include "numeric/number.hpp"
#include "test_framework.hpp"

#include <stdexcept>

namespace mmcal::tests {

void runNumberTests(TestRunner& tests) {
    using numeric::BigInt;
    using numeric::Number;
    using numeric::Rational;
    using numeric::RealNumber;

    const Number real{BigInt{5}};
    tests.expect(real.isReal(), "Number: stores real number");
    tests.expectEqual(real.toString(), "5", "Number: formats real number");

    const Number collapsed = Number::complex(RealNumber{BigInt{5}}, RealNumber{});
    tests.expect(collapsed.isReal(), "Number: reduces zero-imaginary complex to real");

    const Number imaginaryUnit = Number::complex(RealNumber{}, RealNumber{BigInt{1}});
    tests.expectEqual(imaginaryUnit.toString(), "I", "Number: formats imaginary unit");
    tests.expectEqual((-imaginaryUnit).toString(), "-I", "Number: formats negative imaginary unit");

    const Number z = Number::complex(RealNumber{BigInt{3}}, RealNumber{BigInt{-2}});
    tests.expectEqual(z.toString(), "3-2I", "Number: formats complex number");
    tests.expectEqual(z.conjugate().toString(), "3+2I", "Number: complex conjugate");

    const Number lhs = Number::complex(RealNumber{BigInt{1}}, RealNumber{BigInt{2}});
    const Number rhs = Number::complex(RealNumber{BigInt{3}}, RealNumber{BigInt{-4}});

    tests.expectEqual((lhs + rhs).toString(), "4-2I", "Number: complex addition");
    tests.expectEqual((lhs - rhs).toString(), "-2+6I", "Number: complex subtraction");
    tests.expectEqual((lhs * rhs).toString(), "11+2I", "Number: complex multiplication");
    tests.expectEqual((lhs / rhs).toString(), "-1/5+2I/5", "Number: complex division");

    auto cancelImaginary = Number::complex(RealNumber{BigInt{2}}, RealNumber{BigInt{3}});
    cancelImaginary += Number::complex(RealNumber{}, RealNumber{BigInt{-3}});
    tests.expect(cancelImaginary.isReal() && cancelImaginary == Number{BigInt{2}},
        "Number: arithmetic normalizes a zero imaginary component back to real storage");

    auto scaledImaginary = Number::complex(RealNumber{BigInt{4}}, RealNumber{});
    scaledImaginary *= Number{BigInt{2}};
    tests.expect(scaledImaginary.isReal() && scaledImaginary == Number{BigInt{8}},
        "Number: real scaling preserves the zero-imaginary normalization invariant");

    auto selfMultiply = lhs;
    selfMultiply *= selfMultiply;
    tests.expectEqual(selfMultiply.toString(), "-3+4I",
        "Number: optimized complex multiplication is safe under self aliasing");

    auto selfDivide = lhs;
    selfDivide /= selfDivide;
    tests.expectEqual(selfDivide.toString(), "1",
        "Number: optimized complex division is safe under self aliasing");

    auto realFastPath = Number{Rational{BigInt{2}, BigInt{3}}};
    realFastPath += Number{Rational{BigInt{1}, BigInt{6}}};
    realFastPath *= Number{Rational{BigInt{6}, BigInt{5}}};
    tests.expectEqual(realFastPath.toString(), "1",
        "Number: real-real fast path preserves exact rational arithmetic");

    const Number exact = Number::complex(
        RealNumber{Rational{BigInt{1}, BigInt{3}}},
        RealNumber{Rational{BigInt{2}, BigInt{5}}});
    tests.expectEqual(exact.toString(), "1/3+2I/5", "Number: complex number with exact rational parts");

    tests.expect(
        Number{BigInt{2}} == Number{Rational{BigInt{4}, BigInt{2}}},
        "Number: compares integer with equivalent rational");
    tests.expect(
        Number::complex(RealNumber{BigInt{2}}, RealNumber{}) == Number{BigInt{2}},
        "Number: compares zero-imaginary complex with real");

    tests.expectThrows<std::logic_error>(
        [&] { static_cast<void>(z.asReal()); },
        "Number: rejects complex real access");
    tests.expectThrows<std::domain_error>(
        [&] { static_cast<void>(lhs / Number{}); },
        "Number: rejects complex division by zero");
}

} // namespace mmcal::tests
