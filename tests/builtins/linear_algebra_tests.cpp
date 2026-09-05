// 行列・線形代数の回帰テスト
#include "linear_algebra_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "linear_algebra/fraction_free_elimination.hpp"
#include "linear_algebra/modular_linear_algebra.hpp"
#include "numeric/big_int.hpp"
#include "numeric/rational.hpp"
#include "test_framework.hpp"

#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::tests {
namespace {

std::string eval(kernel::KernelSession& session, std::string_view source) {
    return formatting::formatExpr(session.evaluate(source));
}

error::CalcError evalError(kernel::KernelSession& session, std::string_view source) {
    try { static_cast<void>(session.evaluate(source)); }
    catch (const error::CalcError& e) { return e; }
    throw std::logic_error("Expected CalcError");
}

linear_algebra::IntegerMatrixBuffer denseIntegerMatrix(
    std::size_t rows,
    std::size_t columns,
    std::size_t coefficientBits) {
    numeric::BigInt value{1};
    if (coefficientBits > 1)
        value <<= coefficientBits - 1;
    std::vector<numeric::BigInt> elements(rows * columns, value);
    return linear_algebra::IntegerMatrixBuffer{rows, columns, std::move(elements)};
}

std::string denseShiftedIdentity(std::size_t order) {
    std::string result{"{"};
    for (std::size_t row = 0; row < order; ++row) {
        if (row != 0)
            result += ',';
        result += '{';
        for (std::size_t column = 0; column < order; ++column) {
            if (column != 0)
                result += ',';
            result += row == column ? "2" : "1";
        }
        result += '}';
    }
    result += '}';
    return result;
}

std::string repeatedIntegerVector(std::size_t count, std::size_t value) {
    std::string result{"{"};
    for (std::size_t index = 0; index < count; ++index) {
        if (index != 0)
            result += ", ";
        result += std::to_string(value);
    }
    result += '}';
    return result;
}

} // namespace

void runLinearAlgebraTests(TestRunner& tests) {
    using numeric::BigInt;
    using numeric::Rational;

    // modular backendはword-sized imageからCRTで整数を一意復元し，
    // solveではrational reconstruction後に元の整数系でexact verificationする。
    {
        linear_algebra::IntegerMatrixBuffer matrix{3, 3, {
            BigInt::parse("2147483659"), BigInt{2}, BigInt{-3},
            BigInt{5}, BigInt::parse("4294967311"), BigInt{7},
            BigInt{-11}, BigInt{13}, BigInt::parse("8589934609")}};
        const BigInt expected = linear_algebra::bareissDeterminant(matrix);
        linear_algebra::ModularLinearAlgebraStats stats;
        const BigInt actual = linear_algebra::modularDeterminant(matrix, &stats);
        tests.expect(actual == expected && stats.primesAccepted >= 2,
            "Linear algebra: modular determinant reconstructs the exact BigInt through CRT");
    }
    {
        linear_algebra::IntegerMatrixBuffer augmented{3, 4, {
            BigInt{2}, BigInt{}, BigInt{}, BigInt{1},
            BigInt{}, BigInt{3}, BigInt{}, BigInt{1},
            BigInt{}, BigInt{}, BigInt{5}, BigInt{1}}};
        const auto solution = linear_algebra::modularSolve(augmented, 3);
        tests.expect(solution && *solution == std::vector<Rational>{
            Rational{BigInt{1}, BigInt{2}},
            Rational{BigInt{1}, BigInt{3}},
            Rational{BigInt{1}, BigInt{5}}},
            "Linear algebra: modular solve reconstructs and verifies rational solutions");
    }
    {
        constexpr std::int64_t firstPrime = 2'147'483'647LL;
        linear_algebra::IntegerMatrixBuffer augmented{2, 3, {
            BigInt{firstPrime}, BigInt{}, BigInt{firstPrime},
            BigInt{}, BigInt{1}, BigInt{2}}};
        linear_algebra::ModularLinearAlgebraStats stats;
        const auto solution = linear_algebra::modularSolve(augmented, 2, &stats);
        tests.expect(solution && *solution == std::vector<Rational>{
            Rational{BigInt{1}}, Rational{BigInt{2}}}
            && stats.primesTried > stats.primesAccepted,
            "Linear algebra: modular solve skips bad primes without changing the exact solution");
    }
    {
        linear_algebra::IntegerMatrixBuffer matrix{3, 3, {
            BigInt{2}, BigInt{}, BigInt{},
            BigInt{}, BigInt{3}, BigInt{},
            BigInt{}, BigInt{}, BigInt{5}}};
        const auto inverse = linear_algebra::modularInverse(matrix);
        tests.expect(inverse && *inverse == std::vector<Rational>{
            Rational{BigInt{1}, BigInt{2}}, Rational{}, Rational{},
            Rational{}, Rational{BigInt{1}, BigInt{3}}, Rational{},
            Rational{}, Rational{}, Rational{BigInt{1}, BigInt{5}}},
            "Linear algebra: modular inverse reconstructs the adjugate and exact common determinant");
    }

    kernel::KernelSession session;
    tests.expectEqual(eval(session, "dimensions[{{1,2},{3}}]"), std::string{"{2}"},
        "Array: non-rectangular brace values expose only their common dimensions");
    tests.expectEqual(eval(session, "arrayRank[{{1,2},{3}}]"), std::string{"1"},
        "Array: non-rectangular brace rank follows the common rectangular prefix");
    tests.expectEqual(eval(session, "length[{{1,2},{3}}]"), std::string{"2"},
        "Array: length is defined for a general brace value");
    tests.expectEqual(eval(session, "at[{{1,2},{3}},0]"), std::string{"{1, 2}"},
        "Array: at indexes a general brace value");
    tests.expectEqual(eval(session, "transpose[{{1,2,3},{4,5,6}}]"),
        std::string{"{{1, 4}, {2, 5}, {3, 6}}"},
        "Linear algebra: transpose preserves exact elements");
    tests.expectEqual(eval(session, "dimensions[transpose[zeros[0,3]]]"),
        std::string{"{3, 0}"},
        "Linear algebra: transpose handles zero-sized matrices without dimension loops");
    tests.expectEqual(eval(session, "madd[{{1,2},{3,4}},{{5,6},{7,8}}]"),
        std::string{"{{6, 8}, {10, 12}}"},
        "Linear algebra: matrix addition is exact");
    tests.expectEqual(eval(session, "dot[{{1,2},{3,4}},{{5,6},{7,8}}]"),
        std::string{"{{19, 22}, {43, 50}}"},
        "Linear algebra: matrix multiplication is exact through dot");
    tests.expectEqual(eval(session, "dot[{1,2,3},{4,5,6}]"), std::string{"32"},
        "Linear algebra: vector dot product uses dot");
    tests.expectEqual(eval(session, "dot[{{1,2},{3,4}},{5,6}]"),
        std::string{"{17, 39}"},
        "Linear algebra: matrix-vector dot preserves rank");
    tests.expectEqual(eval(session, "dot[{5,6},{{1,2},{3,4}}]"),
        std::string{"{23, 34}"},
        "Linear algebra: vector-matrix dot preserves rank");
    tests.expectEqual(eval(session, "matmul[{{1,2},{3,4}},{{5,6},{7,8}}]"),
        std::string{"{{19, 22}, {43, 50}}"},
        "Linear algebra: legacy matmul aliases dot");
    tests.expectEqual(eval(session, "det[{{1,2},{3,4}}]"), std::string{"-2"},
        "Linear algebra: determinant is exact");
    tests.expectEqual(eval(session, "det[{{1/2,1/3},{2/5,3/7}}]"),
        std::string{"17/210"},
        "Linear algebra: Bareiss clears rational row denominators exactly");
    tests.expectEqual(eval(session, "det[{{0,2,1},{3,4,5},{6,7,8}}]"),
        std::string{"9"},
        "Linear algebra: Bareiss row pivoting preserves determinant sign");

    kernel::KernelSession modularSession;
    const std::string dense48 = denseShiftedIdentity(48);
    tests.expectEqual(eval(modularSession, "det[" + dense48 + "]"), std::string{"49"},
        "Linear algebra: automatic determinant dispatcher preserves the exact dense determinant");
    tests.expect(modularSession.lastEvaluationUsage().modularPrimes > 0,
        "Linear algebra: automatic determinant dispatcher selects CRT past the measured crossover");

    const std::string dense32 = denseShiftedIdentity(32);
    tests.expectEqual(eval(modularSession,
        "solveLinear[" + dense32 + "," + repeatedIntegerVector(32, 33) + "]"),
        repeatedIntegerVector(32, 1),
        "Linear algebra: automatic solve dispatcher preserves exact dense solutions");
    tests.expect(modularSession.lastEvaluationUsage().modularPrimes > 0,
        "Linear algebra: automatic solve dispatcher selects CRT/rational reconstruction");

    const auto conservativeInverse = modularSession.evaluate("inverse[identity[10]]");
    tests.expect(conservativeInverse.isArray()
        && conservativeInverse.asArray().shape == std::vector<std::size_t>{10, 10},
        "Linear algebra: automatic inverse dispatcher returns the exact matrix shape");
    tests.expect(modularSession.lastEvaluationUsage().modularPrimes == 0,
        "Linear algebra: inverse dispatcher keeps Bareiss while modular inverse has no measured crossover");

    tests.expect(!linear_algebra::preferModularDeterminant(denseIntegerMatrix(47, 47, 1)),
        "Linear algebra: determinant dispatcher stays Bareiss one order below the size-only threshold");
    tests.expect(linear_algebra::preferModularDeterminant(denseIntegerMatrix(48, 48, 1)),
        "Linear algebra: determinant dispatcher selects modular at the size-only threshold");
    tests.expect(!linear_algebra::preferModularDeterminant(denseIntegerMatrix(32, 32, 63)),
        "Linear algebra: determinant dispatcher stays Bareiss one bit below the 32x32 height threshold");
    tests.expect(linear_algebra::preferModularDeterminant(denseIntegerMatrix(32, 32, 64)),
        "Linear algebra: determinant dispatcher selects modular at the 32x32 height threshold");
    tests.expect(!linear_algebra::preferModularDeterminant(denseIntegerMatrix(24, 24, 191)),
        "Linear algebra: determinant dispatcher stays Bareiss one bit below the 24x24 height threshold");
    tests.expect(linear_algebra::preferModularDeterminant(denseIntegerMatrix(24, 24, 192)),
        "Linear algebra: determinant dispatcher selects modular at the 24x24 height threshold");

    tests.expect(!linear_algebra::preferModularSolve(denseIntegerMatrix(23, 24, 1), 23),
        "Linear algebra: solve dispatcher stays Bareiss one variable below the size-only threshold");
    tests.expect(linear_algebra::preferModularSolve(denseIntegerMatrix(24, 25, 1), 24),
        "Linear algebra: solve dispatcher selects modular at the size-only threshold");
    tests.expect(!linear_algebra::preferModularSolve(denseIntegerMatrix(12, 13, 95), 12),
        "Linear algebra: solve dispatcher stays Bareiss one bit below the 12-variable height threshold");
    tests.expect(linear_algebra::preferModularSolve(denseIntegerMatrix(12, 13, 96), 12),
        "Linear algebra: solve dispatcher selects modular at the 12-variable height threshold");
    tests.expect(!linear_algebra::preferModularSolve(denseIntegerMatrix(8, 9, 255), 8),
        "Linear algebra: solve dispatcher stays Bareiss one bit below the 8-variable height threshold");
    tests.expect(linear_algebra::preferModularSolve(denseIntegerMatrix(8, 9, 256), 8),
        "Linear algebra: solve dispatcher selects modular at the 8-variable height threshold");
    tests.expect(!linear_algebra::preferModularSolve(denseIntegerMatrix(6, 7, 511), 6),
        "Linear algebra: solve dispatcher stays Bareiss one bit below the 6-variable height threshold");
    tests.expect(linear_algebra::preferModularSolve(denseIntegerMatrix(6, 7, 512), 6),
        "Linear algebra: solve dispatcher selects modular at the 6-variable height threshold");
    tests.expectEqual(eval(session, "det[{{a,b},{c,d}}]"), std::string{"a d-b c"},
        "Linear algebra: symbolic determinant remains exact");
    tests.expect(eval(session,
        "det[{{a,b,c,d,e,f},{g,h,i,j,k,l},{m,n,o,p,q,r},{s,t,u,v,w,x},"
        "{y,z,aa,bb,cc,dd},{ee,ff,gg,hh,ii,jj}}]").starts_with("det["),
        "Linear algebra: dense symbolic determinant stops before factorial expression growth");
    tests.expectEqual(eval(session,
        "det[{{a,b,c,d,e},{0,f,g,h,i},{0,0,j,k,l},{0,0,0,m,n},{0,0,0,0,p}}]"),
        std::string{"a f j m p"},
        "Linear algebra: triangular symbolic determinant bypasses the expansion budget");
    tests.expectEqual(eval(session, "inverse[{{1,2},{3,4}}]"),
        std::string{"{{-2, 1}, {3/2, -1/2}}"},
        "Linear algebra: exact rational matrix inverse");
    tests.expectEqual(eval(session, "inverse[{{1/2,1/3},{2/5,3/7}}]"),
        std::string{"{{90/17, -70/17}, {-84/17, 105/17}}"},
        "Linear algebra: fraction-free augmented elimination preserves rational inverse");
    tests.expectEqual(eval(session, "inverse[{{0,2},{3,4}}]"),
        std::string{"{{-2/3, 1/3}, {1/2, 0}}"},
        "Linear algebra: BigInt Bareiss back substitution preserves row-swapped inverse signs");
    tests.expectEqual(eval(session, "inverse[{{a,b},{c,d}}]"),
        std::string{"{{d/(a d-b c), -b/(a d-b c)}, {-c/(a d-b c), a/(a d-b c)}}"},
        "Linear algebra: symbolic inverse keeps determinant denominators");
    tests.expectEqual(eval(session, "rref[{{1,2},{3,4}}]"),
        std::string{"{{1, 0}, {0, 1}}"},
        "Linear algebra: rref uses exact pivots");
    tests.expectEqual(eval(session, "rref[{{1,0},{0,1},{1,1}}]"),
        std::string{"{{1, 0}, {0, 1}, {0, 0}}"},
        "Linear algebra: full-column-rank RREF skips Rational backward elimination");
    tests.expectEqual(eval(session, "rref[{{1/2,1/3,5/6},{1,2/3,5/3}}]"),
        std::string{"{{1, 2/3, 5/3}, {0, 0, 0}}"},
        "Linear algebra: Bareiss forward elimination feeds canonical rational RREF");
    tests.expectEqual(eval(session, "matrixRank[{{0,1,2},{0,2,4},{3,0,0}}]"),
        std::string{"2"},
        "Linear algebra: Bareiss rank handles skipped columns and row pivots");
    tests.expectEqual(eval(session, "matrixRank[{{1,2},{2,4}}]"), std::string{"1"},
        "Linear algebra: matrixRank is exact for rational matrices");
    tests.expectEqual(eval(session, "solveLinear[{{2,1},{1,-1}},{5,1}]"),
        std::string{"{2, 1}"},
        "Linear algebra: solveLinear returns an exact unique solution");
    tests.expectEqual(eval(session,
        "solveLinear[{{1/2,1/3},{2/5,3/7}},{0,-17/35}]"),
        std::string{"{2, -3}"},
        "Linear algebra: solveLinear clears rational row denominators with Bareiss");
    tests.expectEqual(eval(session,
        "solveLinear[{{1,0},{0,1},{1,1}},{2,3,5}]"),
        std::string{"{2, 3}"},
        "Linear algebra: solveLinear accepts consistent overdetermined systems");
    tests.expectEqual(eval(session,
        "solveLinear[{{1,I},{I,1}},{1+2I,2+I}]"),
        std::string{"{1, 2}"},
        "Linear algebra: solveLinear keeps exact complex Gaussian fallback");
    tests.expectEqual(eval(session, "solveLinear[{{1,0},{0,1}},{x,y}]"),
        std::string{"{x, y}"},
        "Linear algebra: solveLinear handles symbolic right-hand sides with decidable pivots");
    tests.expectEqual(eval(session, "solveLinear[reshape[{}, {0,0}],{}]"),
        std::string{"{}"},
        "Linear algebra: solveLinear preserves the unique empty solution for a 0x0 system");
    tests.expectEqual(eval(session, "solveLinear[{{x,0},{0,1}},{1,2}]"),
        std::string{"solveLinear[{{x, 0}, {0, 1}}, {1, 2}]"},
        "Linear algebra: solveLinear leaves undecidable symbolic pivots unevaluated");
    tests.expectEqual(eval(session, "rank[{{1,2},{2,4}}]"), std::string{"1"},
        "Linear algebra: legacy rank aliases matrixRank");
    tests.expectEqual(eval(session, "matrixRank[{{x,0},{0,1}}]"),
        std::string{"matrixRank[{{x, 0}, {0, 1}}]"},
        "Linear algebra: symbolic matrixRank does not guess an undecidable pivot");
    tests.expectEqual(eval(session, "det[{{1+I,2},{3,4-I}}]"),
        std::string{"-1+3I"},
        "Linear algebra: exact complex matrices keep the Gaussian fallback");
    tests.expectEqual(eval(session, "trace[{{1,2},{3,4}}]"), std::string{"5"},
        "Linear algebra: trace uses the shared exact matrix path");
    tests.expectEqual(eval(session, "norm[{3+4I}]"), std::string{"5"},
        "Linear algebra: norm is Hermitian for complex vectors");
    tests.expectEqual(eval(session, "normalize[{3,4}]"), std::string{"{3/5, 4/5}"},
        "Linear algebra: normalize preserves exact rationals");
    tests.expectEqual(eval(session, "nullSpace[{{1,2},{2,4}}]"),
        std::string{"{{-2, 1}}"},
        "Linear algebra: nullSpace returns the canonical RREF basis");
    tests.expectEqual(eval(session, "nullSpace[{{1,2,3},{2,4,6}}]"),
        std::string{"{{-2, 1, 0}, {-3, 0, 1}}"},
        "Linear algebra: nullSpace handles rectangular rank-deficient matrices");
    tests.expectEqual(eval(session, "nullSpace[{{1,0},{0,1}}]"),
        std::string{"reshape[{}, {0, 2}]"},
        "Linear algebra: full-column-rank nullSpace preserves the empty basis shape");
    tests.expectEqual(eval(session, "nullSpace[reshape[{}, {0,3}]]"),
        std::string{"{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}"},
        "Linear algebra: zero-row matrix nullSpace is the full coordinate space");
    tests.expectEqual(eval(session, "nullSpace[{{1,I},{I,-1}}]"),
        std::string{"{{-I, 1}}"},
        "Linear algebra: nullSpace keeps the exact complex Gaussian fallback");
    tests.expectEqual(eval(session, "nullSpace[{{1,a},{0,0}}]"),
        std::string{"{{-a, 1}}"},
        "Linear algebra: nullSpace handles symbolic free-variable coefficients with a decidable pivot");
    tests.expectEqual(eval(session, "nullSpace[{{x,1}}]"),
        std::string{"nullSpace[{{x, 1}}]"},
        "Linear algebra: nullSpace does not guess an undecidable symbolic pivot");
    tests.expectEqual(eval(session, "N[nullSpace[{{1,Pi}}],12]"),
        std::string{"{{-3.14159265359, 1}}"},
        "Linear algebra: N keeps exact-first pivot structure for nullSpace before approximating its basis");
    tests.expectEqual(eval(session, "nullSpace[N[{{Pi,1}},12]]"),
        std::string{"{{-0.31830988618, 1.0}}"},
        "Linear algebra: approximate nullSpace accepts a pivot structure certified by input information");
    tests.expectEqual(eval(session, "nullSpace[N[{{Pi,0},{0,0}},12]]"),
        std::string{"nullSpace[{{3.14159265359, 0.0}, {0.0, 0.0}}]"},
        "Linear algebra: approximate nullSpace does not recover nullity from finite-precision zero truth");
    tests.expectEqual(eval(session, "luDecomposition[{{0,2},{3,4}}]"),
        std::string{"{{{0, 1}, {1, 0}}, {{1, 0}, {0, 1}}, {{3, 4}, {0, 2}}}"},
        "Linear algebra: LU returns {P,L,U} with PA=LU and exact row pivoting");
    tests.expectEqual(eval(session,
        "dot[at[luDecomposition[{{0,2},{3,4}}],0],{{0,2},{3,4}}]"),
        eval(session,
            "dot[at[luDecomposition[{{0,2},{3,4}}],1],at[luDecomposition[{{0,2},{3,4}}],2]]"),
        "Linear algebra: LU factors reconstruct the permuted exact matrix");
    tests.expectEqual(eval(session, "qrDecomposition[{{3,0},{4,0}}]"),
        std::string{"{{{3/5, 4/5}, {4/5, -3/5}}, {{5, 0}, {0, 0}}}"},
        "Linear algebra: fraction-free exact QR delays normalization and preserves rational factors when the norm closes exactly");
    tests.expectEqual(eval(session,
        "N[dot[at[qrDecomposition[{{1,2},{3,4}}],0],at[qrDecomposition[{{1,2},{3,4}}],1]],8]"),
        std::string{"{{1.0, 2.0}, {3.0, 4.0}}"},
        "Linear algebra: fraction-free exact QR reconstructs A");
    tests.expectEqual(eval(session, "N[qrDecomposition[{{1,2},{3,4}}],8]"),
        std::string{"{{{-0.31622777, -0.94868330}, {-0.94868330, 0.31622777}}, {{-3.1622777, -4.4271887}, {0, -0.63245553}}}"},
        "Linear algebra: N dispatches QR directly to the certified Householder backend");
    tests.expectEqual(eval(session, "dimensions[at[qrDecomposition[{{1,0},{0,1},{0,0}}],0]]"),
        std::string{"{3, 2}"},
        "Linear algebra: reduced QR exposes a rectangular Q through {Q,R}");
    tests.expectEqual(eval(session, "dimensions[at[qrDecomposition[{{1,0},{0,1},{0,0}}],1]]"),
        std::string{"{2, 2}"},
        "Linear algebra: reduced QR exposes a rectangular R through {Q,R}");
    tests.expectEqual(eval(session,
        "N[dot[at[qrDecomposition[{{1,2},{3,4},{5,6}}],0],at[qrDecomposition[{{1,2},{3,4},{5,6}}],1]],8]"),
        std::string{"{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}"},
        "Linear algebra: rectangular reduced QR reconstructs A");
    tests.expectEqual(eval(session, "N[luDecomposition[{{Pi,1},{2,3}}],8]"),
        std::string{"{{{1, 0}, {0, 1}}, {{1, 0}, {0.63661977, 1}}, {{3.1415927, 1}, {0, 2.3633802}}}"},
        "Linear algebra: N dispatches LU directly to the certified decomposition backend");
    tests.expectEqual(eval(session, "qrDecomposition[{{a,b},{0,c}}]"),
        std::string{"{{{1, 0}, {0, 1}}, {{a, b}, {0, c}}}"},
        "Linear algebra: already upper-triangular symbolic matrices avoid unnecessary QR expansion");
    tests.expectEqual(eval(session,
        "dimensions[at[qrDecomposition[{{19,-1,-2,2},{2,19,0,-1},{-1,-2,19,1},{1,0,-1,19}}],0]]"),
        std::string{"{4, 4}"},
        "Linear algebra: fraction-free exact QR no longer has the former 3x3 order cap");
    tests.expectEqual(eval(session,
        "N[dot[at[qrDecomposition[{{19,-1,-2,2},{2,19,0,-1},{-1,-2,19,1},{1,0,-1,19}}],0],at[qrDecomposition[{{19,-1,-2,2},{2,19,0,-1},{-1,-2,19,1},{1,0,-1,19}}],1]],8]"),
        std::string{"{{19.0, -1.0, -2.0, 2.0}, {2.0, 19.0, 0.0, -1.0}, {-1.0, -2.0, 19.0, 1.0}, {1.0, 0.0, -1.0, 19.0}}"},
        "Linear algebra: uncapped fraction-free exact QR reconstructs a dense 4x4 matrix");
    tests.expectEqual(eval(session,
        "dimensions[at[qrDecomposition[" + denseShiftedIdentity(8) + "],0]]"),
        std::string{"{8, 8}"},
        "Linear algebra: fraction-free exact QR remains uncapped beyond the former small-order regime");
    tests.expectEqual(eval(session,
        "N[dot[at[qrDecomposition[{{1/2,1/3},{2/5,3/7}}],0],at[qrDecomposition[{{1/2,1/3},{2/5,3/7}}],1]],8]"),
        std::string{"{{0.50, 0.33333333}, {0.40, 0.42857143}}"},
        "Linear algebra: fraction-free exact QR preserves Rational row-denominator semantics");
    tests.expectEqual(eval(session, "N[det[{{Pi,0},{0,2}}],12]"),
        std::string{"6.28318530718"},
        "Linear algebra: N pushes precision directly into determinant");
    tests.expectEqual(eval(session, "N[dot[{{Pi,0},{0,Pi}},{{1,2},{3,4}}],12]"),
        std::string{"{{3.14159265359, 6.28318530718}, {9.42477796077, 12.5663706144}}"},
        "Linear algebra: N pushes precision directly into dot");
    tests.expectEqual(eval(session, "N[inverse[{{Pi,0},{0,2}}],12]"),
        std::string{"{{0.318309886184, 0.0}, {0.0, 0.50}}"},
        "Linear algebra: N uses the certified inverse backend");
    tests.expectEqual(eval(session, "N[rref[{{Pi,0},{0,2}}],12]"),
        std::string{"{{1.0, 0.0}, {0.0, 1.0}}"},
        "Linear algebra: N uses certified pivots for rref");
    tests.expectEqual(eval(session, "N[norm[{Pi,1}],12]"),
        std::string{"3.29690830948"},
        "Linear algebra: N uses the certified Hermitian norm backend");
    tests.expectEqual(eval(session, "dot[N[{{Pi,0},{0,Pi}},12],{{1,2},{3,4}}]"),
        std::string{"{{3.14159265359, 6.28318530718}, {9.42477796077, 12.5663706144}}"},
        "Linear algebra: approximate Array inputs infer backend precision like FFT");
    tests.expectEqual(eval(session, "N[matrixRank[{{Pi,2Pi},{1,2}}],12]"),
        std::string{"1"},
        "Linear algebra: N displays an exact-source integer rank without a redundant decimal marker");
    tests.expectEqual(eval(session, "matrixRank[N[{{Pi,2Pi},{1,2}},12]]"),
        std::string{"matrixRank[{{3.14159265359, 6.28318530718}, {1.0, 2.0}}]"},
        "Linear algebra: approximate matrixRank does not invent a rank-deficiency threshold");
    const std::string finiteNearSingular =
        "{{N[1,5],N[1,5]},{N[1,5],N[1+1/10^10,5]}}";
    tests.expectEqual(eval(session, "matrixRank[" + finiteNearSingular + "]"),
        std::string{"matrixRank[{{1.0, 1.0}, {1.0, 1.0}}]"},
        "Linear algebra: finite input information cannot prove a hidden near-singular rank");
    tests.expectEqual(eval(session, "inverse[" + finiteNearSingular + "]"),
        std::string{"inverse[{{1.0, 1.0}, {1.0, 1.0}}]"},
        "Linear algebra: inverse does not use hidden certified points to prove pivots");
    tests.expectEqual(eval(session, "precision[det[" + finiteNearSingular + "]]"),
        std::string{"0"},
        "Linear algebra: determinant cancellation does not invent relative precision");
    tests.expectEqual(eval(session, "accuracy[det[" + finiteNearSingular + "]]"),
        std::string{"3"},
        "Linear algebra: determinant cancellation retains only input-supported absolute accuracy");
    tests.expectEqual(eval(session, "N[conditionNumber[" + finiteNearSingular + "],20]"),
        std::string{"conditionNumber[{{1.0, 1.0}, {1.0, 1.0}}]"},
        "Linear algebra: conditionNumber does not infer hidden rank from finite input information");
    tests.expectEqual(eval(session, "N[svd[" + finiteNearSingular + "],20]"),
        std::string{"svd[{{1.0, 1.0}, {1.0, 1.0}}]"},
        "Linear algebra: SVD stays conservative for explicit finite-precision matrices");
    tests.expectEqual(eval(session,
        "N[solveLinear[{{Pi,0},{0,2}},{Pi,4}],12]"),
        std::string{"{1.0, 2.0}"},
        "Linear algebra: N pushes precision directly into solveLinear");
    tests.expectEqual(eval(session, "conjugateTranspose[{{1,I},{2,3I}}]"),
        std::string{"{{1, 2}, {-I, -3I}}"},
        "Linear algebra: conjugateTranspose supplies the Hermitian transpose for complex factors");
    tests.expectEqual(eval(session, "svd[{{3,0},{0,4}}]"),
        std::string{"{{{0, 1}, {1, 0}}, {{4, 0}, {0, 3}}, {{0, 1}, {1, 0}}}"},
        "Linear algebra: exact diagonal SVD sorts singular values without numerical fallback");
    tests.expectEqual(eval(session, "dimensions[at[N[svd[{{1,2},{3,4},{5,6}}],8],0]]"),
        std::string{"{3, 2}"},
        "Linear algebra: reduced SVD returns U with shape m x min(m,n)");
    tests.expectEqual(eval(session, "dimensions[at[N[svd[{{1,2},{3,4},{5,6}}],8],2]]"),
        std::string{"{2, 2}"},
        "Linear algebra: reduced SVD returns V with shape n x min(m,n)");
    tests.expectEqual(eval(session,
        "N[dot[at[svd[{{1,2},{3,4},{5,6}}],0],dot[at[svd[{{1,2},{3,4},{5,6}}],1],transpose[at[svd[{{1,2},{3,4},{5,6}}],2]]]],8]"),
        std::string{"{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}}"},
        "Linear algebra: numerical reduced SVD reconstructs a tall matrix");
    tests.expectEqual(eval(session,
        "N[dot[at[svd[{{1,2,3},{4,5,6}}],0],dot[at[svd[{{1,2,3},{4,5,6}}],1],transpose[at[svd[{{1,2,3},{4,5,6}}],2]]]],8]"),
        std::string{"{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}}"},
        "Linear algebra: numerical reduced SVD reconstructs a wide matrix");
    tests.expectEqual(eval(session,
        "N[dot[at[svd[{{1,I},{0,1}}],0],dot[at[svd[{{1,I},{0,1}}],1],conjugateTranspose[at[svd[{{1,I},{0,1}}],2]]]],8]"),
        std::string{"{{1.0, 1.0I}, {0.0I, 1.0}}"},
        "Linear algebra: complex numerical SVD reconstructs A with V Hermitian-transposed");
    tests.expectEqual(eval(session, "conditionNumber[{{3,0},{0,4}}]"),
        std::string{"4/3"},
        "Linear algebra: conditionNumber preserves exact diagonal singular values");
    tests.expectEqual(eval(session, "conditionNumber[{{1,2},{2,4}}]"),
        std::string{"Infinity"},
        "Linear algebra: conditionNumber detects exact rank deficiency before numerical SVD");
    tests.expectEqual(eval(session, "conditionNumber[{{1,2,3}}]"),
        std::string{"1"},
        "Linear algebra: a nonzero one-singular-value rectangular matrix has condition number one");
    tests.expectEqual(eval(session, "pseudoInverse[{{1,2},{2,4}}]"),
        std::string{"{{1/25, 2/25}, {2/25, 4/25}}"},
        "Linear algebra: pseudoInverse keeps rank-deficient Rational matrices exact");
    tests.expectEqual(eval(session, "pseudoInverse[{{I,0},{0,2I}}]"),
        std::string{"{{-I, 0}, {0, -I/2}}"},
        "Linear algebra: pseudoInverse preserves exact complex arithmetic");
    tests.expectEqual(eval(session, "dimensions[pseudoInverse[zeros[0,3]]]"),
        std::string{"{3, 0}"},
        "Linear algebra: pseudoInverse transposes zero-sized matrix dimensions");
    tests.expectEqual(eval(session, "dimensions[pseudoInverse[zeros[3,0]]]"),
        std::string{"{0, 3}"},
        "Linear algebra: pseudoInverse handles zero-column matrices without storage assumptions");
    tests.expectEqual(eval(session,
        "leastSquares[{{1,0},{0,1},{1,1}},{1,2,4}]"),
        std::string{"{4/3, 7/3}"},
        "Linear algebra: leastSquares returns the exact minimum-norm solution through A^+ b");
    tests.expectEqual(eval(session, "leastSquares[{{1,2},{2,4}},{1,2}]"),
        std::string{"{1/5, 2/5}"},
        "Linear algebra: leastSquares handles exact rank deficiency without an epsilon threshold");
    tests.expectEqual(eval(session, "leastSquares[zeros[0,3],{}]"),
        std::string{"{0, 0, 0}"},
        "Linear algebra: leastSquares defines the empty-observation minimum-norm solution");
    tests.expectEqual(eval(session, "N[conditionNumber[{{1,2},{3,4}}],8]"),
        std::string{"14.933034"},
        "Linear algebra: N conditionNumber uses certified singular values");
    tests.expectEqual(eval(session, "N[pseudoInverse[{{1,2},{3,4}}],8]"),
        std::string{"{{-2.0, 1.0}, {1.50, -0.50}}"},
        "Linear algebra: N pseudoInverse uses the certified SVD backend");
    tests.expectEqual(eval(session,
        "N[leastSquares[{{1,0},{0,1},{1,1}},{1,2,4}],8]"),
        std::string{"{1.3333333, 2.3333333}"},
        "Linear algebra: N leastSquares propagates certified SVD precision");
    tests.expectEqual(eval(session, "eigenvalues[{{1,2},{3,4}}]"),
        std::string{"{(5+sqrt[33])/2, (5-sqrt[33])/2}"},
        "Linear algebra: exact 2x2 eigenvalues preserve radicals");
    tests.expectEqual(eval(session, "eigenvalues[{{0,-1},{1,0}}]"),
        std::string{"{I, -I}"},
        "Linear algebra: exact 2x2 eigenvalues extend naturally to complex roots");
    tests.expectEqual(eval(session, "eigenvalues[{{7/29,-4/29},{2/29,3/29}}]"),
        std::string{"{5/29+2I/29, 5/29-2I/29}"},
        "Linear algebra: exact rational imaginary coefficients are formatted unambiguously");
    tests.expectEqual(eval(session, "eigensystem[{{2,0},{0,3}}]"),
        std::string{"{{2, 3}, {{1, 0}, {0, 1}}}"},
        "Linear algebra: exact diagonal eigensystem returns column eigenvectors");
    tests.expectEqual(eval(session, "eigenvectors[{{1,2},{3,4}}]"),
        std::string{"{{2, 2}, {(5+sqrt[33])/2-1, (5-sqrt[33])/2-1}}"},
        "Linear algebra: distinct exact 2x2 eigenvectors stay symbolic and column-oriented");
    tests.expectEqual(eval(session, "eigenvectors[{{1,1},{0,1}}]"),
        std::string{"eigenvectors[{{1, 1}, {0, 1}}]"},
        "Linear algebra: defective exact 2x2 matrices do not receive duplicate eigenvectors");
    tests.expectEqual(eval(session, "N[eigenvalues[{{1,2},{3,4}}],8]"),
        std::string{"{-0.37228132, 5.3722813}"},
        "Linear algebra: N dispatches eigenvalues directly to the Schur backend");
    tests.expectEqual(eval(session, "N[eigenvalues[{{Pi,1},{0,2}}],8]"),
        std::string{"{3.1415927, 2.0}"},
        "Linear algebra: eigen relation audit uses certified input intervals");
    tests.expectEqual(eval(session, "N[eigenvectors[{{1,1},{0,1}}],8]"),
        std::string{"eigenvectors[{{1.0, 1.0}, {0.0, 1.0}}]"},
        "Linear algebra: defective repeated eigenvectors are not guessed");
    const auto inconsistent = evalError(session,
        "solveLinear[{{1,1},{2,2}},{1,3}]");
    tests.expect(inconsistent.type() == error::CalcErrorType::Domain,
        "Linear algebra: solveLinear rejects inconsistent systems");
    const auto nonUnique = evalError(session,
        "solveLinear[{{1,1},{2,2}},{1,2}]");
    tests.expect(nonUnique.type() == error::CalcErrorType::Domain,
        "Linear algebra: solveLinear rejects systems without a unique solution");
    const auto rhsMismatch = evalError(session,
        "solveLinear[{{1,0},{0,1}},{1,2,3}]");
    tests.expect(rhsMismatch.type() == error::CalcErrorType::Domain,
        "Linear algebra: solveLinear checks right-hand side dimensions");
    const auto singular = evalError(session, "inverse[{{1,2},{2,4}}]");
    tests.expect(singular.type() == error::CalcErrorType::Domain,
        "Linear algebra: singular inverse is a domain error");
    const auto mismatch = evalError(session, "dot[{{1,2}},{{1,2}}]");
    tests.expect(mismatch.type() == error::CalcErrorType::Domain,
        "Linear algebra: incompatible dimensions are rejected");
}

} // namespace mmcal::tests
