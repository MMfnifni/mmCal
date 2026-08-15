// 行列・線形代数の回帰テスト
#include "linear_algebra_tests.hpp"

#include "error/error_message.hpp"
#include "formatting/expr_formatter.hpp"
#include "kernel/kernel_session.hpp"
#include "test_framework.hpp"

#include <stdexcept>
#include <string>
#include <string_view>

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

} // namespace

void runLinearAlgebraTests(TestRunner& tests) {
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
    tests.expectEqual(eval(session, "inverse[{{a,b},{c,d}}]"),
        std::string{"{{d/(a d-b c), -b/(a d-b c)}, {-c/(a d-b c), a/(a d-b c)}}"},
        "Linear algebra: symbolic inverse keeps determinant denominators");
    tests.expectEqual(eval(session, "rref[{{1,2},{3,4}}]"),
        std::string{"{{1, 0}, {0, 1}}"},
        "Linear algebra: rref uses exact pivots");
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
    tests.expectEqual(eval(session, "nullSpace[N[{{Pi,0},{0,0}},12]]"),
        std::string{"{{0, 1}}"},
        "Linear algebra: approximate nullSpace accepts a certified pivot structure");
    tests.expectEqual(eval(session, "luDecomposition[{{0,2},{3,4}}]"),
        std::string{"{{{0, 1}, {1, 0}}, {{1, 0}, {0, 1}}, {{3, 4}, {0, 2}}}"},
        "Linear algebra: LU returns {P,L,U} with PA=LU and exact row pivoting");
    tests.expectEqual(eval(session,
        "dot[at[luDecomposition[{{0,2},{3,4}}],0],{{0,2},{3,4}}]"),
        eval(session,
            "dot[at[luDecomposition[{{0,2},{3,4}}],1],at[luDecomposition[{{0,2},{3,4}}],2]]"),
        "Linear algebra: LU factors reconstruct the permuted exact matrix");
    tests.expectEqual(eval(session, "qrDecomposition[{{3,0},{4,0}}]"),
        std::string{"{{{-3/5, -4/5}, {-4/5, 3/5}}, {{-5, 0}, {0, 0}}}"},
        "Linear algebra: exact Householder QR preserves rational factors when the norm closes exactly");
    tests.expectEqual(eval(session,
        "N[dot[at[qrDecomposition[{{1,2},{3,4}}],0],at[qrDecomposition[{{1,2},{3,4}}],1]],8]"),
        std::string{"{{1, 2}, {3, 4}}"},
        "Linear algebra: exact Householder QR reconstructs A");
    tests.expectEqual(eval(session, "N[qrDecomposition[{{1,2},{3,4}}],8]"),
        std::string{"{{{-0.31622777, -0.9486833}, {-0.9486833, 0.31622777}}, {{-3.1622777, -4.4271887}, {0, -0.63245553}}}"},
        "Linear algebra: N dispatches QR directly to the certified Householder backend");
    tests.expectEqual(eval(session, "dimensions[at[qrDecomposition[{{1,0},{0,1},{0,0}}],0]]"),
        std::string{"{3, 2}"},
        "Linear algebra: reduced QR exposes a rectangular Q through {Q,R}");
    tests.expectEqual(eval(session, "dimensions[at[qrDecomposition[{{1,0},{0,1},{0,0}}],1]]"),
        std::string{"{2, 2}"},
        "Linear algebra: reduced QR exposes a rectangular R through {Q,R}");
    tests.expectEqual(eval(session,
        "N[dot[at[qrDecomposition[{{1,2},{3,4},{5,6}}],0],at[qrDecomposition[{{1,2},{3,4},{5,6}}],1]],8]"),
        std::string{"{{1, 2}, {3, 4}, {5, 6}}"},
        "Linear algebra: rectangular reduced QR reconstructs A");
    tests.expectEqual(eval(session, "N[luDecomposition[{{Pi,1},{2,3}}],8]"),
        std::string{"{{{1, 0}, {0, 1}}, {{1, 0}, {0.63661977, 1}}, {{3.1415927, 1}, {0, 2.3633802}}}"},
        "Linear algebra: N dispatches LU directly to the certified decomposition backend");
    tests.expectEqual(eval(session, "qrDecomposition[{{a,b},{0,c}}]"),
        std::string{"{{{1, 0}, {0, 1}}, {{a, b}, {0, c}}}"},
        "Linear algebra: already upper-triangular symbolic matrices avoid unnecessary QR expansion");
    tests.expectEqual(eval(session,
        "qrDecomposition[{{19,-1,-2,2},{2,19,0,-1},{-1,-2,19,1},{1,0,-1,19}}]"),
        std::string{"qrDecomposition[{{19, -1, -2, 2}, {2, 19, 0, -1}, {-1, -2, 19, 1}, {1, 0, -1, 19}}]"},
        "Linear algebra: general exact QR stops before the observed 4x4 expression explosion");
    tests.expectEqual(eval(session, "N[det[{{Pi,0},{0,2}}],12]"),
        std::string{"6.28318530718"},
        "Linear algebra: N pushes precision directly into determinant");
    tests.expectEqual(eval(session, "N[dot[{{Pi,0},{0,Pi}},{{1,2},{3,4}}],12]"),
        std::string{"{{3.14159265359, 6.28318530718}, {9.42477796077, 12.5663706144}}"},
        "Linear algebra: N pushes precision directly into dot");
    tests.expectEqual(eval(session, "N[inverse[{{Pi,0},{0,2}}],12]"),
        std::string{"{{0.318309886184, 0}, {0, 0.5}}"},
        "Linear algebra: N uses the certified inverse backend");
    tests.expectEqual(eval(session, "N[rref[{{Pi,0},{0,2}}],12]"),
        std::string{"{{1, 0}, {0, 1}}"},
        "Linear algebra: N uses certified pivots for rref");
    tests.expectEqual(eval(session, "N[norm[{Pi,1}],12]"),
        std::string{"3.29690830948"},
        "Linear algebra: N uses the certified Hermitian norm backend");
    tests.expectEqual(eval(session, "dot[N[{{Pi,0},{0,Pi}},12],{{1,2},{3,4}}]"),
        std::string{"{{3.14159265359, 6.28318530718}, {9.42477796077, 12.5663706144}}"},
        "Linear algebra: approximate Array inputs infer backend precision like FFT");
    tests.expectEqual(eval(session, "N[matrixRank[{{Pi,2Pi},{1,2}}],12]"),
        std::string{"1"},
        "Linear algebra: N preserves exact rank when exact elimination can certify dependence");
    tests.expectEqual(eval(session, "matrixRank[N[{{Pi,2Pi},{1,2}},12]]"),
        std::string{"matrixRank[{{3.14159265359, 6.28318530718}, {1, 2}}]"},
        "Linear algebra: approximate matrixRank does not invent a rank-deficiency threshold");
    tests.expectEqual(eval(session,
        "N[solveLinear[{{Pi,0},{0,2}},{Pi,4}],12]"),
        std::string{"{1, 2}"},
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
        std::string{"{{1, 2}, {3, 4}, {5, 6}}"},
        "Linear algebra: numerical reduced SVD reconstructs a tall matrix");
    tests.expectEqual(eval(session,
        "N[dot[at[svd[{{1,2,3},{4,5,6}}],0],dot[at[svd[{{1,2,3},{4,5,6}}],1],transpose[at[svd[{{1,2,3},{4,5,6}}],2]]]],8]"),
        std::string{"{{1, 2, 3}, {4, 5, 6}}"},
        "Linear algebra: numerical reduced SVD reconstructs a wide matrix");
    tests.expectEqual(eval(session,
        "N[dot[at[svd[{{1,I},{0,1}}],0],dot[at[svd[{{1,I},{0,1}}],1],conjugateTranspose[at[svd[{{1,I},{0,1}}],2]]]],8]"),
        std::string{"{{1, I}, {0.000000000000000I, 1}}"},
        "Linear algebra: complex numerical SVD reconstructs A with V Hermitian-transposed");
    tests.expectEqual(eval(session, "eigenvalues[{{1,2},{3,4}}]"),
        std::string{"{(5+sqrt[33])/2, (5-sqrt[33])/2}"},
        "Linear algebra: exact 2x2 eigenvalues preserve radicals");
    tests.expectEqual(eval(session, "eigenvalues[{{0,-1},{1,0}}]"),
        std::string{"{I, -I}"},
        "Linear algebra: exact 2x2 eigenvalues extend naturally to complex roots");
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
        std::string{"{3.1415927, 2}"},
        "Linear algebra: eigen relation audit uses certified input intervals");
    tests.expectEqual(eval(session, "N[eigenvectors[{{1,1},{0,1}}],8]"),
        std::string{"eigenvectors[{{1, 1}, {0, 1}}]"},
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
