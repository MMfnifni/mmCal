#include "repl_help.hpp"

#include "evaluation/builtin_registry.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <ostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcal::cli {
namespace {

using evaluation::BuiltinId;

struct FunctionHelpEntry final {
    BuiltinId id;
    std::string_view summary;
    std::string_view usage;
    std::string_view inputs;
    std::string_view notes;
    std::string_view examples;
};

#define HELP(id, summary, usage, inputs, examples) \
    {BuiltinId::id, summary, usage, inputs, {}, examples}
#define HELP_NOTE(id, summary, usage, inputs, notes, examples) \
    {BuiltinId::id, summary, usage, inputs, notes, examples}

// Names, aliases, arity, and source-callability remain owned by
// BuiltinRegistry. This catalog supplies the user-facing semantics that cannot
// be inferred from an arity alone. Every source-callable BuiltinId is covered;
// the exhaustive test in kernel_session_tests.cpp guards that invariant.
constexpr FunctionHelpEntry functionHelpEntries[] = {
    HELP_NOTE(Power,
        "Raises a base to an exponent using exact arithmetic and principal complex branches.",
        "pow[base, exponent]",
        "base: An exact, symbolic, or certified numeric value.\n"
        "exponent: The power; integer and rational powers simplify exactly when safe.",
        "Undefined exponent values propagate as Indeterminate. In particular, 0^0, a unit-magnitude base raised to Infinity, and every base raised to 1/0 are indeterminate.",
        "pow[2, 10]  ->  1024\n"
        "N[pow[-8, 1/3], 16]\n"
        "0^I  ->  Indeterminate\n"
        "x^(1/0)  ->  Indeterminate"),
    HELP(Factorial,
        "Computes an exact factorial.",
        "fact[n]",
        "n: A non-negative integer.",
        "fact[10]  ->  3628800"),

    HELP(Derivative,
        "Differentiates an expression symbolically.",
        "D[expression, variable]\n"
        "D[expression, {variable, order}]\n"
        "D[expression, specification1, specification2, ...]",
        "expression: Held while differentiation variables are bound.\n"
        "variable: A symbol.\n"
        "order: A non-negative integer; specifications are applied left to right.",
        "D[exp[x^2], x]  ->  2x exp[x^2]\n"
        "D[x^2 y^3, x, y]  ->  6x y^2\n"
        "D[sin[x], {x, 4}]  ->  sin[x]"),
    HELP_NOTE(SymbolicIntegral,
        "Integrates an expression symbolically, optionally over exact bounds and assumptions.",
        "integrate[expression, variable]\n"
        "integrate[expression, {variable, a, b}]\n"
        "integrate[expression, variable, assumptions]\n"
        "integrate[expression, {variable, a, b}, assumptions]",
        "expression: Held while the integration variable is locally bound.\n"
        "variable: A symbol, or {symbol, lower, upper} for a definite integral.\n"
        "assumptions: A predicate such as x >= 0 or element[x, Real].",
        "Indefinite results omit the arbitrary additive constant. Singularities are checked over a definite interval; Cauchy principal values are not assumed.",
        "integrate[x^2, x]  ->  x^3/3\n"
        "integrate[sin[x], {x, 0, Pi}]  ->  2\n"
        "integrate[abs[x], x, x >= 0]  ->  x^2/2"),
    HELP(Limit,
        "Computes an exact or symbolic limit, including one-sided limits.",
        "limit[expression, variable, point]\n"
        "limit[expression, variable, point, direction]\n"
        "limit[expression, {variable, point, direction}]",
        "expression: Held while variable is locally bound.\n"
        "variable: The symbol tending to point.\n"
        "direction: Optional -1 for the left or 1 for the right; omission requests a two-sided limit. The compact brace form is equivalent to the four-argument form. Proven periodic oscillation with no single limiting value returns Indeterminate.",
        "limit[sin[x]/x, x, 0]  ->  1\n"
        "limit[sin[1/x], x, 0, 1]  ->  Indeterminate\n"
        "limit[x sin[1/x], x, 0]  ->  0\n"
        "limit[1/x, x, 0, 1]  ->  Infinity\n"
        "limit[li[x], {x, 1, 1}]  ->  -Infinity"),

    HELP(Floor,
        "Returns the greatest integer not greater than a real value.",
        "floor[x]",
        "x: An exact or certifiably ordered real value.",
        "floor[-3/2]  ->  -2"),
    HELP(Ceil,
        "Returns the least integer not less than a real value.",
        "ceil[x]",
        "x: An exact or certifiably ordered real value.",
        "ceil[-3/2]  ->  -1"),
    HELP(Trunc,
        "Removes the fractional part by rounding toward zero.",
        "trunc[x]",
        "x: An exact or certified real value.",
        "trunc[-3/2]  ->  -1"),
    HELP(Round,
        "Rounds a real value to nearest-even, optionally at a decimal place.",
        "round[x]\n"
        "round[x, digits]",
        "x: An exact or certified real value.\n"
        "digits: An integer; the rounding quantum is 10^(-digits), so negative values round left of the decimal point.",
        "round[5/2]  ->  2\n"
        "round[135, -1]  ->  140"),
    HELP(Frac,
        "Returns the non-negative fractional part x-floor[x].",
        "frac[x]",
        "x: An exact or certified real value.",
        "frac[-3/2]  ->  1/2"),
    HELP_NOTE(BitAnd,
        "Computes the bitwise AND of two or more arbitrary-size integers.",
        "bitand[a, b, ...]",
        "a, b, ...: Integers.",
        "Negative integers use infinite two's-complement sign extension.",
        "bitand[-1, 5]  ->  5"),
    HELP_NOTE(BitOr,
        "Computes the bitwise OR of two or more arbitrary-size integers.",
        "bitor[a, b, ...]",
        "a, b, ...: Integers.",
        "Negative integers use infinite two's-complement sign extension.",
        "bitor[-8, 3]  ->  -5"),
    HELP_NOTE(BitXor,
        "Computes the bitwise exclusive OR of two or more arbitrary-size integers.",
        "bitxor[a, b, ...]",
        "a, b, ...: Integers.",
        "Negative integers use infinite two's-complement sign extension.",
        "bitxor[-1, 5]  ->  -6"),
    HELP_NOTE(BitNot,
        "Computes the bitwise complement of an arbitrary-size integer.",
        "bitnot[n]",
        "n: An integer.",
        "Uses infinite two's-complement semantics, so bitnot[n] equals -n-1.",
        "bitnot[5]  ->  -6"),
    HELP(BitShiftLeft,
        "Shifts an integer left by a non-negative bit count.",
        "bitshiftl[n, count]",
        "n: An integer.\n"
        "count: A non-negative integer.",
        "bitshiftl[5, 3]  ->  40"),
    HELP_NOTE(BitShiftRight,
        "Performs an arithmetic right shift.",
        "bitshiftr[n, count]",
        "n: An integer.\n"
        "count: A non-negative integer.",
        "Negative values are sign-extended.",
        "bitshiftr[-3, 1]  ->  -2"),
    HELP(BitLength,
        "Returns the number of bits needed for an integer magnitude.",
        "bitlength[n]",
        "n: An integer.",
        "bitlength[255]  ->  8"),
    HELP(BitCount,
        "Counts one bits in a non-negative integer.",
        "bitcount[n]",
        "n: A non-negative integer; negative infinite two's-complement values are rejected.",
        "bitcount[15]  ->  4"),
    HELP_NOTE(BitGet,
        "Returns one bit of an integer.",
        "bitget[n, index]",
        "n: An integer.\n"
        "index: A zero-based non-negative bit position.",
        "Negative integers are sign-extended.",
        "bitget[-2, 100]  ->  1"),
    HELP(Gcd,
        "Returns the greatest common divisor of two or more integers.",
        "gcd[a, b, ...]",
        "a, b, ...: Integers.",
        "gcd[84, 126, 210]  ->  42"),
    HELP(Lcm,
        "Returns the least common multiple of two or more integers.",
        "lcm[a, b, ...]",
        "a, b, ...: Integers.",
        "lcm[6, 8, 9]  ->  72"),
    HELP_NOTE(Mod,
        "Returns the modulus associated with floor division.",
        "mod[a, m]",
        "a: An integer dividend.\n"
        "m: A nonzero integer modulus.",
        "Unlike rem, mod uses a floor quotient.",
        "mod[-5, 3]  ->  1"),
    HELP_NOTE(Rem,
        "Returns the remainder associated with truncation toward zero.",
        "rem[a, b]",
        "a: An integer dividend.\n"
        "b: A nonzero integer divisor.",
        "Unlike mod, rem uses a truncate-toward-zero quotient.",
        "rem[-5, 3]  ->  -2"),
    HELP(Quotient,
        "Returns the integer quotient rounded toward zero.",
        "quotient[a, b]",
        "a: An integer dividend.\n"
        "b: A nonzero integer divisor.",
        "quotient[-5, 3]  ->  -1"),
    HELP_NOTE(IsPrime,
        "Tests primality exactly within the supported unsigned 64-bit proof range.",
        "isprime[n]",
        "n: A non-negative integer no greater than 2^64-1.",
        "Larger inputs remain unevaluated rather than returning a probable-prime answer.",
        "isprime[97]  ->  True"),
    HELP(NextPrime,
        "Returns the least prime strictly greater than an integer.",
        "nextprime[n]",
        "n: An integer whose search result fits the supported exact-primality range.",
        "nextprime[14]  ->  17"),
    HELP(PreviousPrime,
        "Returns the greatest prime strictly less than an integer.",
        "prevprime[n]",
        "n: An integer greater than 2.",
        "prevprime[14]  ->  13"),
    HELP_NOTE(FactorInteger,
        "Factors an integer into a flat list of exact prime factors.",
        "factorint[n]",
        "n: A nonzero integer whose magnitude is within the supported unsigned 64-bit factorization range.",
        "A negative input prefixes -1; zero is a DomainError.",
        "factorint[360]  ->  {2, 2, 2, 3, 3, 5}\n"
        "factorint[-12]  ->  {-1, 2, 2, 3}"),
    HELP(Totient,
        "Computes Euler's totient: the count of residues coprime to n.",
        "totient[n]",
        "n: A positive integer within the supported exact factorization range.",
        "totient[9]  ->  6"),
    HELP(Permutation,
        "Counts ordered selections of r items from n.",
        "perm[n, r]",
        "n: A non-negative integer.\n"
        "r: A non-negative integer; r>n returns 0.",
        "perm[10, 3]  ->  720\n"
        "perm[5, 6]  ->  0"),
    HELP(Combination,
        "Computes the exact binomial count 'n choose r'.",
        "comb[n, r]",
        "n: A non-negative integer.\n"
        "r: A non-negative integer; r>n returns 0.",
        "comb[10, 3]  ->  120\n"
        "comb[5, 6]  ->  0"),
    HELP(Fibonacci,
        "Computes an exact Fibonacci number by fast doubling.",
        "fib[n]",
        "n: A non-negative integer.",
        "fib[100]  ->  354224848179261915075"),

    HELP_NOTE(DiscreteFourierTransform,
        "Computes the exact-first discrete Fourier transform directly.",
        "dft[data]",
        "data: A rank-1 array of exact, symbolic, or certified numeric values; an empty array returns {}.",
        "The phase is always in radians and is independent of angleMode[].",
        "dft[{1, 2, 3, 4}]  ->  {10, -2+2I, -2, -2-2I}"),
    HELP_NOTE(FastFourierTransform,
        "Computes an exact-first discrete Fourier transform with optimized algorithms where applicable.",
        "fft[data]\n"
        "N[fft[data], digits]",
        "data: A rank-1 array; an empty array returns {}. Exact input stays exact, while approximate input uses certified interval arithmetic.\n"
        "digits: When wrapped in N, the requested significant decimal digits.",
        "Power-of-two exact lengths use radix-2; supported non-power-of-two exact data use cyclotomic arithmetic. Certified non-power-of-two transforms select direct DFT or Bluestein by the configured threshold.",
        "fft[{1, 0, 0, 0}]  ->  {1, 1, 1, 1}\n"
        "N[fft[{1, 2, 3}], 20]"),
    HELP_NOTE(InverseFourierTransform,
        "Computes the inverse discrete Fourier transform.",
        "ifft[data]",
        "data: A rank-1 array of exact, symbolic, or certified numeric values; an empty array returns {}.",
        "The inverse includes the 1/n normalization and uses radian phase semantics.",
        "ifft[fft[{1, 2, 3, 4}]]  ->  {1, 2, 3, 4}"),
    HELP(Convolution,
        "Computes the full discrete linear convolution of two sequences.",
        "convolve[a, b]",
        "a: A rank-1 array.\n"
        "b: A rank-1 array. If either input is empty, the result is {}.",
        "convolve[{1, 2}, {3, 4}]  ->  {3, 10, 8}"),

    HELP(Transpose,
        "Transposes a dense rectangular matrix.",
        "transpose[matrix]",
        "matrix: A rank-2 rectangular array; ragged braces are not accepted.",
        "transpose[{{1, 2, 3}, {4, 5, 6}}]  ->  {{1, 4}, {2, 5}, {3, 6}}"),
    HELP_NOTE(ConjugateTranspose,
        "Computes the Hermitian transpose (complex conjugate plus transpose).",
        "conjugateTranspose[array]",
        "array: A rank-1 vector or rank-2 rectangular matrix.",
        "A vector is conjugated componentwise; a matrix is also transposed.",
        "conjugateTranspose[{{1, I}, {2, 3}}]  ->  {{1, 2}, {-I, 3}}"),
    HELP(MatrixAdd,
        "Adds two or more matrices element by element.",
        "madd[matrix1, matrix2, ...]",
        "matrix1, matrix2, ...: Rectangular matrices with identical dimensions.",
        "madd[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]  ->  {{6, 8}, {10, 12}}"),
    HELP_NOTE(MatrixMultiply,
        "Contracts vectors and matrices: vector dot product or matrix product.",
        "dot[a, b]",
        "a, b: Rank-1 or rank-2 arrays with compatible inner dimensions.",
        "Array-by-Array * is not matrix multiplication; use dot explicitly.",
        "dot[{1, 2, 3}, {4, 5, 6}]  ->  32\n"
        "dot[{{1, 2}, {3, 4}}, {5, 6}]  ->  {17, 39}\n"
        "dot[{{1, 2}, {3, 4}}, {{5, 6}, {7, 8}}]  ->  {{19, 22}, {43, 50}}"),
    HELP(Determinant,
        "Computes an exact-first matrix determinant.",
        "det[matrix]",
        "matrix: A square rectangular matrix. Exact rational and complex entries remain exact.",
        "det[{{1, 2}, {3, 4}}]  ->  -2"),
    HELP(Inverse,
        "Computes an exact-first inverse of a nonsingular square matrix.",
        "inverse[matrix]",
        "matrix: A square rectangular matrix with a provably nonzero determinant.",
        "inverse[{{1, 2}, {3, 4}}]  ->  {{-2, 1}, {3/2, -1/2}}"),
    HELP(Rref,
        "Returns reduced row-echelon form using exact or certified pivots.",
        "rref[matrix]",
        "matrix: A rectangular matrix; a symbolic pivot is used only when nonzero status is provable.",
        "rref[{{1, 2}, {3, 4}}]  ->  {{1, 0}, {0, 1}}"),
    HELP(Rank,
        "Computes matrix rank without a floating epsilon guess.",
        "matrixRank[matrix]",
        "matrix: A rectangular matrix. Approximate pivots must be certified; undecidable symbolic cases remain unevaluated.",
        "matrixRank[{{1, 2}, {2, 4}}]  ->  1"),
    HELP_NOTE(SolveLinear,
        "Solves a linear system A x = b when the solution is unique.",
        "solveLinear[A, b]",
        "A: An m by n rectangular matrix.\n"
        "b: A vector of length m.",
        "Consistent overdetermined full-column-rank systems are accepted. Inconsistent or underdetermined systems are DomainErrors rather than parametric answers.",
        "solveLinear[{{2, 1}, {1, -1}}, {5, 1}]  ->  {2, 1}\n"
        "solveLinear[{{1, 0}, {0, 1}, {1, 1}}, {2, 3, 5}]  ->  {2, 3}"),
    HELP_NOTE(NullSpace,
        "Returns a canonical row-basis for the right null space of a matrix.",
        "nullSpace[matrix]",
        "matrix: An m by n rectangular matrix with exact or certifiable pivot structure.",
        "Basis rows follow free columns in ascending order; result dimensions are {nullity, n}.",
        "nullSpace[{{1, 2}, {2, 4}}]  ->  {{-2, 1}}\n"
        "nullSpace[{{1, 0}, {0, 1}}]  ->  reshape[{}, {0, 2}]"),
    HELP_NOTE(LuDecomposition,
        "Computes a pivoted LU decomposition of a square matrix.",
        "luDecomposition[matrix]\n"
        "N[luDecomposition[matrix], digits]",
        "matrix: A square exact or certified numeric matrix.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Returns {P, L, U} with P matrix = L U. Extract factors with at[result, 0], at[result, 1], and at[result, 2].",
        "lu := luDecomposition[{{0, 1}, {2, 3}}]\n"
        "dot[at[lu, 0], {{0, 1}, {2, 3}}] == dot[at[lu, 1], at[lu, 2]]  ->  True"),
    HELP_NOTE(QrDecomposition,
        "Computes a reduced QR decomposition.",
        "qrDecomposition[matrix]\n"
        "N[qrDecomposition[matrix], digits]",
        "matrix: An m by n rectangular real or complex matrix.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "With k=min(m,n), returns {Q,R}, where Q is m by k, R is k by n, and matrix = Q R. Exact real matrices use fraction-free orthogonalization and delay square roots until Q/R materialization; N dispatches directly to the certified Householder backend.",
        "qr := qrDecomposition[{{1, 0}, {0, 1}, {1, 1}}]\n"
        "dot[at[qr, 0], at[qr, 1]]  ->  {{1, 0}, {0, 1}, {1, 1}}\n"
        "N[qrDecomposition[{{1, Pi}, {2, 3}}], 16]"),
    HELP_NOTE(SingularValueDecomposition,
        "Computes a reduced singular value decomposition.",
        "svd[matrix]\n"
        "N[svd[matrix], digits]",
        "matrix: An m by n real or complex matrix.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Returns {U,S,V} with k=min(m,n). Real input satisfies A=U S transpose[V]; complex input uses conjugateTranspose[V]. General numeric results are returned only after reconstruction and orthogonality checks pass.",
        "svd[{{3, 0}, {0, 4}}]\n"
        "N[svd[{{1, 2}, {3, 4}}], 16]\n"
        "s := svd[{{3, 0}, {0, 4}}]\n"
        "dot[dot[at[s, 0], at[s, 1]], transpose[at[s, 2]]]"),
    HELP_NOTE(ConditionNumber,
        "Computes the spectral condition number sigma_max/sigma_min.",
        "conditionNumber[matrix]\n"
        "N[conditionNumber[matrix], digits]",
        "matrix: A non-empty rectangular real or complex matrix.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Exact rank deficiency returns Infinity. Exact closed cases remain exact; general numerical evaluation uses the certified SVD backend and does not guess rank from a floating epsilon.",
        "conditionNumber[{{3, 0}, {0, 4}}]  ->  4/3\n"
        "conditionNumber[{{1, 2}, {2, 4}}]  ->  Infinity\n"
        "N[conditionNumber[{{1, 2}, {3, 4}}], 12]"),
    HELP_NOTE(PseudoInverse,
        "Computes the Moore-Penrose pseudoinverse.",
        "pseudoInverse[matrix]\n"
        "N[pseudoInverse[matrix], digits]",
        "matrix: A rectangular real or complex matrix.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Exact numeric matrices use exact rank factorization, including rank-deficient cases. General full-rank numerical matrices use certified SVD. A finite-precision rank deficiency is never inferred from an arbitrary threshold.",
        "pseudoInverse[{{1, 2}, {2, 4}}]  ->  {{1/25, 2/25}, {2/25, 4/25}}\n"
        "pseudoInverse[{{I, 0}, {0, 2I}}]  ->  {{-I, 0}, {0, -I/2}}\n"
        "N[pseudoInverse[{{1, 2}, {3, 4}}], 12]"),
    HELP_NOTE(LeastSquares,
        "Returns the minimum-norm least-squares solution of A x approximately equals b.",
        "leastSquares[A, b]\n"
        "N[leastSquares[A, b], digits]",
        "A: An m by n rectangular real or complex matrix.\n"
        "b: A vector of length m.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "The result is A^+ b. Exact numeric inputs remain exact, including consistent or rank-deficient systems. General full-rank numerical inputs use the certified SVD pseudoinverse; uncertain numerical rank remains unevaluated.",
        "leastSquares[{{1, 0}, {0, 1}, {1, 1}}, {1, 2, 4}]  ->  {4/3, 7/3}\n"
        "leastSquares[{{1, 2}, {2, 4}}, {1, 2}]  ->  {1/5, 2/5}\n"
        "N[leastSquares[{{1, 0}, {0, 1}, {1, 1}}, {1, 2, 4}], 12]"),
    HELP_NOTE(Eigenvalues,
        "Computes the eigenvalues of a square matrix.",
        "eigenvalues[matrix]\n"
        "N[eigenvalues[matrix], digits]",
        "matrix: A square exact or certified numeric matrix.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Exact closed cases remain exact; general N uses a certified Schur/eigenpair backend.",
        "eigenvalues[{{2, 0}, {0, 3}}]  ->  {2, 3}\n"
        "N[eigenvalues[{{1, 2}, {3, 4}}], 16]"),
    HELP_NOTE(Eigenvectors,
        "Computes matrix eigenvectors and returns them as columns.",
        "eigenvectors[matrix]\n"
        "N[eigenvectors[matrix], digits]",
        "matrix: A square exact or certified numeric matrix.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Repeated or defective cases remain unevaluated when a stable independent basis cannot be certified.",
        "eigenvectors[{{2, 0}, {0, 3}}]  ->  {{1, 0}, {0, 1}}\n"
        "N[eigenvectors[{{1, 2}, {3, 4}}], 16]"),
    HELP_NOTE(Eigensystem,
        "Computes paired eigenvalues and column eigenvectors.",
        "eigensystem[matrix]\n"
        "N[eigensystem[matrix], digits]",
        "matrix: A square exact or certified numeric matrix.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Returns {values, vectors}; column j of vectors corresponds to values[[j]] (extract with at).",
        "eigensystem[{{2, 0}, {0, 3}}]  ->  {{2, 3}, {{1, 0}, {0, 1}}}\n"
        "N[eigensystem[{{1, 2}, {3, 4}}], 16]"),

    HELP_NOTE(NumericDerivative,
        "Evaluates a symbolic derivative at a point with certified decimal arithmetic.",
        "diff[expression, variable, point]\n"
        "diff[expression, variable, point, digits]",
        "expression: Held while variable is bound and first differentiated by D.\n"
        "variable: A symbol.\n"
        "point: An exact value at which to evaluate.\n"
        "digits: Optional significant decimal digits; default 16.",
        "This is not a finite-difference formula: mmCal constructs D[expression,variable] first.",
        "diff[x^2, x, 3]  ->  6.0\n"
        "diff[sin[x], x, 0, 20]  ->  1.0000000000000000000"),
    HELP_NOTE(NumericIntegral,
        "Computes a certified numerical definite integral.",
        "nintegrate[expression, {variable, a, b}]\n"
        "nintegrate[expression, {variable, a, b}, digits]",
        "expression: Held while variable is locally bound.\n"
        "{variable,a,b}: A three-item iterator with finite exact bounds.\n"
        "digits: Optional significant decimal digits; default 16.",
        "Certified quadrature and derivative bounds must prove the final rounding. The whole interval is preflighted for singularities.",
        "nintegrate[x^2, {x, 0, 1}, 12]  ->  0.333333333333\n"
        "nintegrate[sin[x], {x, 0, Pi}, 12]  ->  2.00000000000"),

    HELP(Cbrt,
        "Computes the real cube root of a real value.",
        "cbrt[x]",
        "x: A real exact or certified numeric value; negative values are allowed.",
        "cbrt[-8]  ->  -2"),
    HELP(Hypot,
        "Computes sqrt[x^2+y^2] without discarding exact structure.",
        "hypot[x, y]",
        "x, y: Real numeric values.",
        "hypot[3, 4]  ->  5"),
    HELP_NOTE(Fma,
        "Computes a*b+c as one fused semantic operation.",
        "fma[a, b, c]",
        "a, b, c: Exact or certified numeric scalars.",
        "Exact operands stay exact; certified operands avoid an intermediate DecimalApproximation rounding.",
        "fma[2, 3, 4]  ->  10"),
    HELP(Clamp,
        "Clamps a real value to an inclusive interval.",
        "clamp[x, lower, upper]",
        "x: A real value.\n"
        "lower, upper: Real bounds with lower <= upper; approximate ordering must be certified.",
        "clamp[7, 0, 5]  ->  5"),
    HELP_NOTE(Proj,
        "Applies the complex projective projection helper.",
        "proj[z]",
        "z: A finite exact or certified complex value.",
        "For the currently supported finite values this is the identity; full Riemann-sphere infinity semantics are deferred.",
        "proj[3+4I]  ->  3+4I"),
    HELP(Cis,
        "Constructs the unit complex number cos[x]+I sin[x].",
        "cis[x]",
        "x: An angle; a Deg, Rad, or Grad suffix overrides angleMode[].",
        "cis[Pi/3]  ->  1/2+I sqrt[3]/2"),
    HELP(Polar,
        "Constructs a complex value from magnitude and angle as r*cis[theta].",
        "polar[r, theta]",
        "r: A real magnitude.\n"
        "theta: An angle; explicit Deg, Rad, or Grad overrides angleMode[].",
        "polar[2, Pi/3]  ->  1+I sqrt[3]"),
    HELP(NextPow2,
        "Returns the smallest exponent n for which 2^n is at least x.",
        "nextpow2[x]",
        "x: A positive real numeric value.",
        "nextpow2[9]  ->  4"),
    HELP(DegreeToRadian,
        "Converts a degree measure to radians.",
        "DtoR[x]",
        "x: A real degree value; do not append an angle-unit suffix.",
        "DtoR[180]  ->  Pi"),
    HELP(DegreeToGradian,
        "Converts a degree measure to gradians.",
        "DtoG[x]",
        "x: A real degree value; do not append an angle-unit suffix.",
        "DtoG[90]  ->  100"),
    HELP(RadianToDegree,
        "Converts a radian measure to degrees.",
        "RtoD[x]",
        "x: A real radian value; do not append an angle-unit suffix.",
        "RtoD[Pi]  ->  180"),
    HELP(RadianToGradian,
        "Converts a radian measure to gradians.",
        "RtoG[x]",
        "x: A real radian value; do not append an angle-unit suffix.",
        "RtoG[Pi]  ->  200"),
    HELP(GradianToDegree,
        "Converts a gradian measure to degrees.",
        "GtoD[x]",
        "x: A real gradian value; do not append an angle-unit suffix.",
        "GtoD[200]  ->  180"),
    HELP(GradianToRadian,
        "Converts a gradian measure to radians.",
        "GtoR[x]",
        "x: A real gradian value; do not append an angle-unit suffix.",
        "GtoR[200]  ->  Pi"),

    HELP_NOTE(Sum,
        "Adds scalar arguments or all elements of one array.",
        "sum[]\n"
        "sum[x1, x2, ...]\n"
        "sum[array]",
        "x1, x2, ...: Scalar expressions, or supply one array to aggregate its elements.",
        "Symbolic iterator syntax sum[f,{k,a,b}] is not implemented; generate a finite table first.",
        "sum[]  ->  0\n"
        "sum[1, 2, 3]  ->  6\n"
        "sum[{1, 2, 3}]  ->  6"),
    HELP_NOTE(Product,
        "Multiplies scalar arguments or all elements of one array.",
        "prod[]\n"
        "prod[x1, x2, ...]\n"
        "prod[array]",
        "x1, x2, ...: Scalar expressions, or supply one array to aggregate its elements.",
        "Symbolic iterator syntax prod[f,{k,a,b}] is not implemented; generate a finite table first.",
        "prod[]  ->  1\n"
        "prod[2, 3, 4]  ->  24\n"
        "prod[{2, 3, 4}]  ->  24"),
    HELP_NOTE(Map,
        "Applies a function explicitly to every scalar leaf of an array or brace.",
        "map[function, value]",
        "function: A builtin or user-function symbol, written without brackets.\n"
        "value: A dense array or general ragged brace value.",
        "Dense shape and ragged brace structure are preserved. Ordinary calls such as exp[array] are not implicitly element-wise.",
        "map[sin, {0, Pi/2, Pi}]  ->  {0, 1, 0}"),
    HELP(Range,
        "Builds an inclusive exact arithmetic progression.",
        "range[n]\n"
        "range[a, b]\n"
        "range[a, b, step]",
        "n: An exact endpoint, producing 1 through n.\n"
        "a, b: Exact Integer or Rational endpoints.\n"
        "step: A nonzero exact Integer or Rational step; default 1.",
        "range[5]  ->  {1, 2, 3, 4, 5}\n"
        "range[0, 1, 1/3]  ->  {0, 1/3, 2/3, 1}\n"
        "range[5, 1, -2]  ->  {5, 3, 1}"),
    HELP_NOTE(Table,
        "Evaluates a held expression over an exact local iterator.",
        "table[expression, {variable, count}]\n"
        "table[expression, {variable, a, b}]\n"
        "table[expression, {variable, a, b, step}]",
        "expression: Held and reevaluated for each iterator value.\n"
        "variable: A locally scoped symbol.\n"
        "count, a, b, step: Exact Integer or Rational iterator controls; step must be nonzero.",
        "The iterator shadows an outer definition only inside the table and is restored even if evaluation fails.",
        "table[i^2, {i, 5}]  ->  {1, 4, 9, 16, 25}\n"
        "table[i/2, {i, 0, 2, 1/2}]  ->  {0, 1/4, 1/2, 3/4, 1}"),
    HELP(Min,
        "Returns the least provably ordered real argument or array element.",
        "min[x1, x2, ...]\n"
        "min[array]",
        "Provide one or more real values, or one array. Symbolic order must be provable.",
        "min[3, 1, 2]  ->  1\n"
        "min[x, 3]  ->  min[x, 3]"),
    HELP(Max,
        "Returns the greatest provably ordered real argument or array element.",
        "max[x1, x2, ...]\n"
        "max[array]",
        "Provide one or more real values, or one array. Symbolic order must be provable.",
        "max[3, 1, 2]  ->  3"),
    HELP(Mean,
        "Computes the exact arithmetic mean of arguments or one array.",
        "mean[x1, x2, ...]\n"
        "mean[array]",
        "Provide one or more numeric values, or one non-empty array.",
        "mean[1, 2, 4]  ->  7/3"),

    HELP(Median,
        "Returns the middle value, averaging the middle pair for even data.",
        "median[x1, x2, ...]\n"
        "median[array]",
        "Provide one or more exact real observations, or one non-empty rank-1 array.",
        "median[{1, 2, 3, 4}]  ->  5/2"),
    HELP(Mode,
        "Returns the most frequent value or all tied modes.",
        "mode[x1, x2, ...]\n"
        "mode[array]",
        "Provide one or more comparable observations, or one non-empty rank-1 array.",
        "mode[1, 1, 2, 2]  ->  {1, 2}\n"
        "mode[1, 2, 3]  ->  {}"),
    HELP_NOTE(Quantile,
        "Computes a Type-7 quantile with exact interpolation.",
        "quantile[p, x1, x2, ...]\n"
        "quantile[p, array]",
        "p: An exact real probability in [0,1].\n"
        "data: One or more exact real observations, directly or as one rank-1 array.",
        "The probability is the first argument.",
        "quantile[1/4, 1, 2, 3, 4, 5, 6, 7]  ->  5/2"),
    HELP_NOTE(Percentile,
        "Computes a Type-7 percentile using a 0-to-100 scale.",
        "percentile[p, x1, x2, ...]\n"
        "percentile[p, array]",
        "p: An exact real percentage in [0,100].\n"
        "data: One or more exact real observations, directly or as one rank-1 array.",
        "percentile[p,...] is quantile[p/100,...].",
        "percentile[50, {1, 3, 5}]  ->  3"),
    HELP(VariancePopulation,
        "Computes population variance, dividing by n.",
        "var[x1, x2, ...]\n"
        "var[array]",
        "Provide one or more exact real observations, or one rank-1 array.",
        "var[1, 2, 3]  ->  2/3"),
    HELP(VarianceSample,
        "Computes unbiased sample variance, dividing by n-1.",
        "vars[x1, x2, ...]\n"
        "vars[array]",
        "Provide at least two exact real observations, or one rank-1 array of that length.",
        "vars[{1, 2, 3}]  ->  1"),
    HELP(StddevPopulation,
        "Computes the square root of population variance.",
        "stddev[x1, x2, ...]\n"
        "stddev[array]",
        "Provide one or more exact real observations, or one rank-1 array.",
        "stddev[1, 2, 3]  ->  sqrt[6]/3"),
    HELP(StddevSample,
        "Computes the square root of unbiased sample variance.",
        "stddevs[x1, x2, ...]\n"
        "stddevs[array]",
        "Provide at least two exact real observations, or one rank-1 array of that length.",
        "stddevs[{1, 2, 3}]  ->  1"),
    HELP(GeometricMean,
        "Computes the real geometric mean while preserving exact roots.",
        "geomean[x1, x2, ...]\n"
        "geomean[array]",
        "Provide one or more non-negative exact real observations, or one rank-1 array.",
        "geomean[1, 4, 1/32]  ->  1/2"),
    HELP(HarmonicMean,
        "Computes the harmonic mean n/sum[1/x].",
        "harmmean[x1, x2, ...]\n"
        "harmmean[array]",
        "Provide one or more nonzero exact real observations, or one rank-1 array.",
        "harmmean[1, 2, 6]  ->  9/5"),
    HELP(Rms,
        "Computes the root mean square.",
        "rms[x1, x2, ...]\n"
        "rms[array]",
        "Provide one or more exact real observations, or one rank-1 array.",
        "rms[1, -1, 1, -1]  ->  1"),
    HELP(MedianAbsoluteDeviation,
        "Computes the median absolute deviation from the median.",
        "mad[x1, x2, ...]\n"
        "mad[array]",
        "Provide one or more exact real observations, or one rank-1 array.",
        "mad[1, 1, 2, 2, 4]  ->  1"),
    HELP(MeanAbsoluteDeviation,
        "Computes the mean absolute deviation from the arithmetic mean.",
        "madR[x1, x2, ...]\n"
        "madR[array]",
        "Provide one or more exact real observations, or one rank-1 array.",
        "madR[1, 2, 3]  ->  2/3"),
    HELP(Skewness,
        "Computes population-moment skewness.",
        "skew[x1, x2, ...]\n"
        "skew[array]",
        "Provide exact real observations with nonzero population variance.",
        "skew[1, 2, 3, 4, 5]  ->  0"),
    HELP(KurtosisPopulation,
        "Computes population excess kurtosis.",
        "kurtp[x1, x2, ...]\n"
        "kurtp[array]",
        "Provide exact real observations with nonzero population variance.",
        "kurtp[-2, -1, 0, 1, 2]  ->  -13/10"),
    HELP(KurtosisSample,
        "Computes unbiased Fisher sample excess kurtosis.",
        "kurts[x1, x2, ...]\n"
        "kurts[array]",
        "Provide a sufficiently large exact real sample with nonzero variance.",
        "kurts[-2, -1, 0, 1, 2]  ->  -6/5"),
    HELP(CoefficientVariation,
        "Computes population standard deviation divided by the mean.",
        "cv[x1, x2, ...]\n"
        "cv[array]",
        "Provide at least two exact real observations with a nonzero mean.",
        "cv[10, 10, 10]  ->  0"),
    HELP(StandardError,
        "Computes sample standard deviation divided by sqrt[n].",
        "stderr[x1, x2, ...]\n"
        "stderr[array]",
        "Provide at least two exact real observations, or one rank-1 array.",
        "stderr[{1, 2, 3}]  ->  sqrt[3]/3"),
    HELP(ZScore,
        "Standardizes a value as (x-mean)/sigma.",
        "zscore[x, mean, sigma]",
        "x: A real value.\n"
        "mean: The population mean.\n"
        "sigma: A positive real standard deviation.",
        "zscore[5, 3, 1]  ->  2"),
    HELP(Iqr,
        "Computes the Type-7 interquartile range Q3-Q1.",
        "iqr[x1, x2, ...]\n"
        "iqr[array]",
        "Provide at least two exact real observations, or one rank-1 array of that length.",
        "iqr[1, 2, 3, 4]  ->  3/2"),
    HELP_NOTE(TrimMean,
        "Computes a mean after removing an equal fraction from each tail.",
        "trimmean[p, x1, x2, ...]\n"
        "trimmean[p, array]",
        "p: An exact fraction; floor[p*n] observations are removed from each tail.\n"
        "data: Exact real observations.",
        "The trim fraction is the first argument.",
        "trimmean[1/5, 1, 2, 100, 3, 4]  ->  3"),
    HELP_NOTE(WinsorMean,
        "Computes the mean after symmetric winsorization.",
        "winsor[p, x1, x2, ...]\n"
        "winsor[p, array]",
        "p: An exact tail fraction.\n"
        "data: Exact real observations.",
        "Tail observations are clipped to the retained boundary values instead of removed.",
        "winsor[1/5, 1, 2, 100, 3, 4]  ->  3"),
    HELP_NOTE(Winsorized,
        "Returns symmetrically winsorized data in its original order.",
        "winsorR[p, x1, x2, ...]\n"
        "winsorR[p, array]",
        "p: An exact tail fraction.\n"
        "data: Exact real observations.",
        "Unlike winsor, this returns the clipped observations rather than their mean.",
        "winsorR[1/5, 1, 2, 100, 3, 4]  ->  {2, 2, 4, 3, 4}"),
    HELP_NOTE(Covariance,
        "Computes exact population covariance of two datasets.",
        "cov[xData, yData]\n"
        "cov[x1, x2, ..., y1, y2, ...]",
        "xData, yData: Equal-length rank-1 arrays.\n"
        "Legacy scalar form: An even count split into equal first and second halves.",
        "The result divides by n, not n-1.",
        "cov[{1, 2, 3}, {2, 4, 6}]  ->  4/3"),
    HELP_NOTE(Correlation,
        "Computes Pearson product-moment correlation.",
        "corr[xData, yData]\n"
        "corr[x1, x2, ..., y1, y2, ...]",
        "xData, yData: Equal-length rank-1 arrays with nonzero variance.\n"
        "Legacy scalar form: An even count split into equal halves.",
        "Exact input preserves exact radicals and rational simplifications.",
        "corr[{1, 2, 3}, {2, 4, 6}]  ->  1"),
    HELP_NOTE(SpearmanCorrelation,
        "Computes Spearman rank correlation with exact average ranks for ties.",
        "corrspearman[xData, yData]\n"
        "corrspearman[x1, x2, ..., y1, y2, ...]",
        "xData, yData: Equal-length rank-1 arrays with nonzero rank variance.\n"
        "Legacy scalar form: An even count split into equal halves.",
        "Tied observations receive average ranks.",
        "corrspearman[{1, 1, 2}, {10, 10, 20}]  ->  1"),
    HELP_NOTE(PercentRank,
        "Returns the interpolated percentile rank of a value in a dataset.",
        "percentrank[x, x1, x2, ...]\n"
        "percentrank[x, array]",
        "x: The exact real value to rank.\n"
        "data: Exact real observations, directly or as one rank-1 array.",
        "Observed ties use their average rank; values between observations are linearly interpolated.",
        "percentrank[3, 1, 2, 3, 4, 5]  ->  1/2\n"
        "percentrank[5/2, 1, 2, 3, 4]  ->  1/2"),

    HELP_NOTE(Dimensions,
        "Returns the rectangular dimensions known for an array or brace.",
        "dimensions[value]",
        "value: A dense array or general brace value.",
        "For a ragged brace, only the common rectangular prefix is reported.",
        "dimensions[{{1, 2, 3}, {4, 5, 6}}]  ->  {2, 3}\n"
        "dimensions[{{1, 2}, {3}}]  ->  {2}"),
    HELP_NOTE(ArrayRank,
        "Returns the number of common rectangular dimensions.",
        "arrayRank[value]",
        "value: A dense array or general brace value.",
        "This is structural array rank, not linear-algebra matrixRank.",
        "arrayRank[{{1, 2}, {3, 4}}]  ->  2\n"
        "arrayRank[{{1, 2}, {3}}]  ->  1"),
    HELP(Length,
        "Returns the outer element count of an array or brace.",
        "length[value]",
        "value: An array or general brace, including a ragged brace.",
        "length[{{1, 2}, {3}}]  ->  2"),
    HELP_NOTE(ArrayGet,
        "Extracts an element or prefix slice with zero-based indices.",
        "at[array, index1, index2, ...]",
        "array: A dense array or general brace.\n"
        "indices: One or more zero-based integers within the corresponding dimensions.",
        "A prefix shorter than the array rank returns the remaining subarray.",
        "at[{{1, 2}, {3, 4}}, 1]  ->  {3, 4}\n"
        "at[{{1, 2}, {3, 4}}, 1, 0]  ->  3"),
    HELP_NOTE(Reshape,
        "Changes array dimensions while preserving row-major element order.",
        "reshape[array, {d1, d2, ...}]",
        "array: A dense array or flat/general brace with the required element count.\n"
        "dimensions: A brace of non-negative integer sizes whose product matches the element count.",
        "reshape also preserves trailing shape after a zero-length leading dimension.",
        "reshape[{1, 2, 3, 4}, {2, 2}]  ->  {{1, 2}, {3, 4}}\n"
        "reshape[{}, {0, 3}]"),
    HELP(Identity,
        "Constructs an exact square identity matrix.",
        "identity[n]",
        "n: A non-negative integer matrix order.",
        "identity[2]  ->  {{1, 0}, {0, 1}}"),
    HELP(Zeros,
        "Constructs an exact rectangular zero matrix.",
        "zeros[rows, columns]",
        "rows, columns: Non-negative integer dimensions.",
        "zeros[2, 3]  ->  {{0, 0, 0}, {0, 0, 0}}\n"
        "zeros[0, 3]  ->  reshape[{}, {0, 3}]"),
    HELP(Trace,
        "Sums the main diagonal of a matrix.",
        "trace[matrix]",
        "matrix: A square rectangular matrix.",
        "trace[{{1, 2}, {3, 4}}]  ->  5"),
    HELP(Rows,
        "Returns the row count of a matrix.",
        "rows[matrix]",
        "matrix: A rank-2 rectangular array.",
        "rows[{{1, 2, 3}, {4, 5, 6}}]  ->  2"),
    HELP(Cols,
        "Returns the column count of a matrix.",
        "cols[matrix]",
        "matrix: A rank-2 rectangular array.",
        "cols[{{1, 2, 3}, {4, 5, 6}}]  ->  3"),
    HELP(Diag,
        "Extracts the main diagonal of a rectangular matrix.",
        "diag[matrix]",
        "matrix: A rank-2 rectangular array.",
        "diag[{{1, 2, 3}, {4, 5, 6}}]  ->  {1, 5}"),
    HELP(VectorAdd,
        "Adds two equal-length vectors element by element.",
        "vadd[a, b]",
        "a, b: Rank-1 arrays of the same length.",
        "vadd[{1, 2}, {3, 4}]  ->  {4, 6}"),
    HELP(VectorSubtract,
        "Subtracts two equal-length vectors element by element.",
        "vsub[a, b]",
        "a, b: Rank-1 arrays of the same length.",
        "vsub[{4, 6}, {1, 2}]  ->  {3, 4}"),
    HELP(VectorScale,
        "Multiplies every vector component by a scalar.",
        "vscalar[vector, scalar]",
        "vector: A rank-1 array.\n"
        "scalar: A scalar expression, not an array.",
        "vscalar[{1, 2}, 3]  ->  {3, 6}"),
    HELP(VectorCross,
        "Computes the exact three-dimensional cross product.",
        "vcross[a, b]",
        "a, b: Rank-1 arrays of exactly three components.",
        "vcross[{1, 0, 0}, {0, 1, 0}]  ->  {0, 0, 1}"),
    HELP_NOTE(VectorNorm,
        "Computes the Euclidean/Hermitian norm of a vector.",
        "norm[vector]",
        "vector: A rank-1 real or complex array.",
        "Complex components use conjugate products, so the result is a real magnitude.",
        "norm[{3, 4}]  ->  5\n"
        "norm[{3+4I}]  ->  5"),
    HELP(VectorManhattan,
        "Computes Manhattan (L1) distance between two real vectors.",
        "vmanhattan[a, b]",
        "a, b: Equal-length rank-1 real arrays.",
        "vmanhattan[{1, 2}, {4, 6}]  ->  7"),
    HELP(VectorEuclidean,
        "Computes Euclidean distance between two vectors.",
        "veuclidean[a, b]",
        "a, b: Equal-length rank-1 real or complex arrays.",
        "veuclidean[{1, 2}, {4, 6}]  ->  5"),
    HELP(VectorNormalize,
        "Divides a nonzero vector by its Hermitian norm.",
        "normalize[vector]",
        "vector: A nonzero rank-1 real or complex array.",
        "normalize[{3, 4}]  ->  {3/5, 4/5}"),
    HELP(VectorProject,
        "Projects the first vector onto the direction of the second.",
        "vproject[vector, onto]",
        "vector, onto: Equal-length rank-1 arrays; onto must be nonzero.",
        "vproject[{1, 2}, {0, 1}]  ->  {0, 2}"),
    HELP(VectorAngle,
        "Computes the angle between two nonzero vectors.",
        "vangle[a, b]",
        "a, b: Equal-length nonzero rank-1 vectors.",
        "vangle[{1, 0}, {0, 1}]  ->  Pi/2"),
    HELP(VectorReflect,
        "Reflects a vector across a hyperplane specified by its normal.",
        "vreflect[vector, normal]",
        "vector, normal: Equal-length rank-1 arrays; normal must be nonzero.",
        "vreflect[{1, 1}, {0, 1}]  ->  {1, -1}"),
    HELP(VectorReflectAxis,
        "Reflects a vector about an axis direction.",
        "vreflect_axis[vector, axis]",
        "vector, axis: Equal-length rank-1 arrays; axis must be nonzero.",
        "vreflect_axis[{1, 1}, {0, 1}]  ->  {-1, 1}"),
    HELP(VectorSum,
        "Sums all components of a vector.",
        "vsum[vector]",
        "vector: A rank-1 array.",
        "vsum[{1, 2, 3}]  ->  6"),

    HELP(Expm1,
        "Computes exp[x]-1 with certified stability near zero.",
        "expm1[x]",
        "x: A real or complex exact/certified numeric value.",
        "expm1[0]  ->  0\n"
        "N[expm1[1/10^30], 50]"),
    HELP(Log1p,
        "Computes the principal log[1+x] with certified stability near zero.",
        "log1p[x]",
        "x: A real or complex exact/certified value; x=-1 is a branch singularity.",
        "log1p[0]  ->  0\n"
        "N[log1p[1/10^30], 50]"),
    HELP_NOTE(Sinc,
        "Computes the cardinal sine sin[x]/x with the removable value filled at zero.",
        "sinc[x]",
        "x: An angle; Deg, Rad, or Grad overrides angleMode[].",
        "The denominator is the angle converted to radians.",
        "sinc[0]  ->  1\n"
        "sinc[Pi/2]  ->  2/Pi\n"
        "sinc[90 Deg]  ->  2/Pi"),
    HELP_NOTE(Cosc,
        "Computes the cardinal cosine (1-cos[x])/x with its limit at zero.",
        "cosc[x]",
        "x: An angle; Deg, Rad, or Grad overrides angleMode[].",
        "The denominator is the angle converted to radians.",
        "cosc[0]  ->  0\n"
        "N[cosc[Pi/3], 20]  ->  0.47746482927568600731"),
    HELP_NOTE(Tanc,
        "Computes tan[x]/x with its removable value at zero.",
        "tanc[x]",
        "x: An angle away from tangent poles; Deg, Rad, or Grad overrides angleMode[].",
        "The denominator is the angle converted to radians.",
        "tanc[0]  ->  1\n"
        "tanc[Pi/4]  ->  4/Pi"),
    HELP(Sinhc,
        "Computes sinh[x]/x with the removable value 1 at zero.",
        "sinhc[x]",
        "x: A real or complex exact/certified numeric value.",
        "sinhc[0]  ->  1"),
    HELP(Tanhc,
        "Computes tanh[x]/x with the removable value 1 at zero.",
        "tanhc[x]",
        "x: A real or complex exact/certified value away from complex tanh poles.",
        "tanhc[0]  ->  1"),
    HELP(Expc,
        "Computes (exp[x]-1)/x with the removable value 1 at zero.",
        "expc[x]",
        "x: A real or complex exact/certified numeric value.",
        "expc[0]  ->  1"),
    HELP(Log2,
        "Computes the principal base-2 logarithm.",
        "log2[x]",
        "x: A nonzero real or complex value; positive real input gives a real result.",
        "log2[8]  ->  3"),
    HELP(Log10,
        "Computes the principal base-10 logarithm.",
        "log10[x]",
        "x: A nonzero real or complex value; positive real input gives a real result.",
        "log10[1000]  ->  3"),
    HELP_NOTE(Gamma,
        "Computes Euler's gamma function exactly when possible.",
        "gamma[x]\n"
        "N[gamma[x], digits]",
        "x: A real or complex value that is not a non-positive integer pole.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "The certified numerical backend currently emphasizes real arguments.",
        "gamma[5]  ->  24\n"
        "gamma[1/2]  ->  sqrt[Pi]\n"
        "N[gamma[1/3], 20]  ->  2.6789385347077476337"),
    HELP_NOTE(LogGamma,
        "Computes log[abs[gamma[x]]] on the real axis.",
        "lgamma[x]",
        "x: A real value that is not a non-positive integer pole.",
        "This is a real log-magnitude function, kept distinct from a general complex LogGamma.",
        "lgamma[5]  ->  log[24]"),
    HELP_NOTE(LambertW,
        "Returns a branch of Lambert W, the inverse of w*exp[w].",
        "lambertw[z]\n"
        "lambertw[branch, z]",
        "z: An exact or certified numeric value.\n"
        "branch: An integer branch index; omission selects branch 0.",
        "Certified numerical evaluation currently covers the real branches 0 and -1 on their real domains.",
        "lambertw[0]  ->  0\n"
        "lambertw[E]  ->  1\n"
        "N[lambertw[-1, -1/10], 20]  ->  -3.5771520639572972184"),
    HELP(Erf,
        "Computes the Gaussian error function.",
        "erf[x]\n"
        "N[erf[x], digits]",
        "x: A real or complex value.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "erf[0]  ->  0\n"
        "N[erf[1], 20]  ->  0.84270079294971486934"),
    HELP(Erfc,
        "Computes the complementary error function 1-erf[x].",
        "erfc[x]\n"
        "N[erfc[x], digits]",
        "x: A real or complex value.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "erfc[0]  ->  1"),
    HELP_NOTE(FresnelC,
        "Computes the Fresnel cosine integral from 0 to x.",
        "fresnelc[x]\n"
        "N[fresnelc[x], digits]",
        "x: A real value for the current certified backend.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "The definition is integral_0^x cos[Pi*t^2/2] dt and never depends on angleMode[].",
        "fresnelc[0]  ->  0\n"
        "N[fresnelc[1], 20]  ->  0.77989340037682282947"),
    HELP_NOTE(FresnelS,
        "Computes the Fresnel sine integral from 0 to x.",
        "fresnels[x]\n"
        "N[fresnels[x], digits]",
        "x: A real value for the current certified backend.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "The definition is integral_0^x sin[Pi*t^2/2] dt and never depends on angleMode[].",
        "fresnels[0]  ->  0\n"
        "N[fresnels[1], 20]  ->  0.43825914739035476608"),
    HELP_NOTE(Hypergeometric1F1,
        "Computes Kummer's confluent hypergeometric function 1F1.",
        "hypergeometric1F1[a, b, z]\n"
        "N[hypergeometric1F1[a, b, z], digits]",
        "a, b: Exact parameters; b=0,-1,-2,... is generally a pole.\n"
        "z: The argument. The certified real backend currently accepts exact Rational a,b,z.",
        "Exact terminating series and safe degenerations are simplified.",
        "hypergeometric1F1[0, 3, 2]  ->  1\n"
        "hypergeometric1F1[2, 2, 1]  ->  E\n"
        "N[hypergeometric1F1[1/6, 7/6, 1], 20]"),
    HELP_NOTE(Hypergeometric2F1,
        "Computes the principal Gauss hypergeometric function 2F1.",
        "hypergeometric2F1[a, b, c, z]\n"
        "N[hypergeometric2F1[a, b, c, z], digits]",
        "a, b, c: Exact parameters; c=0,-1,-2,... is generally a pole.\n"
        "z: The argument. Certified real evaluation currently requires exact Rational parameters and |z|<1.",
        "Exact terminating series and safe degenerations are simplified.",
        "hypergeometric2F1[-2, 1, 3, 1/2]  ->  17/24\n"
        "hypergeometric2F1[0, 2, 3, x]  ->  1\n"
        "N[hypergeometric2F1[1/2, 1/2, 3/2, 1/4], 20]"),
    HELP_NOTE(EllipticF,
        "Computes Legendre's incomplete elliptic integral of the first kind.",
        "ellipticF[phi, m]\n"
        "N[ellipticF[phi, m], digits]",
        "phi: The amplitude, always interpreted in radians.\n"
        "m: The elliptic parameter (not the modulus). Certified real evaluation currently requires exact Rational phi and |m|<1.",
        "angleMode[] does not affect phi.",
        "ellipticF[x, 0]  ->  x\n"
        "N[ellipticF[1/2, 1/3], 20]  ->  0.5068477562654311092"),
    HELP_NOTE(EllipticE,
        "Computes Legendre's incomplete elliptic integral of the second kind.",
        "ellipticE[phi, m]\n"
        "N[ellipticE[phi, m], digits]",
        "phi: The amplitude, always interpreted in radians.\n"
        "m: The elliptic parameter. Certified real evaluation currently requires exact Rational phi and |m|<1.",
        "angleMode[] does not affect phi.",
        "ellipticE[x, 0]  ->  x\n"
        "N[ellipticE[1/2, 1/3], 20]  ->  0.49331536201475850521"),
    HELP_NOTE(EllipticPi,
        "Computes Legendre's incomplete elliptic integral of the third kind.",
        "ellipticPi[n, phi, m]\n"
        "N[ellipticPi[n, phi, m], digits]",
        "n: The characteristic.\n"
        "phi: The amplitude, always interpreted in radians.\n"
        "m: The elliptic parameter. Certified real evaluation currently requires exact Rational inputs with |n|<1 and |m|<1.",
        "Principal branches are used; angleMode[] does not affect phi.",
        "ellipticPi[0, x, 0]  ->  x\n"
        "N[ellipticPi[1/5, 1/2, 1/3], 20]"),
    HELP_NOTE(ExponentialIntegralEi,
        "Computes the principal exponential integral Ei.",
        "Ei[x]\n"
        "N[Ei[x], digits]",
        "x: An exact value; the certified backend covers supported real regions.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Ei has branch structure; unsupported certified regions remain symbolic. Real-axis limits at 0 and +/-Infinity are known exactly.",
        "N[Ei[1], 20]  ->  1.8951178163559367555\n"
        "D[Ei[x], x]  ->  exp[x]/x\n"
        "limit[Ei[x], x, 0]  ->  -Infinity\n"
        "limit[Ei[x], x, -Infinity]  ->  0"),
    HELP_NOTE(SineIntegralSi,
        "Computes the entire sine integral Si.",
        "Si[x]\n"
        "N[Si[x], digits]",
        "x: An exact value; the certified backend covers supported real regions.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Si is odd.",
        "Si[0]  ->  0\n"
        "Si[-1]  ->  -Si[1]\n"
        "N[Si[1], 20]  ->  0.94608307036718301494"),
    HELP_NOTE(CosineIntegralCi,
        "Computes the principal cosine integral Ci.",
        "Ci[x]\n"
        "N[Ci[x], digits]",
        "x: An exact value; the certified backend covers supported real regions.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "Ci has logarithmic branch structure. On the principal branch, Ci[x] -> I Pi along the negative real axis at -Infinity.",
        "N[Ci[1], 20]  ->  0.33740392290096813466\n"
        "D[Ci[x], x]  ->  cos[x]/x\n"
        "limit[Ci[x], x, 0]  ->  -Infinity\n"
        "limit[Ci[x], x, -Infinity]  ->  I Pi"),
    HELP_NOTE(LogarithmicIntegralLi,
        "Computes the principal logarithmic integral li.",
        "li[x]\n"
        "N[li[x], digits]",
        "x: An exact value away from the logarithmic singularity; certified evaluation covers supported real regions.\n"
        "digits: Optional significant digits supplied by an enclosing N.",
        "This is the offset-independent principal logarithmic integral with derivative 1/log[x]. li[0] is exactly 0 and li[1] is singular.",
        "li[0]  ->  0\n"
        "N[li[2], 20]  ->  1.0451637801174927848\n"
        "D[li[x], x]  ->  1/log[x]"),
    HELP_NOTE(Polylog,
        "Computes the principal polylogarithm Li_s[z].",
        "polylog[s, z]\n"
        "N[polylog[s, z], digits]",
        "s: The order; the certified backend currently requires supported exact parameters.\n"
        "z: The argument on the principal branch.",
        "Several low-order and special-point identities simplify exactly.",
        "polylog[0, z]  ->  z/(1-z)\n"
        "polylog[1, z]  ->  -log[1-z]\n"
        "polylog[2, 1]  ->  Pi^2/6\n"
        "N[polylog[2, 1/2], 20]"),
    HELP_NOTE(Beta,
        "Computes Euler's beta function for positive real arguments.",
        "beta[a, b]",
        "a, b: Positive real values.",
        "The implementation avoids unsafe unconditional Gamma-ratio expansion near cancellations.",
        "beta[2, 3]  ->  1/12\n"
        "beta[1/2, 1/2]  ->  Pi"),
    HELP_NOTE(BetaLog,
        "Computes the logarithm of Euler's beta function on its positive-real domain.",
        "betaln[a, b]",
        "a, b: Positive real values.",
        "This form is useful when beta itself would be very large or small.",
        "betaln[1/2, 1/2]  ->  log[Pi]"),
    HELP_NOTE(Zeta,
        "Computes the Riemann zeta function with exact special values where known.",
        "zeta[s]\n"
        "N[zeta[s], digits]",
        "s: A real or complex value; the current certified real backend covers s>1.",
        "General complex analytic continuation is not yet certified.",
        "zeta[0]  ->  -1/2\n"
        "zeta[-2]  ->  0\n"
        "zeta[2]  ->  Pi^2/6\n"
        "N[zeta[3], 20]"),
    HELP_NOTE(Digamma,
        "Computes the logarithmic derivative of gamma.",
        "digamma[x]\n"
        "N[digamma[x], digits]",
        "x: A real or complex value that is not a non-positive integer pole.",
        "digamma[x] is D[lgamma[x],x] on the real lgamma domain. Certified N covers the positive real axis and complex inputs that can be shifted into the right half-plane with a certified asymptotic remainder.",
        "N[digamma[1], 20]  ->  -0.57721566490153286061\n"
        "N[digamma[1+I], 20]"),
    HELP_NOTE(Trigamma,
        "Computes the derivative of digamma (the order-1 polygamma function).",
        "trigamma[x]\n"
        "N[trigamma[x], digits]",
        "x: A real or complex value that is not a non-positive integer pole.",
        "Positive integer arguments have exact reductions. Certified complex N uses recurrence plus an Euler-Maclaurin enclosure in the right half-plane.",
        "trigamma[2]  ->  Pi^2/6-1\n"
        "N[trigamma[1], 20]  ->  1.6449340668482264365\n"
        "N[trigamma[1+I], 20]"),
    HELP_NOTE(IncompleteBeta,
        "Computes the regularized incomplete beta function I_x(a,b).",
        "ibeta[a, b, x]\n"
        "N[ibeta[a, b, x], digits]",
        "a, b: Positive exact real parameters.\n"
        "x: A real value in [0,1]. The certified backend currently uses exact Rational a,b and certified real x.",
        "Positive-integer a,b with exact Rational x can reduce to a finite binomial sum.",
        "ibeta[1, 1, 1/4]  ->  1/4\n"
        "ibeta[2, 3, 1/2]  ->  11/16\n"
        "N[ibeta[1/3, 2/3, 1/4], 20]"),
    HELP(GeneralizedBinomial,
        "Computes a generalized binomial coefficient by a finite product.",
        "binom[x, n]",
        "x: An exact scalar.\n"
        "n: A non-negative integer order.",
        "binom[1/2, 2]  ->  -1/8"),
    HELP(FallingFactorial,
        "Computes x*(x-1)*...*(x-n+1) as an exact finite product.",
        "fallingfact[x, n]",
        "x: An exact scalar.\n"
        "n: A non-negative integer order.",
        "fallingfact[5, 3]  ->  60"),
    HELP(RisingFactorial,
        "Computes x*(x+1)*...*(x+n-1) as an exact finite product.",
        "risingfact[x, n]",
        "x: An exact scalar.\n"
        "n: A non-negative integer order.",
        "risingfact[5, 3]  ->  210"),

    HELP_NOTE(RandSeed,
        "Sets or refreshes the current session's pseudorandom seed.",
        "randSeed[]\n"
        "randSeed[seed]",
        "seed: Optional integer. With no argument, entropy selects and returns a reproducible integer seed.",
        "Random state belongs to each KernelSession and is not cryptographically secure.",
        "randSeed[42]  ->  42\n"
        "a := rand[]\n"
        "randSeed[42]  ->  42\n"
        "rand[] == a  ->  True"),
    HELP_NOTE(Rand,
        "Generates an exact uniform real sample on a 53-bit dyadic lattice.",
        "rand[]\n"
        "rand[upper]\n"
        "rand[lower, upper]",
        "upper: A non-negative real bound, giving [0,upper).\n"
        "lower, upper: Real bounds with lower <= upper, giving [lower,upper).",
        "The result is an exact Rational and consumes session RNG state.",
        "randSeed[42]  ->  42\n"
        "rand[]  ->  227930101193189/1125899906842624"),
    HELP_NOTE(RandInt,
        "Generates an unbiased uniform integer from an inclusive range.",
        "randint[]\n"
        "randint[bound]\n"
        "randint[lower, upper]",
        "With no input the set is {0,1}.\n"
        "bound >= 0 gives [0,bound]; bound < 0 gives [bound,0].\n"
        "lower, upper are inclusive arbitrary-size integer bounds with lower <= upper.",
        "Rejection sampling avoids modulo bias; the function consumes session RNG state.",
        "randint[1, 6]\n"
        "randint[-5]"),
    HELP(Choice,
        "Selects one element uniformly from supplied choices.",
        "choice[x1, x2, ...]\n"
        "choice[array]",
        "Provide one or more values directly, or one non-empty rank-1 array.",
        "choice[2, 3, 5, 7]\n"
        "choice[{2, 3, 5, 7}]"),
    HELP_NOTE(RandN,
        "Generates a normally distributed symbolic sample with Box-Muller.",
        "randn[]\n"
        "randn[mean]\n"
        "randn[mean, sigma]",
        "mean: Optional real center; default 0.\n"
        "sigma: Optional non-negative real standard deviation; default 1.",
        "Uses explicit radians, consumes session RNG state, and is not cryptographically secure.",
        "randn[5, 0]  ->  5\n"
        "N[randn[], 8]"),

    HELP(Sqrt,
        "Computes the principal square root exactly when possible.",
        "sqrt[x]",
        "x: An exact, symbolic, or certified real/complex value.",
        "sqrt[72]  ->  6sqrt[2]\n"
        "sqrt[-4]  ->  2I"),
    HELP(Abs,
        "Computes real absolute value or complex magnitude.",
        "abs[x]",
        "x: A real or complex exact/certified value.",
        "abs[-3]  ->  3\n"
        "abs[3+4I]  ->  5"),
    HELP(Sign,
        "Returns real sign, or z/abs[z] for a nonzero complex value.",
        "sign[x]",
        "x: A real or complex exact/certified value.",
        "sign[-3]  ->  -1\n"
        "sign[3+4I]  ->  3/5+4/5I"),
    HELP(Re,
        "Returns the real part of a complex value.",
        "re[z]",
        "z: A real or complex exact/certified value.",
        "re[3+4I]  ->  3"),
    HELP(Im,
        "Returns the imaginary part of a complex value.",
        "im[z]",
        "z: A real or complex exact/certified value.",
        "im[3+4I]  ->  4"),
    HELP(Conj,
        "Returns the complex conjugate.",
        "conj[z]",
        "z: A real or complex exact/certified value.",
        "conj[3+4I]  ->  3-4I"),

    HELP_NOTE(Sin,
        "Computes sine with exact special-angle simplification where possible.",
        "sin[x]",
        "x: A real or complex angle. A bare real uses angleMode[]; Deg, Rad, or Grad overrides it.",
        "The default session angle mode is Rad.",
        "sin[Pi/6]  ->  1/2\n"
        "sin[30 Deg]  ->  1/2\n"
        "sin[100 Grad]  ->  1"),
    HELP_NOTE(Cos,
        "Computes cosine with exact special-angle simplification where possible.",
        "cos[x]",
        "x: A real or complex angle. A bare real uses angleMode[]; Deg, Rad, or Grad overrides it.",
        "The default session angle mode is Rad.",
        "cos[Pi/3]  ->  1/2\n"
        "cos[60 Deg]  ->  1/2"),
    HELP_NOTE(Tan,
        "Computes tangent with exact special-angle simplification and pole checks.",
        "tan[x]",
        "x: A real or complex angle away from tangent poles. A bare real uses angleMode[]; an explicit unit overrides it.",
        "Exact poles produce a DomainError rather than a fabricated finite value.",
        "tan[Pi/4]  ->  1\n"
        "tan[50 Grad]  ->  1"),
    HELP_NOTE(Cot,
        "Computes cotangent as the reciprocal of tangent, with pole checks.",
        "cot[x]",
        "x: A real or complex angle away from sine zeros. A bare real uses angleMode[].",
        "Exact poles produce a DomainError.",
        "cot[Pi/4]  ->  1"),
    HELP_NOTE(Sec,
        "Computes secant as the reciprocal of cosine, with pole checks.",
        "sec[x]",
        "x: A real or complex angle away from cosine zeros. A bare real uses angleMode[].",
        "Exact poles produce a DomainError.",
        "sec[Pi/3]  ->  2"),
    HELP_NOTE(Csc,
        "Computes cosecant as the reciprocal of sine, with pole checks.",
        "csc[x]",
        "x: A real or complex angle away from sine zeros. A bare real uses angleMode[].",
        "Exact poles produce a DomainError.",
        "csc[Pi/6]  ->  2"),
    HELP_NOTE(Asin,
        "Computes the principal inverse sine.",
        "asin[x]",
        "x: A real or complex value. Real results are expressed in the current angle mode.",
        "The real principal range is [-quarter-turn, quarter-turn].",
        "asin[1/2]  ->  Pi/6\n"
        "angleMode[Deg]  ->  Deg\n"
        "asin[1/2]  ->  30"),
    HELP_NOTE(Acos,
        "Computes the principal inverse cosine.",
        "acos[x]",
        "x: A real or complex value. Real results are expressed in the current angle mode.",
        "The real principal range is [0, half-turn].",
        "acos[1/2]  ->  Pi/3\n"
        "angleMode[Deg]  ->  Deg\n"
        "acos[1/2]  ->  60"),
    HELP_NOTE(Atan,
        "Computes the principal inverse tangent.",
        "atan[x]",
        "x: A real or complex value. Real results are expressed in the current angle mode.",
        "The real principal range is (-quarter-turn, quarter-turn).",
        "atan[1]  ->  Pi/4\n"
        "angleMode[Deg]  ->  Deg\n"
        "atan[1]  ->  45"),
    HELP_NOTE(Atan2,
        "Computes the quadrant-aware angle of the Cartesian point (x,y).",
        "atan2[y, x]",
        "y: The vertical coordinate.\n"
        "x: The horizontal coordinate. Both must be real; the result uses angleMode[].",
        "The argument order is y first, x second.",
        "atan2[1, -1]  ->  3Pi/4\n"
        "angleMode[Deg]  ->  Deg\n"
        "atan2[1, -1]  ->  135"),
    HELP(Sinh,
        "Computes the hyperbolic sine on its entire complex domain.",
        "sinh[x]",
        "x: A real or complex exact/certified value; angleMode[] does not apply.",
        "sinh[0]  ->  0"),
    HELP(Cosh,
        "Computes the hyperbolic cosine on its entire complex domain.",
        "cosh[x]",
        "x: A real or complex exact/certified value; angleMode[] does not apply.",
        "cosh[0]  ->  1"),
    HELP(Tanh,
        "Computes the hyperbolic tangent on its principal meromorphic domain.",
        "tanh[x]",
        "x: A real or complex exact/certified value away from complex poles.",
        "tanh[0]  ->  0"),
    HELP(Asinh,
        "Computes the principal inverse hyperbolic sine.",
        "asinh[x]",
        "x: A real or complex exact/certified value; angleMode[] does not apply.",
        "asinh[0]  ->  0"),
    HELP(Acosh,
        "Computes the principal inverse hyperbolic cosine.",
        "acosh[x]",
        "x: A real or complex exact/certified value; real output requires x>=1.",
        "acosh[1]  ->  0"),
    HELP(Atanh,
        "Computes the principal inverse hyperbolic tangent.",
        "atanh[x]",
        "x: A real or complex value away from branch singularities at -1 and 1.",
        "atanh[0]  ->  0"),
    HELP(Csch,
        "Computes the reciprocal of hyperbolic sine.",
        "csch[x]",
        "x: A real or complex value away from zeros of sinh.",
        "N[csch[1], 16]"),
    HELP(Sech,
        "Computes the reciprocal of hyperbolic cosine.",
        "sech[x]",
        "x: A real or complex value away from complex zeros of cosh.",
        "sech[0]  ->  1"),
    HELP(Coth,
        "Computes the reciprocal of hyperbolic tangent.",
        "coth[x]",
        "x: A real or complex value away from zeros of tanh.",
        "N[coth[1], 16]"),
    HELP_NOTE(Arg,
        "Returns the principal argument of a nonzero complex value.",
        "arg[z]",
        "z: A nonzero real or complex exact/certified value.",
        "The result carries an explicit angle unit, so it is not reinterpreted by a later angleMode[] change.",
        "arg[-1]  ->  Pi Rad\n"
        "N[arg[-1], 20]  ->  3.1415926535897932385 Rad"),
    HELP_NOTE(Log,
        "Computes the principal natural logarithm or a logarithm in an explicit base.",
        "log[x]\n"
        "log[base, x]",
        "x: A nonzero real or complex value.\n"
        "base: A nonzero value other than 1; it is the first argument in the two-argument form.",
        "Negative real input uses the principal complex branch.",
        "log[E]  ->  1\n"
        "log[-1]  ->  I Pi\n"
        "log[10, 1000]  ->  3"),
    HELP(Exp,
        "Computes the entire complex exponential.",
        "exp[x]",
        "x: A real or complex exact/certified value.",
        "exp[0]  ->  1\n"
        "exp[1]  ->  E"),

    HELP_NOTE(NumericalApproximation,
        "Returns a certified decimal approximation while preserving unsupported symbolic parts.",
        "N[expression]\n"
        "N[expression, digits]",
        "expression: Held until the precision context is established; arrays are traversed recursively.\n"
        "digits: Positive significant decimal digits; default 16. This is not digits after the decimal point.",
        "Supported operations prove a containing interval and unique final rounding. Existing low-precision values cannot regain hidden guard digits through an outer N.",
        "N[Pi, 20]  ->  3.1415926535897932385\n"
        "N[sqrt[2], 30]  ->  1.41421356237309504880168872421\n"
        "N[x+Pi, 20]  ->  3.1415926535897932385+x"),
    HELP_NOTE(Precision,
        "Returns a lower bound on guaranteed relative decimal digits.",
        "precision[x]",
        "x: An exact expression or a certified decimal approximation.",
        "Exact values return Infinity. If the information interval contains zero, relative precision is 0.",
        "precision[1/3]  ->  Infinity\n"
        "precision[N[1/3, 20]]  ->  19"),
    HELP_NOTE(Accuracy,
        "Returns a lower bound on guaranteed absolute decimal digits.",
        "accuracy[x]",
        "x: An exact expression or a certified decimal approximation.",
        "Exact values return Infinity; approximate values use their InformationEnclosure, not hidden guard digits.",
        "accuracy[1/3]  ->  Infinity\n"
        "accuracy[N[1/3, 20]]  ->  20"),
    HELP_NOTE(Explain,
        "Reports metadata already carried by an evaluated value or registered function.",
        "explain[value]\n"
        "explain[value, \"internal\"]",
        "value: Evaluated normally before inspection.\n"
        "\"internal\": Optional development mode adding representation diagnostics.",
        "explain is lightweight: it does not compute determinant, rank, eigenvalues, or other derived properties. Internal property names are not compatibility-stable.",
        "explain[Pi]\n"
        "explain[{{1, 2}, {3, 4}}]\n"
        "explain[sin]"),
    HELP_NOTE(Rationalize,
        "Converts certified decimal information to a simple exact Rational.",
        "rationalize[x]\n"
        "rationalize[x, tolerance]",
        "x: A decimal approximation, array/expression containing approximations, or an already exact value.\n"
        "tolerance: Optional non-negative exact real radius around the displayed value.",
        "Without tolerance, chooses the smallest-denominator Rational in the InformationEnclosure. With tolerance 0, converts the displayed finite decimal exactly.",
        "rationalize[N[1/3, 20]]  ->  1/3\n"
        "rationalize[N[Pi, 20], 1/1000]  ->  201/64"),
    HELP_NOTE(Root,
        "Represents a certified exact algebraic root without forcing a radical expansion.",
        "root[{a0, a1, ..., an}, k]\n"
        "root[{a0, a1, ..., an}, k, Complex]",
        "coefficients: Exact Rational coefficients of a0+a1*x+...+an*x^n, in ascending power order.\n"
        "k: A 1-based distinct-root index. Real roots are increasing; Complex roots use a deterministic certified order.\n"
        "Complex: Optional domain selector; omission selects real roots.",
        "The defining polynomial is normalized to monic square-free form. N[root[...,k],digits] refines its certified isolating interval or disk.",
        "root[{-2, 0, 1}, 2]  ->  root[{-2, 0, 1}, 2]\n"
        "N[root[{-2, 0, 1}, 2], 30]  ->  1.41421356237309504880168872421\n"
        "N[root[{1, 0, 1}, 2, Complex], 30]  ->  I"),
    HELP_NOTE(Simplify,
        "Applies conservative exact simplification, optionally under assumptions.",
        "simplify[expression]\n"
        "simplify[expression, assumptions]",
        "expression: Any expression.\n"
        "assumptions: A predicate, array of predicates, or And-like condition set.",
        "Rewrites that would change principal branches or remove domain holes are not applied without sufficient assumptions.",
        "simplify[sin[x]^2+cos[x]^2]  ->  1\n"
        "simplify[sqrt[x^2], element[x, Real]]  ->  abs[x]\n"
        "simplify[sqrt[x^2], x >= 0]  ->  x"),
    HELP_NOTE(FullSimplify,
        "Searches a bounded set of equivalent forms for a lower-cost expression.",
        "fullSimplify[expression]\n"
        "fullSimplify[expression, assumptions]",
        "expression: Any expression.\n"
        "assumptions: A predicate, array of predicates, or And-like condition set.",
        "A shorter form is rejected when it changes the domain or principal-branch meaning.",
        "fullSimplify[x^2+2x+1]  ->  (1+x)^2\n"
        "fullSimplify[(x^2-1)/(x-1), x != 1]  ->  1+x"),
    HELP(Expand,
        "Expands products and non-negative integer powers.",
        "expand[expression]",
        "expression: An algebraic expression; expansion is bounded by resource limits.",
        "expand[(x+1)^3]  ->  x^3+3x^2+3x+1"),
    HELP(Factor,
        "Factors a polynomial expression over exact supported coefficients.",
        "factor[expression]",
        "expression: A polynomial or rational-algebraic expression in supported exact domains.",
        "factor[x^2-1]  ->  (x-1)(x+1)"),
    HELP(Collect,
        "Collects like powers of a selected symbol.",
        "collect[expression, variable]",
        "expression: An algebraic expression.\n"
        "variable: The symbol whose powers should be grouped.",
        "collect[a*x+b*x+c, x]  ->  c+(a+b)x"),
    HELP_NOTE(GroebnerBasis,
        "Computes a reduced Groebner basis of an exact rational polynomial ideal.",
        "groebnerBasis[{p1, ...}, {x1, ...}]\n"
        "groebnerBasis[{p1, ...}, {x1, ...}, order]",
        "p1, ...: Exact Rational-coefficient multivariate polynomials.\n"
        "x1, ...: The declared polynomial-ring variables.\n"
        "order: Optional Lex, GrLex, or GrevLex; default GrevLex.",
        "The implementation uses bounded Buchberger reduction over Q[x1,...,xn]. Approximate coefficients are rejected rather than silently rationalized.",
        "groebnerBasis[{x y-1, y^2-x}, {x, y}, Lex]  ->  {x-y^2, y^3-1}"),
    HELP_NOTE(PolynomialReduce,
        "Divides an exact multivariate polynomial by an ordered divisor list.",
        "polynomialReduce[p, {g1, ...}, {x1, ...}]\n"
        "polynomialReduce[p, {g1, ...}, {x1, ...}, order]",
        "p: Exact Rational-coefficient dividend.\n"
        "g1, ...: Exact Rational-coefficient divisors.\n"
        "x1, ...: Polynomial-ring variables.\n"
        "order: Optional Lex, GrLex, or GrevLex; default GrevLex.",
        "Returns {{q1,...}, r} such that p=sum(qi gi)+r under the selected multivariate division order.",
        "polynomialReduce[x^2+y^2, {x-y, y^2-1}, {x, y}, Lex]  ->  {{x+y, 2}, 2}"),
    HELP_NOTE(Solve,
        "Solves equations, inequalities, or equation systems while preserving unresolved cases.",
        "solve[equation, variable]\n"
        "solve[equation, variable, domainOrConstraint]\n"
        "solve[equation, domain]\n"
        "solve[{equation1, ...}, {variable1, ...}]",
        "equation: An equality, ordered inequality, or brace of equations.\n"
        "variable: A user symbol, or a brace of symbols for a system.\n"
        "domain: Integer, Rational, Real, or Complex; default equation domain is Complex. A constraint may be used instead.",
        "The domain-only shorthand infers a variable only when exactly one unknown user symbol exists. Unsupported families remain unresolved rather than being reported as empty.",
        "solve[x^2==1, x]  ->  {x==1, x==-1}\n"
        "solve[x^2+1==0, x, Real]  ->  {}\n"
        "solve[{2x+3y==5, x-2y==9}, {x, y}]  ->  {{x==37/7, y==-13/7}}\n"
        "solve[sin[x]==0, x, Real]  ->  {x==Pi k where k in Integer}"),

    HELP(Element,
        "Tests or states membership in a numeric domain.",
        "element[value, domain]",
        "value: An exact, symbolic, or certified value.\n"
        "domain: Integer, Rational, Real, or Complex.",
        "element[1/2, Integer]  ->  False\n"
        "element[Pi, Rational]  ->  False\n"
        "element[x, Real]  ->  element[x, Real]"),
    HELP_NOTE(If,
        "Selects exactly one branch after evaluating a Boolean condition.",
        "if[condition, trueExpression, falseExpression]",
        "condition: Must evaluate to True or False.\n"
        "trueExpression: Evaluated only when condition is True.\n"
        "falseExpression: Evaluated only when condition is False.",
        "if is a held special form; errors and random-number consumption in the unselected branch do not occur.",
        "if[True, 1, 1/0]  ->  1\n"
        "if[False, 1/0, 2]  ->  2"),
    HELP_NOTE(Cases,
        "Represents a scalar mathematical value by ordered conditional branches.",
        "cases[value1 if condition1; value2 if condition2; ...]\n"
        "cases[value if condition; defaultValue]",
        "value: A scalar expression; an unselected value is held.\n"
        "condition: A Boolean predicate.\n"
        "defaultValue: Optional final branch without a condition.",
        "Unlike if, cases is a first-class mathematical expression. Unknown conditions are preserved, while simplify can eliminate branches proven by assumptions. N approximates values without numericalizing predicates.",
        "cases[1/x if x != 0; 0 if x == 0]\n"
        "simplify[cases[1/x if x != 0; 0 if x == 0], x != 0]  ->  1/x"),
    HELP_NOTE(InputHistory,
        "Retrieves a previous input expression and reevaluates it in the current session.",
        "In[index]",
        "index: A nonzero integer. Positive values are absolute prompt numbers; negative values count prior input slots relative to the current input.",
        "In uses current definitions and RNG state. @, @@, @@@, ... abbreviate In[-1], In[-2], In[-3], ... .",
        "In[1]\n"
        "In[-1]\n"
        "N[@, 30]"),
    HELP_NOTE(OutputHistory,
        "Returns a stored successful output snapshot without reevaluation.",
        "Out[index]",
        "index: A nonzero integer. Positive values are absolute input numbers; negative values count successful outputs from the most recent.",
        "%, %%, %%%, ... abbreviate Out[-1], Out[-2], Out[-3], ... . Failed evaluations have no positive Out snapshot.",
        "Out[1]\n"
        "Out[-1]\n"
        "%"),
    HELP(Exit,
        "Leaves the interactive calculator.",
        "Exit[]",
        "No arguments.",
        "Exit[]"),
    HELP_NOTE(Clear,
        "Clears user definitions and all input/output history in the current session.",
        "Clear[]",
        "No arguments.",
        "The next interactive prompt returns to In [1]. Builtins and protected constants are not removed.",
        "Clear[]"),
    HELP(Definitions,
        "Lists current global user variables and user-function definitions.",
        "Defs[]",
        "No arguments.",
        "x := 3\n"
        "f[t] := t^2\n"
        "Defs[]"),
    HELP_NOTE(Undefine,
        "Removes one or more user variable or user-function definitions.",
        "UnDef[name1, name2, ...]",
        "name1, name2, ...: Unevaluated user symbols; builtins and protected predefined symbols are rejected.",
        "Returns the number of supplied names whose definitions were actually removed.",
        "x := 3\n"
        "f[t] := t^2\n"
        "UnDef[x, f]  ->  2"),
    HELP_NOTE(AngleMode,
        "Shows or changes the session's default angle unit.",
        "angleMode[]\n"
        "angleMode[Rad]\n"
        "angleMode[Deg]\n"
        "angleMode[Grad]",
        "unit: Optional Rad, Deg, or Grad. With no argument, reports the current mode.",
        "The default is Rad. An explicit unit suffix on a trigonometric argument always overrides this session setting.",
        "angleMode[]  ->  Rad\n"
        "angleMode[Deg]  ->  Deg\n"
        "sin[30]  ->  1/2"),
};

#undef HELP
#undef HELP_NOTE

struct ConstantHelpEntry final {
    std::string_view name;
    std::string_view summary;
    std::string_view notes;
    std::string_view examples;
};

constexpr ConstantHelpEntry constantHelpEntries[] = {
    {"Pi", "The exact circle constant: circumference divided by diameter.",
        "Pi is protected, real, positive, irrational, and transcendental.",
        "Pi\nN[Pi, 20]  ->  3.1415926535897932385\nsin[Pi/6]  ->  1/2"},
    {"E", "The exact base of the natural logarithm.",
        "E is protected, real, positive, irrational, and transcendental.",
        "E\nlog[E]  ->  1\nexp[1]  ->  E"},
    {"Phi", "The exact golden ratio (1+sqrt[5])/2.",
        "Phi is a protected positive real algebraic constant.",
        "Phi\nN[Phi, 20]  ->  1.6180339887498948482"},
    {"I", "The exact imaginary unit satisfying I^2=-1.",
        "I is protected and is lowered to an exact complex number during evaluation.",
        "I^2  ->  -1\nabs[3+4I]  ->  5"},
    {"True", "The Boolean truth value.",
        "True is protected and is accepted by held conditionals such as if.",
        "if[True, 1, 2]  ->  1"},
    {"False", "The Boolean false value.",
        "False is protected and is accepted by held conditionals such as if.",
        "if[False, 1, 2]  ->  2"},
    {"Integer", "The domain of exact integers.",
        "Use it with element, simplify assumptions, or solve domain selection.",
        "element[3, Integer]  ->  True\nsolve[x^2==4, x, Integer]"},
    {"Rational", "The domain of exact rational numbers.",
        "Every Integer is Rational; known irrational constants are not.",
        "element[1/2, Rational]  ->  True\nelement[Pi, Rational]  ->  False"},
    {"Real", "The real-number domain.",
        "Use it with element, assumptions, and solve to select real branches and ordered inequalities.",
        "element[Pi, Real]  ->  True\nsolve[x^2+1==0, x, Real]  ->  {}"},
    {"Complex", "The complex-number domain.",
        "Complex is the default ambient domain for equation solving.",
        "element[I, Complex]  ->  True\nsolve[x^2+1==0, x, Complex]  ->  {x==I, x==-I}"},
    {"Infinity", "A positive extended-real infinity sentinel.",
        "It is used by limits and by precision/accuracy for exact values. Undefined combinations such as Infinity-Infinity return Indeterminate.",
        "limit[1/x, x, 0, 1]  ->  Infinity\n"
        "Infinity/Infinity  ->  Indeterminate\n"
        "precision[Pi]  ->  Infinity"},
    {"ComplexInfinity", "An infinite magnitude whose real or complex direction is undetermined.",
        "A nonzero exact value divided by exact zero produces ComplexInfinity. It is protected and is distinct from Indeterminate.",
        "1/0  ->  ComplexInfinity\n"
        "I/0  ->  ComplexInfinity\n"
        "ComplexInfinity^0  ->  Indeterminate"},
    {"Indeterminate", "A protected result for a numerical value that is not unambiguously defined.",
        "Arithmetic and scalar mathematical functions propagate it. It is not a real or complex number and is not equal to itself.",
        "0/0  ->  Indeterminate\n"
        "0^0  ->  Indeterminate\n"
        "sin[Indeterminate]  ->  Indeterminate"},
    {"Rad", "The radian angle-unit symbol and default angle mode.",
        "Append it to an angle expression or pass it to angleMode. A full turn is 2Pi Rad.",
        "sin[Pi/6 Rad]  ->  1/2\nangleMode[Rad]  ->  Rad"},
    {"Deg", "The degree angle-unit symbol.",
        "Append it to an angle expression or pass it to angleMode. A full turn is 360 Deg.",
        "sin[30 Deg]  ->  1/2\nangleMode[Deg]  ->  Deg"},
    {"Grad", "The gradian angle-unit symbol.",
        "Append it to an angle expression or pass it to angleMode. A full turn is 400 Grad.",
        "sin[100 Grad]  ->  1\nangleMode[Grad]  ->  Grad"},
};

[[nodiscard]] std::string_view trimHelpWhitespace(std::string_view text) noexcept {
    while (!text.empty() && (text.front() == ' ' || text.front() == '\t'))
        text.remove_prefix(1);
    while (!text.empty() && (text.back() == ' ' || text.back() == '\t'))
        text.remove_suffix(1);
    return text;
}

[[nodiscard]] const FunctionHelpEntry* findFunctionHelp(BuiltinId id) noexcept {
    for (const FunctionHelpEntry& entry : functionHelpEntries) {
        if (entry.id == id)
            return &entry;
    }
    return nullptr;
}

[[nodiscard]] const ConstantHelpEntry* findConstantHelp(std::string_view name) noexcept {
    for (const ConstantHelpEntry& entry : constantHelpEntries) {
        if (entry.name == name)
            return &entry;
    }
    return nullptr;
}

void printIndentedBlock(
    std::string_view heading,
    std::string_view contents,
    std::ostream& output) {
    if (contents.empty())
        return;

    output << heading << ":\n";
    while (true) {
        const std::size_t newline = contents.find('\n');
        output << "  " << contents.substr(0, newline) << '\n';
        if (newline == std::string_view::npos)
            break;
        contents.remove_prefix(newline + 1);
    }
}

void printArgumentCount(
    const evaluation::BuiltinDefinition& definition,
    std::ostream& output) {
    output << "Arguments: ";
    if (definition.minimumArguments == definition.maximumArguments) {
        output << definition.minimumArguments;
    }
    else if (definition.maximumArguments == evaluation::BuiltinDefinition::unlimited) {
        output << definition.minimumArguments << " or more";
    }
    else {
        output << definition.minimumArguments << " to " << definition.maximumArguments;
    }
    output << '\n';
}

[[nodiscard]] std::string argumentName(std::size_t index) {
    constexpr std::string_view commonNames[] = {"x", "y", "z"};
    if (index < 3)
        return std::string{commonNames[index]};
    return "arg" + std::to_string(index + 1);
}

[[nodiscard]] std::string generatedUsage(
    std::string_view name,
    const evaluation::BuiltinDefinition& definition) {
    std::string result{name};
    result.push_back('[');
    for (std::size_t index = 0; index < definition.minimumArguments; ++index) {
        if (index != 0)
            result += ", ";
        result += argumentName(index);
    }
    if (definition.maximumArguments != definition.minimumArguments) {
        if (definition.minimumArguments != 0)
            result += ", ";
        result += "...";
    }
    result.push_back(']');
    return result;
}

void printFunctionIndex(
    const evaluation::BuiltinRegistry& registry,
    std::ostream& output) {
    const auto available = registry.sourceFunctionNames();
    std::vector<std::string> names{available.begin(), available.end()};
    std::sort(names.begin(), names.end());

    output << "Functions (canonical names and callable aliases):\n  ";
    std::size_t column = 2;
    for (std::size_t index = 0; index < names.size(); ++index) {
        const std::size_t needed = names[index].size() + (index == 0 ? 0 : 2);
        if (index != 0 && column + needed > 78) {
            output << "\n  ";
            column = 2;
        }
        else if (index != 0) {
            output << ", ";
            column += 2;
        }
        output << names[index];
        column += names[index].size();
    }
    output << "\nUse :help <function> for description, input rules, and examples.\n";
}

void printConstantIndex(std::ostream& output) {
    output << "Constants and predefined symbols:\n  ";
    for (std::size_t index = 0; index < std::size(constantHelpEntries); ++index) {
        if (index != 0)
            output << ", ";
        output << constantHelpEntries[index].name;
    }
    output << "\nUse :help <name> for details.\n";
}

void printReplHelpSummary(std::ostream& output) {
    output
        << "REPL help:\n"
        << "  :help <function|constant>  show description, input rules, and examples\n"
        << "  :help functions            list built-in function names\n"
        << "  :help constants            list predefined symbols\n"
        << "  :fix <0..1000>|off         set fixed display digits\n"
        << "  :status                    show session status\n"
        << "  Exit[]                     leave the calculator\n";
}

[[nodiscard]] std::size_t damerauLevenshtein(
    std::string_view left,
    std::string_view right) {
    const std::size_t width = right.size() + 1;
    std::vector<std::size_t> distances((left.size() + 1) * width);
    for (std::size_t i = 0; i <= left.size(); ++i)
        distances[i * width] = i;
    for (std::size_t j = 0; j <= right.size(); ++j)
        distances[j] = j;

    for (std::size_t i = 1; i <= left.size(); ++i) {
        for (std::size_t j = 1; j <= right.size(); ++j) {
            const std::size_t substitution = distances[(i - 1) * width + j - 1]
                + (left[i - 1] == right[j - 1] ? 0U : 1U);
            const std::size_t deletion = distances[(i - 1) * width + j] + 1;
            const std::size_t insertion = distances[i * width + j - 1] + 1;
            std::size_t best = std::min({substitution, deletion, insertion});
            if (i > 1 && j > 1 && left[i - 1] == right[j - 2]
                && left[i - 2] == right[j - 1]) {
                best = std::min(best, distances[(i - 2) * width + j - 2] + 1);
            }
            distances[i * width + j] = best;
        }
    }
    return distances[left.size() * width + right.size()];
}

[[nodiscard]] std::vector<std::string> helpTopicNames(
    const evaluation::BuiltinRegistry& registry) {
    const auto available = registry.sourceFunctionNames();
    std::vector<std::string> names{available.begin(), available.end()};
    names.reserve(names.size() + std::size(constantHelpEntries));
    for (const ConstantHelpEntry& entry : constantHelpEntries)
        names.emplace_back(entry.name);
    std::sort(names.begin(), names.end());
    return names;
}

[[nodiscard]] std::string suggestHelpTopic(
    std::string_view requested,
    const evaluation::BuiltinRegistry& registry) {
    if (requested.empty() || requested.size() > 128)
        return {};

    // These are discoverability shorthands, not callable aliases.
    constexpr std::pair<std::string_view, std::string_view> commonShorthands[] = {
        {"qr", "qrDecomposition"},
        {"lu", "luDecomposition"},
    };
    for (const auto& [shortName, fullName] : commonShorthands) {
        if (requested == shortName)
            return std::string{fullName};
    }

    const std::vector<std::string> names = helpTopicNames(registry);
    std::string prefixMatch;
    if (requested.size() >= 2) {
        for (const std::string& name : names) {
            if (name.starts_with(requested)
                && (prefixMatch.empty() || name.size() < prefixMatch.size())) {
                prefixMatch = name;
            }
        }
    }
    if (!prefixMatch.empty())
        return prefixMatch;

    std::string bestName;
    std::size_t bestDistance = std::numeric_limits<std::size_t>::max();
    for (const std::string& name : names) {
        const std::size_t distance = damerauLevenshtein(requested, name);
        if (distance < bestDistance
            || (distance == bestDistance && (bestName.empty() || name.size() < bestName.size()))
            || (distance == bestDistance && name.size() == bestName.size() && name < bestName)) {
            bestDistance = distance;
            bestName = name;
        }
    }

    const std::size_t threshold = requested.size() <= 4 ? 1
        : requested.size() <= 8 ? 2 : 3;
    return bestDistance <= threshold ? bestName : std::string{};
}

void printConstantHelp(const ConstantHelpEntry& entry, std::ostream& output) {
    output << entry.name << '\n' << entry.summary << '\n';
    printIndentedBlock("Usage", entry.name, output);
    printIndentedBlock("Inputs", "This is a predefined symbol; use it without brackets.", output);
    printIndentedBlock("Notes", entry.notes, output);
    printIndentedBlock("Examples", entry.examples, output);
}

} // namespace

bool handleReplHelpCommand(
    std::string_view line,
    const evaluation::BuiltinRegistry& registry,
    std::ostream& output) {
    line = trimHelpWhitespace(line);
    if (!line.starts_with(":help"))
        return false;
    if (line.size() > 5 && line[5] != ' ' && line[5] != '\t')
        return false;

    const std::string_view requested = trimHelpWhitespace(line.substr(5));
    if (requested.empty()) {
        printReplHelpSummary(output);
        return true;
    }
    if (requested == "functions") {
        printFunctionIndex(registry, output);
        return true;
    }
    if (requested == "constants") {
        printConstantIndex(output);
        return true;
    }
    if (requested == ":fix" || requested == "fix") {
        output << "Usage: :fix <0..1000>|off\n";
        return true;
    }
    if (requested == ":status" || requested == "status") {
        output << "Usage: :status\n";
        return true;
    }
    if (const ConstantHelpEntry* constant = findConstantHelp(requested)) {
        printConstantHelp(*constant, output);
        return true;
    }

    const evaluation::BuiltinDefinition* definition = registry.find(requested);
    if (!definition || !definition->sourceCallable) {
        output << "No help for '" << requested << "'.";
        const std::string suggestion = suggestHelpTopic(requested, registry);
        if (!suggestion.empty())
            output << " Did you mean '" << suggestion << "'?";
        output << "\nUse :help functions or :help constants to list available topics.\n";
        return true;
    }

    const expression::Symbol& canonicalSymbol = registry.symbol(definition->id);
    const evaluation::BuiltinDefinition* canonical = registry.find(canonicalSymbol);
    const std::string_view primaryName = canonical && canonical->sourceCallable
        ? canonicalSymbol.view()
        : requested;

    output << primaryName;
    if (requested != primaryName)
        output << " (alias: " << requested << ')';
    output << '\n';

    const FunctionHelpEntry* entry = findFunctionHelp(definition->id);
    if (!entry) {
        output << "Detailed help is not yet available for this registered builtin.\n";
        printIndentedBlock("Usage", generatedUsage(primaryName, *definition), output);
        printArgumentCount(*definition, output);
        return true;
    }

    output << entry->summary << '\n';
    printIndentedBlock("Usage", entry->usage, output);
    printArgumentCount(*definition, output);
    printIndentedBlock("Inputs", entry->inputs, output);
    printIndentedBlock("Notes", entry->notes, output);
    printIndentedBlock("Examples", entry->examples, output);
    return true;
}

} // namespace mmcal::cli
