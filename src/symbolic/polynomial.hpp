#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/assumption_set.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/rational.hpp"

#include <cstddef>
#include <optional>
#include <span>
#include <string_view>
#include <vector>

namespace mmcal::symbolic {

// Solverが高速に係数へアクセスするための一変数・有理係数多項式。
// coefficients[i] が x^i の係数。内部計算用なので表示順とは独立している。
class RationalPolynomial final {
public:
    RationalPolynomial();
    explicit RationalPolynomial(std::vector<numeric::Rational> coefficients);

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] std::size_t degree() const noexcept;
    [[nodiscard]] const numeric::Rational& coefficient(std::size_t exponent) const noexcept;
    [[nodiscard]] const std::vector<numeric::Rational>& coefficients() const noexcept;

private:
    std::vector<numeric::Rational> coefficients_;
    void normalize();
};


// 指定した1変数についてだけ多項式とみなし、それ以外の式を係数として保持する。
// 係数は x を含まない任意のExprでよく、sin[y], a+b, log[z] などもそのまま保持できる。
// これにより collect/factor/solve が「他の記号はすべて変数」という制約から独立できる。
class ExpressionPolynomial final {
public:
    explicit ExpressionPolynomial(std::vector<expression::Expr> coefficients);

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] std::size_t degree() const noexcept;
    [[nodiscard]] const expression::Expr& coefficient(std::size_t exponent) const noexcept;
    [[nodiscard]] std::span<const expression::Expr> coefficients() const noexcept;

private:
    std::vector<expression::Expr> coefficients_;
    void normalize();
};

struct ExpressionPolynomialConversionResult final {
    ExpressionPolynomial polynomial;
    // 元の式が定義されるために必要な条件。solverが代数変形で分母等を
    // 消しても、この条件をSolutionSetへ運ぶために変換結果と一体で保持する。
    mathematics::AssumptionSet domainConditions;
};

// 多変数多項式の単項式。factorはSymbol名順に正規化され、指数0は保持しない。
struct MonomialFactor final {
    expression::Symbol variable;
    std::size_t exponent = 0;

    [[nodiscard]] bool operator==(const MonomialFactor&) const = default;
};

class Monomial final {
public:
    Monomial() = default;
    explicit Monomial(std::vector<MonomialFactor> factors);

    [[nodiscard]] bool isOne() const noexcept;
    [[nodiscard]] std::size_t totalDegree() const noexcept;
    [[nodiscard]] std::size_t exponentOf(const expression::Symbol& variable) const noexcept;
    [[nodiscard]] std::span<const MonomialFactor> factors() const noexcept;
    [[nodiscard]] Monomial without(const expression::Symbol& variable) const;
    [[nodiscard]] bool operator==(const Monomial&) const = default;

private:
    std::vector<MonomialFactor> factors_;
};

enum class MonomialOrder {
    Lex,
    GrLex,
    GrevLex
};

[[nodiscard]] std::optional<MonomialOrder> parseMonomialOrder(std::string_view name) noexcept;
[[nodiscard]] std::string_view monomialOrderName(MonomialOrder order) noexcept;

// Q[x1,...,xn] の変数順序とterm orderを一体で保持する。
// Gröbner算法では「変数集合」と「leading termの意味」を暗黙global stateにしない。
class PolynomialRing final {
public:
    PolynomialRing(
        std::vector<expression::Symbol> variables,
        MonomialOrder order = MonomialOrder::GrevLex);

    [[nodiscard]] std::span<const expression::Symbol> variables() const noexcept;
    [[nodiscard]] MonomialOrder order() const noexcept;
    [[nodiscard]] bool contains(const expression::Symbol& variable) const noexcept;
    [[nodiscard]] bool contains(const Monomial& monomial) const noexcept;
    // -1: lhs<rhs, 0: equal, +1: lhs>rhs。
    [[nodiscard]] int compare(const Monomial& lhs, const Monomial& rhs) const;

private:
    std::vector<expression::Symbol> variables_;
    MonomialOrder order_ = MonomialOrder::GrevLex;
};

[[nodiscard]] bool monomialDivides(const Monomial& divisor, const Monomial& dividend) noexcept;
[[nodiscard]] Monomial multiplyMonomials(const Monomial& lhs, const Monomial& rhs);
[[nodiscard]] Monomial leastCommonMultiple(const Monomial& lhs, const Monomial& rhs);
[[nodiscard]] std::optional<Monomial> divideMonomials(
    const Monomial& dividend,
    const Monomial& divisor);

struct PolynomialTerm final {
    Monomial monomial;
    numeric::Rational coefficient;

    [[nodiscard]] bool operator==(const PolynomialTerm&) const = default;
};

// Expand / collect / factor / 多変数solve が共有する多変数・有理係数表現。
// termsは同じ単項式を持たず、係数0の項を持たない。
class MultivariateRationalPolynomial final {
public:
    MultivariateRationalPolynomial();
    explicit MultivariateRationalPolynomial(std::vector<PolynomialTerm> terms);

    [[nodiscard]] bool isZero() const noexcept;
    [[nodiscard]] std::size_t termCount() const noexcept;
    [[nodiscard]] std::size_t totalDegree() const noexcept;
    [[nodiscard]] std::size_t degree(const expression::Symbol& variable) const noexcept;
    [[nodiscard]] std::span<const PolynomialTerm> terms() const noexcept;
    [[nodiscard]] std::vector<expression::Symbol> variables() const;
    [[nodiscard]] bool belongsTo(const PolynomialRing& ring) const noexcept;
    [[nodiscard]] std::optional<PolynomialTerm> leadingTerm(const PolynomialRing& ring) const;

private:
    std::vector<PolynomialTerm> terms_;
    void normalize();
};

struct PolynomialConversionOptions final {
    // 巨大な式展開を多項式変換の副作用として起こさないための防壁。
    std::size_t maximumDegree = 4096;
    std::size_t maximumTerms = 4096;
};

[[nodiscard]] MultivariateRationalPolynomial negatePolynomial(
    const MultivariateRationalPolynomial& value);
[[nodiscard]] MultivariateRationalPolynomial addPolynomials(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs);
[[nodiscard]] MultivariateRationalPolynomial subtractPolynomials(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs);
[[nodiscard]] std::optional<MultivariateRationalPolynomial> multiplyPolynomials(
    const MultivariateRationalPolynomial& lhs,
    const MultivariateRationalPolynomial& rhs,
    PolynomialConversionOptions options = {});
[[nodiscard]] MultivariateRationalPolynomial multiplyPolynomialByTerm(
    const MultivariateRationalPolynomial& polynomial,
    const PolynomialTerm& term);
[[nodiscard]] MultivariateRationalPolynomial monicPolynomial(
    const MultivariateRationalPolynomial& polynomial,
    const PolynomialRing& ring);

[[nodiscard]] std::optional<MultivariateRationalPolynomial> toMultivariateRationalPolynomial(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    PolynomialConversionOptions options = {});

[[nodiscard]] std::optional<MultivariateRationalPolynomial> toMultivariateRationalPolynomial(
    const expression::Expr& expression,
    const PolynomialRing& ring,
    const evaluation::BuiltinRegistry& builtins,
    PolynomialConversionOptions options = {});

[[nodiscard]] std::optional<RationalPolynomial> toRationalPolynomial(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    PolynomialConversionOptions options = {});

// variable を含まない部分式を「係数」として扱う一変数多項式へ変換する。
// 係数の四則演算はSimplifierでexactに正規化するため、MathRegistry/AngleSemanticsを受け取る。
[[nodiscard]] std::optional<ExpressionPolynomial> toExpressionPolynomial(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    PolynomialConversionOptions options = {});

// solver向け。記号分母やlog/arg等の定義条件を捨てずに変換する。
[[nodiscard]] std::optional<ExpressionPolynomialConversionResult>
toExpressionPolynomialWithConditions(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    PolynomialConversionOptions options = {});

// solverのsymbolic coefficient解析向け。式が係数として扱えるかを判定し、
// Divide/Log/Arg/負整数Powerなどの定義条件を返す。
[[nodiscard]] std::optional<mathematics::AssumptionSet> scalarExpressionDomainConditions(
    const expression::Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics);

[[nodiscard]] expression::Expr expressionPolynomialToCollectedExpr(
    const ExpressionPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

// 指定Symbolが式木のどこかに自由に現れるかを構造的に調べる。
// builtin headはSymbol式ではないため、引数だけを辿る。
[[nodiscard]] bool containsSymbol(
    const expression::Expr& expression,
    const expression::Symbol& symbol);

// 通常の数式表示に近い降冪順でExprへ戻す。
[[nodiscard]] expression::Expr polynomialToExpandedExpr(
    const RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins);

[[nodiscard]] expression::Expr polynomialToExpandedExpr(
    const MultivariateRationalPolynomial& polynomial,
    const evaluation::BuiltinRegistry& builtins);

// 指定変数の各次数で係数をまとめる。係数は他の変数を含む多項式Exprになり得る。
[[nodiscard]] expression::Expr polynomialToCollectedExpr(
    const MultivariateRationalPolynomial& polynomial,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins);

[[nodiscard]] numeric::Rational evaluatePolynomial(
    const RationalPolynomial& polynomial,
    const numeric::Rational& value);

// (x-root) でexactに割り切れる場合だけ商を返す。
[[nodiscard]] std::optional<RationalPolynomial> divideByLinearFactor(
    const RationalPolynomial& polynomial,
    const numeric::Rational& root);

struct RationalRootSearchOptions final {
    // Rational Root Theoremの約数列挙は巨大整数では高価なので有限の防壁を置く。
    // 防壁に達した場合 complete=false として「根なし」とは断定しない。
    std::size_t maximumTrialDivisor = 1'000'000;
    std::size_t maximumCandidates = 200'000;
};

struct RationalRootSearchResult final {
    std::optional<numeric::Rational> root;
    bool complete = false;
};

[[nodiscard]] RationalRootSearchResult findRationalRoot(
    const RationalPolynomial& polynomial,
    RationalRootSearchOptions options = {});

struct RationalPolynomialFactorOptions final {
    // 一般Q[x]はsquare-free -> finite-field factorization -> Hensel ->
    // Zassenhausを主経路とする。degreeだけの一律gateにせず，各phaseが
    // 実際に消費するresourceへ個別の防壁を置く。
    std::size_t maximumFiniteFieldDegree = 512;
    std::size_t maximumBerlekampMatrixEntries = 4096;
    std::size_t maximumPrimeTrials = 16;
    std::size_t maximumGoodPrimes = 8;
    // 2 prime以上を交差後，候補subsetがこの数以下ならliftへ進む。
    std::size_t preferredRecombinationCandidates = 1024;
    std::size_t maximumCantorZassenhausAttempts = 64;
    std::size_t maximumLiftBits = 4096;
    std::size_t maximumHenselWorkingBits = 8192;
    std::size_t maximumCombinations = 200'000;

    // modular factorが多くsubset列挙が膨らむ場合だけ，CLD latticeで
    // rational factorのblockを推定する。LLLはheuristic acceleratorであり，
    // 候補は必ず元のZ[x]でexact divisionする。
    std::size_t minimumCldModularFactors = 8;
    std::size_t maximumCldLiftBits = 4096;
    std::size_t maximumCldLattices = 4;
    std::size_t maximumCldCandidates = 256;
    std::size_t maximumCldCoefficientColumns = 4;
    std::size_t maximumLllRank = 24;
    std::size_t maximumLllColumns = 192;
    std::size_t maximumLllIntermediateBits = 8192;
    std::size_t maximumLllSwaps = 128;
    std::size_t maximumLllSizeReductions = 512;

    // exact sparse preprocessing用。候補は必ずQ[x]で割り切れることを
    // 検証するため，budget到達は「未分解」であって誤因子にはならない。
    std::size_t maximumSparseTerms = 24;
    std::size_t maximumSparseCandidates = 8192;

    // bounded Kronecker fallback用。
    std::size_t maximumKroneckerDegree = 64;
    std::size_t maximumSampleRadius = 24;
    std::size_t maximumDivisorsPerSample = 4096;
};

struct RationalPolynomialFactorization final {
    // polynomial == scalar * product(factors)。各factorはmonicかつ非定数である。
    numeric::Rational scalar;
    std::vector<RationalPolynomial> factors;
    // trueならmodular/Kroneckerを含む全leafについてQ上既約まで証明済み。
    // falseでもscalar*product(factors)は常に入力とexactに一致し，証明済みの
    // partial factorsを捨てずに返す。
    bool complete = false;
};

// 一変数Q[x]をbounded modular factorizationでexactに分解する。
// budget外ではKronecker補間へfallbackする。数値近似rootは用いず，
// 候補は必ずQ[x]のexact divisionで検証する。
[[nodiscard]] RationalPolynomialFactorization factorRationalPolynomialOverQ(
    const RationalPolynomial& polynomial,
    RationalPolynomialFactorOptions options = {});

// 式中の通常Symbolを構造的に収集する。builtin headは引数ではないため対象外。
[[nodiscard]] std::vector<expression::Symbol> collectSymbols(
    const expression::Expr& expression);

} // namespace mmcal::symbolic
