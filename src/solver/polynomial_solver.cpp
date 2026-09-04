// 多項式方程式solver
#include "polynomial_solver.hpp"

#include "mathematics/exact_algebra.hpp"
#include "mathematics/math_ids.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/predicate.hpp"
#include "mathematics/relation_builtin.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "solver/solve_constraints.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/groebner.hpp"
#include "symbolic/substitution.hpp"
#include "symbolic/algebraic_expression.hpp"
#include "symbolic/algebraic_number.hpp"
#include "symbolic/algebra_transforms.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using mathematics::RelationKind;

[[nodiscard]] Rational rational(std::int64_t numerator, std::int64_t denominator = 1) {
    return Rational{BigInt{numerator}, BigInt{denominator}};
}

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return expression.isCall()
        && expression.asCall().head.sameIdentity(builtins.symbol(id));
}

[[nodiscard]] Expr simplify(
    Expr expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplification::Simplifier{}.simplify(
        expression,
        simplification::SimplificationContext{builtins, mathematics, angles});
}

[[nodiscard]] SolutionBranch branch(
    const expression::Symbol& variable,
    Expr value,
    std::optional<std::size_t> multiplicity = std::nullopt) {
    return SolutionBranch{{SolutionBinding{variable, std::move(value)}}, {}, multiplicity, {}, std::nullopt};
}


[[nodiscard]] Expr zeroExpr() { return Expr{Number{BigInt{0}}}; }

[[nodiscard]] mathematics::Predicate equalZero(const Expr& value) {
    return mathematics::relation(mathematics::RelationKind::Equal, value, zeroExpr());
}

[[nodiscard]] mathematics::Predicate notEqualZero(const Expr& value) {
    return mathematics::relation(mathematics::RelationKind::NotEqual, value, zeroExpr());
}

[[nodiscard]] mathematics::TruthValue proveZero(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    const mathematics::AssumptionSet assumptions;
    return mathematics::KnowledgeContext{builtins, mathematics, assumptions}.prove(equalZero(value));
}

[[nodiscard]] mathematics::TruthValue proveNonZero(
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    const mathematics::AssumptionSet assumptions;
    return mathematics::KnowledgeContext{builtins, mathematics, assumptions}.prove(notEqualZero(value));
}

[[nodiscard]] mathematics::AssumptionSet conditions(
    std::initializer_list<mathematics::Predicate> predicates) {
    return mathematics::AssumptionSet{std::vector<mathematics::Predicate>{predicates}};
}

[[nodiscard]] mathematics::AssumptionSet appendConditions(
    const mathematics::AssumptionSet& base,
    std::initializer_list<mathematics::Predicate> predicates) {
    std::vector<mathematics::Predicate> result(base.predicates().begin(), base.predicates().end());
    result.insert(result.end(), predicates.begin(), predicates.end());
    return mathematics::AssumptionSet{std::move(result)};
}

[[nodiscard]] Expr symbolicQuadraticDiscriminant(
    const Expr& a,
    const Expr& b,
    const Expr& c,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr bSquared = Expr::call(
        builtins.symbol(BuiltinId::Power), {b, Expr{Number{BigInt{2}}}});
    Expr fourAC = Expr::call(
        builtins.symbol(BuiltinId::Multiply), {Expr{Number{BigInt{4}}}, a, c});
    return simplify(
        Expr::call(builtins.symbol(BuiltinId::Subtract), {std::move(bSquared), std::move(fourAC)}),
        builtins, mathematics, angles);
}

[[nodiscard]] Expr symbolicQuadraticRoot(
    const Expr& a,
    const Expr& b,
    const Expr& discriminant,
    bool plus,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr minusB = simplify(
        Expr::call(builtins.symbol(BuiltinId::Negate), {b}),
        builtins, mathematics, angles);
    Expr root = simplify(
        Expr::call(builtins.symbol(BuiltinId::Sqrt), {discriminant}),
        builtins, mathematics, angles);
    Expr numerator = Expr::call(
        builtins.symbol(plus ? BuiltinId::Add : BuiltinId::Subtract),
        {std::move(minusB), std::move(root)});
    Expr denominator = Expr::call(
        builtins.symbol(BuiltinId::Multiply), {Expr{Number{BigInt{2}}}, a});
    return simplify(
        Expr::call(builtins.symbol(BuiltinId::Divide),
            {std::move(numerator), std::move(denominator)}),
        builtins, mathematics, angles);
}

[[nodiscard]] Expr symbolicRepeatedQuadraticRoot(
    const Expr& a,
    const Expr& b,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr minusB = simplify(
        Expr::call(builtins.symbol(BuiltinId::Negate), {b}),
        builtins, mathematics, angles);
    Expr denominator = Expr::call(
        builtins.symbol(BuiltinId::Multiply), {Expr{Number{BigInt{2}}}, a});
    return simplify(
        Expr::call(builtins.symbol(BuiltinId::Divide),
            {std::move(minusB), std::move(denominator)}),
        builtins, mathematics, angles);
}

[[nodiscard]] Expr symbolicLinearRoot(
    const Expr& a,
    const Expr& b,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (b.isNumber() && b.asNumber().isZero())
        return zeroExpr();
    Expr minusB = simplify(
        Expr::call(builtins.symbol(BuiltinId::Negate), {b}),
        builtins, mathematics, angles);
    return simplify(
        Expr::call(builtins.symbol(BuiltinId::Divide), {std::move(minusB), a}),
        builtins, mathematics, angles);
}

[[nodiscard]] SolutionCase finiteCase(
    mathematics::AssumptionSet caseConditions,
    std::vector<SolutionBranch> branches) {
    return SolutionCase{
        std::move(caseConditions), SolutionSetKind::Finite, std::move(branches)};
}

[[nodiscard]] SolutionCase universalCase(mathematics::AssumptionSet caseConditions) {
    return SolutionCase{
        std::move(caseConditions), SolutionSetKind::Universal, {}};
}

[[nodiscard]] SolutionCase emptyCase(mathematics::AssumptionSet caseConditions) {
    return SolutionCase{
        std::move(caseConditions), SolutionSetKind::Empty, {}};
}

void appendDegenerateLinearCases(
    std::vector<SolutionCase>& cases,
    const mathematics::AssumptionSet& baseConditions,
    const Expr& b,
    const Expr& c,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto bZero = proveZero(b, builtins, mathematics);
    const auto bNonZero = proveNonZero(b, builtins, mathematics);
    if (bNonZero == mathematics::TruthValue::True) {
        cases.push_back(finiteCase(
            baseConditions,
            {branch(variable, symbolicLinearRoot(b, c, builtins, mathematics, angles))}));
        return;
    }
    if (bZero == mathematics::TruthValue::True) {
        const auto cZero = proveZero(c, builtins, mathematics);
        const auto cNonZero = proveNonZero(c, builtins, mathematics);
        if (cZero == mathematics::TruthValue::True) {
            cases.push_back(universalCase(baseConditions));
            return;
        }
        if (cNonZero == mathematics::TruthValue::True) {
            cases.push_back(emptyCase(baseConditions));
            return;
        }
        cases.push_back(universalCase(appendConditions(baseConditions, {equalZero(c)})));
        cases.push_back(emptyCase(appendConditions(baseConditions, {notEqualZero(c)})));
        return;
    }

    cases.push_back(finiteCase(
        appendConditions(baseConditions, {notEqualZero(b)}),
        {branch(variable, symbolicLinearRoot(b, c, builtins, mathematics, angles))}));
    const mathematics::AssumptionSet bIsZero = appendConditions(baseConditions, {equalZero(b)});
    const auto cZero = proveZero(c, builtins, mathematics);
    const auto cNonZero = proveNonZero(c, builtins, mathematics);
    if (cZero == mathematics::TruthValue::True)
        cases.push_back(universalCase(bIsZero));
    else if (cNonZero == mathematics::TruthValue::True)
        cases.push_back(emptyCase(bIsZero));
    else {
        cases.push_back(universalCase(appendConditions(bIsZero, {equalZero(c)})));
        cases.push_back(emptyCase(appendConditions(bIsZero, {notEqualZero(c)})));
    }
}

[[nodiscard]] SolutionSet solveExpressionPolynomial(
    const symbolic::ExpressionPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Complex}};
    if (polynomial.degree() > 2)
        return SolutionSet::unresolved(variables);

    const Expr c = polynomial.coefficient(0);
    if (polynomial.degree() == 0) {
        const auto cZero = proveZero(c, builtins, mathematics);
        if (cZero == mathematics::TruthValue::True)
            return SolutionSet::universal(variables);
        if (proveNonZero(c, builtins, mathematics) == mathematics::TruthValue::True)
            return SolutionSet::empty(variables);
        return SolutionSet::conditional(
            variables,
            {
                universalCase(conditions({equalZero(c)})),
                emptyCase(conditions({notEqualZero(c)}))
            });
    }

    const Expr b = polynomial.coefficient(1);
    if (polynomial.degree() == 1) {
        const auto bNonZero = proveNonZero(b, builtins, mathematics);
        if (bNonZero == mathematics::TruthValue::True)
            return SolutionSet::finite(
                variables,
                {branch(variable, symbolicLinearRoot(b, c, builtins, mathematics, angles))});

        std::vector<SolutionCase> cases;
        appendDegenerateLinearCases(
            cases, {}, b, c, variable, builtins, mathematics, angles);
        if (cases.size() == 1) {
            const SolutionCase& only = cases.front();
            if (only.outcome == SolutionSetKind::Finite)
                return SolutionSet::finite(variables, only.branches);
            if (only.outcome == SolutionSetKind::Universal)
                return SolutionSet::universal(variables, only.conditions);
            if (only.outcome == SolutionSetKind::Empty)
                return SolutionSet::empty(variables);
        }
        return SolutionSet::conditional(variables, std::move(cases));
    }

    const Expr a = polynomial.coefficient(2);
    const Expr discriminant = symbolicQuadraticDiscriminant(
        a, b, c, builtins, mathematics, angles);
    const auto aNonZero = proveNonZero(a, builtins, mathematics);
    const auto dZero = proveZero(discriminant, builtins, mathematics);
    const auto dNonZero = proveNonZero(discriminant, builtins, mathematics);

    auto quadraticBranches = [&](bool repeated) {
        if (repeated) {
            return std::vector<SolutionBranch>{branch(
                variable,
                symbolicRepeatedQuadraticRoot(a, b, builtins, mathematics, angles),
                std::size_t{2})};
        }
        return std::vector<SolutionBranch>{
            branch(variable, symbolicQuadraticRoot(
                a, b, discriminant, true, builtins, mathematics, angles)),
            branch(variable, symbolicQuadraticRoot(
                a, b, discriminant, false, builtins, mathematics, angles))};
    };

    if (aNonZero == mathematics::TruthValue::True) {
        if (dZero == mathematics::TruthValue::True)
            return SolutionSet::finite(variables, quadraticBranches(true));
        if (dNonZero == mathematics::TruthValue::True)
            return SolutionSet::finite(variables, quadraticBranches(false));
        return SolutionSet::conditional(
            variables,
            {
                finiteCase(conditions({notEqualZero(discriminant)}), quadraticBranches(false)),
                finiteCase(conditions({equalZero(discriminant)}), quadraticBranches(true))
            });
    }

    std::vector<SolutionCase> cases;
    mathematics::AssumptionSet aIsNonZero = conditions({notEqualZero(a)});
    if (dZero == mathematics::TruthValue::True)
        cases.push_back(finiteCase(aIsNonZero, quadraticBranches(true)));
    else if (dNonZero == mathematics::TruthValue::True)
        cases.push_back(finiteCase(aIsNonZero, quadraticBranches(false)));
    else {
        cases.push_back(finiteCase(
            appendConditions(aIsNonZero, {notEqualZero(discriminant)}),
            quadraticBranches(false)));
        cases.push_back(finiteCase(
            appendConditions(aIsNonZero, {equalZero(discriminant)}),
            quadraticBranches(true)));
    }

    appendDegenerateLinearCases(
        cases,
        conditions({equalZero(a)}),
        b, c, variable, builtins, mathematics, angles);
    return SolutionSet::conditional(variables, std::move(cases));
}

[[nodiscard]] Expr sqrtExpr(
    const Rational& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplify(
        Expr::call(builtins.symbol(BuiltinId::Sqrt), {Expr{Number{value}}}),
        builtins,
        mathematics,
        angles);
}

[[nodiscard]] Expr quadraticRoot(
    const Rational& a,
    const Rational& b,
    const Rational& discriminant,
    bool plus,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Expr minusB{Number{-b}};
    const Expr root = sqrtExpr(discriminant, builtins, mathematics, angles);
    Expr numerator = Expr::call(
        builtins.symbol(plus ? BuiltinId::Add : BuiltinId::Subtract),
        {minusB, root});
    Expr result = Expr::call(
        builtins.symbol(BuiltinId::Divide),
        {std::move(numerator), Expr{Number{rational(2) * a}}});
    return simplify(std::move(result), builtins, mathematics, angles);
}

[[nodiscard]] std::vector<SolutionBranch> solveQuadratic(
    const symbolic::RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const Rational a = polynomial.coefficient(2);
    const Rational b = polynomial.coefficient(1);
    const Rational c = polynomial.coefficient(0);
    const Rational discriminant = b * b - rational(4) * a * c;

    if (discriminant.isZero()) {
        Expr root{Number{(-b) / (rational(2) * a)}};
        return {branch(variable, std::move(root), std::size_t{2})};
    }

    return {
        branch(variable, quadraticRoot(a, b, discriminant, true, builtins, mathematics, angles)),
        branch(variable, quadraticRoot(a, b, discriminant, false, builtins, mathematics, angles))
    };
}

[[nodiscard]] std::optional<BigInt> exactPositiveIntegerNthRoot(
    const BigInt& value,
    std::uint64_t degree) {
    if (value.isNegative() || degree == 0)
        return std::nullopt;
    if (value.isZero() || value == BigInt{1} || degree == 1)
        return value;

    const std::size_t bitGuess = (value.bitLength() + static_cast<std::size_t>(degree) - 1)
        / static_cast<std::size_t>(degree) + 1;
    BigInt low{0};
    BigInt high = BigInt{1} << bitGuess;
    const BigInt one{1};
    while (high - low > one) {
        const BigInt middle = (low + high) / BigInt{2};
        const BigInt powered = numeric::pow(middle, degree);
        if (powered <= value)
            low = middle;
        else
            high = middle;
    }
    return numeric::pow(low, degree) == value
        ? std::optional<BigInt>{std::move(low)}
        : std::nullopt;
}

[[nodiscard]] std::optional<Rational> perfectPositiveRationalNthRoot(
    const Rational& value,
    std::uint64_t degree) {
    if (value < rational(0))
        return std::nullopt;
    const auto numerator = exactPositiveIntegerNthRoot(value.numerator(), degree);
    const auto denominator = exactPositiveIntegerNthRoot(value.denominator(), degree);
    if (!numerator || !denominator)
        return std::nullopt;
    return Rational{*numerator, *denominator};
}

[[nodiscard]] const expression::Symbol& piSymbol(const mathematics::MathRegistry& mathematics) {
    const auto* definition = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!definition)
        throw std::logic_error("Pi is not registered");
    return definition->symbol;
}

[[nodiscard]] Expr imaginaryUnitExpr() {
    return Expr{Number::complex(
        numeric::RealNumber{BigInt{0}},
        numeric::RealNumber{BigInt{1}})};
}

[[nodiscard]] Expr principalMagnitudeRoot(
    const Rational& magnitude,
    std::size_t degree,
    const evaluation::BuiltinRegistry& builtins) {
    if (const auto exact = perfectPositiveRationalNthRoot(
        magnitude, static_cast<std::uint64_t>(degree)))
        return Expr{Number{*exact}};
    return Expr::call(
        builtins.symbol(BuiltinId::Power),
        {Expr{Number{magnitude}}, Expr{Number{Rational{BigInt{1}, BigInt::parse(std::to_string(degree))}}}});
}

[[nodiscard]] std::vector<SolutionBranch> solveBinomial(
    const symbolic::RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::size_t degree = polynomial.degree();
    const Rational rhs = -polynomial.coefficient(0) / polynomial.coefficient(degree);
    if (rhs.isZero())
        return {branch(variable, Expr{Number{BigInt{0}}}, degree)};

    const bool negative = rhs < rational(0);
    const Rational magnitude = negative ? -rhs : rhs;
    const Expr radial = principalMagnitudeRoot(magnitude, degree, builtins);
    const Expr iPi = Expr::call(
        builtins.symbol(BuiltinId::Multiply),
        {imaginaryUnitExpr(), Expr{piSymbol(mathematics)}});

    std::vector<SolutionBranch> result;
    result.reserve(degree);
    for (std::size_t k = 0; k < degree; ++k) {
        const BigInt numerator = BigInt::parse(std::to_string(2 * k + (negative ? 1 : 0)));
        const BigInt denominator = BigInt::parse(std::to_string(degree));
        const Rational phaseCoefficient{numerator, denominator};
        Expr exponent = mathematics::scaleExactExpression(phaseCoefficient, iPi, builtins);
        Expr phase = simplify(
            Expr::call(builtins.symbol(BuiltinId::Exp), {std::move(exponent)}),
            builtins,
            mathematics,
            angles);
        Expr root = simplify(
            Expr::call(builtins.symbol(BuiltinId::Multiply), {radial, std::move(phase)}),
            builtins,
            mathematics,
            angles);
        result.push_back(branch(variable, std::move(root)));
    }

    // exactな実根が得られた場合は先頭へ置く。
    // 集合の数学的意味は順序に依存しないが、x^3=-8 -> {-2, ...} のような読みやすい表示を維持する。
    std::stable_partition(result.begin(), result.end(), [](const SolutionBranch& item) {
        return !item.bindings.empty()
            && item.bindings.front().value.isNumber()
            && item.bindings.front().value.asNumber().isReal();
    });
    return result;
}

[[nodiscard]] bool isBinomial(const symbolic::RationalPolynomial& polynomial) {
    if (polynomial.degree() < 1 || polynomial.coefficient(polynomial.degree()).isZero())
        return false;
    for (std::size_t exponent = 1; exponent < polynomial.degree(); ++exponent)
        if (!polynomial.coefficient(exponent).isZero())
            return false;
    return true;
}

[[nodiscard]] std::optional<std::vector<SolutionBranch>> solvePerfectCubeBinomial(
    const symbolic::RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (polynomial.degree() != 3 || !isBinomial(polynomial))
        return std::nullopt;

    const Rational rhs = -polynomial.coefficient(0) / polynomial.coefficient(3);
    const bool negative = rhs < rational(0);
    const Rational magnitude = negative ? -rhs : rhs;
    const auto positiveRoot = perfectPositiveRationalNthRoot(magnitude, 3);
    if (!positiveRoot)
        return std::nullopt;
    const Rational realCubeRoot = negative ? -*positiveRoot : *positiveRoot;

    Expr sqrt3 = sqrtExpr(rational(3), builtins, mathematics, angles);
    const Expr iSqrt3 = Expr::call(
        builtins.symbol(BuiltinId::Multiply), {imaginaryUnitExpr(), sqrt3});
    const Rational realPartCoefficient = -realCubeRoot / rational(2);
    const Rational imaginaryCoefficient = realCubeRoot / rational(2);
    const Expr realPart{Number{realPartCoefficient}};
    const Expr imaginaryPart = mathematics::scaleExactExpression(
        imaginaryCoefficient, iSqrt3, builtins);

    Expr second = simplify(
        Expr::call(builtins.symbol(BuiltinId::Add), {realPart, imaginaryPart}),
        builtins, mathematics, angles);
    Expr third = simplify(
        Expr::call(builtins.symbol(BuiltinId::Subtract), {realPart, imaginaryPart}),
        builtins, mathematics, angles);
    return std::vector<SolutionBranch>{
        branch(variable, Expr{Number{realCubeRoot}}),
        branch(variable, std::move(second)),
        branch(variable, std::move(third))};
}

[[nodiscard]] std::optional<std::vector<SolutionBranch>> solveByRationalDeflation(
    symbolic::RationalPolynomial polynomial,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<SolutionBranch> branches;

    while (polynomial.degree() > 2) {
        const symbolic::RationalRootSearchResult search = symbolic::findRationalRoot(polynomial);
        if (!search.root)
            break;
        const Rational root = *search.root;
        std::size_t multiplicity = 0;
        while (true) {
            auto quotient = symbolic::divideByLinearFactor(polynomial, root);
            if (!quotient)
                break;
            polynomial = std::move(*quotient);
            ++multiplicity;
            if (polynomial.degree() == 0)
                break;
        }
        branches.push_back(branch(variable, Expr{Number{root}}, multiplicity));
    }

    if (polynomial.degree() == 0)
        return branches;
    if (polynomial.degree() == 1) {
        const Rational root = -polynomial.coefficient(0) / polynomial.coefficient(1);
        branches.push_back(branch(variable, Expr{Number{root}}));
        return branches;
    }
    if (polynomial.degree() == 2) {
        auto quadratic = solveQuadratic(polynomial, variable, builtins, mathematics, angles);
        branches.insert(
            branches.end(),
            std::make_move_iterator(quadratic.begin()),
            std::make_move_iterator(quadratic.end()));
        return branches;
    }
    if (isBinomial(polynomial) && polynomial.degree() <= 256) {
        auto binomial = solveBinomial(polynomial, variable, builtins, mathematics, angles);
        branches.insert(
            branches.end(),
            std::make_move_iterator(binomial.begin()),
            std::make_move_iterator(binomial.end()));
        return branches;
    }
    return std::nullopt;
}


[[nodiscard]] std::optional<RelationKind> relationKindOf(
    const Expr& relation,
    const evaluation::BuiltinRegistry& builtins) {
    if (!relation.isCall() || relation.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto* definition = builtins.find(relation.asCall().head);
    if (!definition)
        return std::nullopt;
    return mathematics::relationKindForBuiltin(definition->id);
}


struct RationalFunctionPolynomial final {
    symbolic::RationalPolynomial numerator;
    symbolic::RationalPolynomial denominator;
};

[[nodiscard]] symbolic::RationalPolynomial polynomialOne() {
    return symbolic::RationalPolynomial{{rational(1)}};
}

[[nodiscard]] symbolic::RationalPolynomial addPolynomials(
    const symbolic::RationalPolynomial& lhs,
    const symbolic::RationalPolynomial& rhs) {
    const std::size_t size = std::max(lhs.coefficients().size(), rhs.coefficients().size());
    std::vector<Rational> coefficients(size, rational(0));
    for (std::size_t i = 0; i < lhs.coefficients().size(); ++i)
        coefficients[i] += lhs.coefficient(i);
    for (std::size_t i = 0; i < rhs.coefficients().size(); ++i)
        coefficients[i] += rhs.coefficient(i);
    return symbolic::RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] symbolic::RationalPolynomial negatePolynomial(
    const symbolic::RationalPolynomial& value) {
    std::vector<Rational> coefficients(value.coefficients().begin(), value.coefficients().end());
    for (Rational& coefficient : coefficients)
        coefficient = -coefficient;
    return symbolic::RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] symbolic::RationalPolynomial multiplyPolynomials(
    const symbolic::RationalPolynomial& lhs,
    const symbolic::RationalPolynomial& rhs) {
    if (lhs.isZero() || rhs.isZero())
        return symbolic::RationalPolynomial{};
    std::vector<Rational> coefficients(
        lhs.degree() + rhs.degree() + 1, rational(0));
    for (std::size_t i = 0; i <= lhs.degree(); ++i)
        for (std::size_t j = 0; j <= rhs.degree(); ++j)
            coefficients[i + j] += lhs.coefficient(i) * rhs.coefficient(j);
    return symbolic::RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] symbolic::RationalPolynomial powerPolynomial(
    symbolic::RationalPolynomial base,
    std::uint64_t exponent) {
    symbolic::RationalPolynomial result = polynomialOne();
    while (exponent != 0) {
        if ((exponent & 1U) != 0U)
            result = multiplyPolynomials(result, base);
        exponent >>= 1U;
        if (exponent != 0)
            base = multiplyPolynomials(base, base);
    }
    return result;
}

[[nodiscard]] std::optional<std::int64_t> smallExactInteger(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    try {
        return std::stoll(expression.asNumber().asReal().asInteger().toString());
    }
    catch (...) {
        return std::nullopt;
    }
}

[[nodiscard]] std::optional<RationalFunctionPolynomial> toRationalFunctionPolynomial(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    if (const auto polynomial = symbolic::toRationalPolynomial(expression, variable, builtins))
        return RationalFunctionPolynomial{*polynomial, polynomialOne()};

    if (!expression.isCall())
        return std::nullopt;
    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;

    if (definition->id == BuiltinId::Negate && arguments.size() == 1) {
        auto inner = toRationalFunctionPolynomial(arguments[0], variable, builtins);
        if (!inner)
            return std::nullopt;
        inner->numerator = negatePolynomial(inner->numerator);
        return inner;
    }

    if ((definition->id == BuiltinId::Add || definition->id == BuiltinId::Multiply)
        && !arguments.empty()) {
        auto result = toRationalFunctionPolynomial(arguments[0], variable, builtins);
        if (!result)
            return std::nullopt;
        for (std::size_t i = 1; i < arguments.size(); ++i) {
            auto rhs = toRationalFunctionPolynomial(arguments[i], variable, builtins);
            if (!rhs)
                return std::nullopt;
            if (definition->id == BuiltinId::Add) {
                const auto numerator = addPolynomials(
                    multiplyPolynomials(result->numerator, rhs->denominator),
                    multiplyPolynomials(rhs->numerator, result->denominator));
                result->denominator = multiplyPolynomials(result->denominator, rhs->denominator);
                result->numerator = numerator;
            }
            else {
                result->numerator = multiplyPolynomials(result->numerator, rhs->numerator);
                result->denominator = multiplyPolynomials(result->denominator, rhs->denominator);
            }
        }
        return result;
    }

    if ((definition->id == BuiltinId::Subtract || definition->id == BuiltinId::Divide)
        && arguments.size() == 2) {
        auto lhs = toRationalFunctionPolynomial(arguments[0], variable, builtins);
        auto rhs = toRationalFunctionPolynomial(arguments[1], variable, builtins);
        if (!lhs || !rhs)
            return std::nullopt;
        if (definition->id == BuiltinId::Subtract) {
            lhs->numerator = addPolynomials(
                multiplyPolynomials(lhs->numerator, rhs->denominator),
                negatePolynomial(multiplyPolynomials(rhs->numerator, lhs->denominator)));
            lhs->denominator = multiplyPolynomials(lhs->denominator, rhs->denominator);
            return lhs;
        }
        if (rhs->numerator.isZero())
            return std::nullopt;
        lhs->numerator = multiplyPolynomials(lhs->numerator, rhs->denominator);
        lhs->denominator = multiplyPolynomials(lhs->denominator, rhs->numerator);
        return lhs;
    }

    if (definition->id == BuiltinId::Power && arguments.size() == 2) {
        const auto exponent = smallExactInteger(arguments[1]);
        if (!exponent || *exponent < -64 || *exponent > 64)
            return std::nullopt;
        auto base = toRationalFunctionPolynomial(arguments[0], variable, builtins);
        if (!base)
            return std::nullopt;
        const auto magnitude = static_cast<std::uint64_t>(
            *exponent < 0 ? -*exponent : *exponent);
        if (*exponent < 0) {
            if (base->numerator.isZero())
                return std::nullopt;
            std::swap(base->numerator, base->denominator);
        }
        base->numerator = powerPolynomial(base->numerator, magnitude);
        base->denominator = powerPolynomial(base->denominator, magnitude);
        return base;
    }

    return std::nullopt;
}

[[nodiscard]] bool acceptsSign(RelationKind relation, int sign) noexcept {
    switch (relation) {
    case RelationKind::Less:
    case RelationKind::LessEqual:
        return sign < 0;
    case RelationKind::Greater:
    case RelationKind::GreaterEqual:
        return sign > 0;
    case RelationKind::Equal:
    case RelationKind::NotEqual:
        return false;
    }
    return false;
}

[[nodiscard]] bool acceptsZero(RelationKind relation) noexcept {
    return relation == RelationKind::LessEqual || relation == RelationKind::GreaterEqual;
}

struct OrderedRealRoot final {
    Expr value;
    std::size_t multiplicity = 1;
};

struct RealRangePiece final {
    std::optional<Expr> lower;
    bool lowerInclusive = false;
    std::optional<Expr> upper;
    bool upperInclusive = false;
};

void mergeRangePieces(std::vector<RealRangePiece>& pieces) {
    if (pieces.empty())
        return;

    std::vector<RealRangePiece> merged;
    merged.reserve(pieces.size());
    for (RealRangePiece piece : pieces) {
        if (merged.empty()) {
            merged.push_back(std::move(piece));
            continue;
        }

        RealRangePiece& previous = merged.back();
        // piecesは既に数直線順に生成される。境界が同じで、少なくとも片方がその点を
        // 含むならunionに穴はない。open interval + singleton + open intervalも
        // この規則を2回適用すれば1区間へまとまる。
        if (previous.upper && piece.lower
            && *previous.upper == *piece.lower
            && (previous.upperInclusive || piece.lowerInclusive)) {
            previous.upper = std::move(piece.upper);
            previous.upperInclusive = piece.upperInclusive;
            continue;
        }
        merged.push_back(std::move(piece));
    }
    pieces = std::move(merged);
}

[[nodiscard]] SolutionBranch realRangeBranch(
    const expression::Symbol& variable,
    RealRangePiece piece) {
    mathematics::AssumptionSet predicates;
    if (piece.lower) {
        predicates.add(mathematics::relation(
            piece.lowerInclusive ? RelationKind::GreaterEqual : RelationKind::Greater,
            Expr{variable}, *piece.lower));
    }
    if (piece.upper) {
        predicates.add(mathematics::relation(
            piece.upperInclusive ? RelationKind::LessEqual : RelationKind::Less,
            Expr{variable}, *piece.upper));
    }
    return SolutionBranch{
        {},
        std::move(predicates),
        std::nullopt,
        {SolverVariable{variable, mathematics::NumericDomain::Real}},
        std::nullopt};
}

[[nodiscard]] SolutionSet solutionFromSignChart(
    const expression::Symbol& variable,
    std::vector<OrderedRealRoot> roots,
    int leadingSign,
    RelationKind relation) {
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Real}};
    if (leadingSign == 0)
        return acceptsZero(relation)
            ? SolutionSet::universal(variables)
            : SolutionSet::empty(variables);

    const std::size_t n = roots.size();
    std::vector<int> intervalSigns(n + 1, leadingSign);
    // 右端(+infinity側)ではleading coefficientの符号。根を左へ跨ぐたび、
    // odd multiplicityなら符号反転、even multiplicityなら維持する。
    for (std::size_t i = n; i > 0; --i) {
        intervalSigns[i - 1] = intervalSigns[i];
        if ((roots[i - 1].multiplicity & 1U) != 0U)
            intervalSigns[i - 1] = -intervalSigns[i - 1];
    }

    std::vector<RealRangePiece> pieces;
    pieces.reserve(2 * n + 1);
    for (std::size_t interval = 0; interval <= n; ++interval) {
        if (!acceptsSign(relation, intervalSigns[interval]))
            continue;
        RealRangePiece piece;
        if (interval != 0)
            piece.lower = roots[interval - 1].value;
        if (interval != n)
            piece.upper = roots[interval].value;
        pieces.push_back(std::move(piece));
    }

    if (acceptsZero(relation)) {
        for (const OrderedRealRoot& root : roots) {
            pieces.push_back(RealRangePiece{root.value, true, root.value, true});
        }
    }

    // interval, singletonを数直線順へ戻す。roots自体は昇順と分かっているが、
    // Expr endpointを一般比較する必要はないよう、roots順からrankを求める。
    auto endpointRank = [&](const std::optional<Expr>& endpoint, bool upper) -> std::size_t {
        if (!endpoint)
            return upper ? 2 * n + 1 : 0;
        for (std::size_t i = 0; i < n; ++i)
            if (*endpoint == roots[i].value)
                return 2 * i + 1;
        return 2 * n + 2;
    };
    std::stable_sort(pieces.begin(), pieces.end(), [&](const RealRangePiece& lhs, const RealRangePiece& rhs) {
        const std::size_t l = endpointRank(lhs.lower, false);
        const std::size_t r = endpointRank(rhs.lower, false);
        if (l != r)
            return l < r;
        // 同じlowerならsingleton/閉区間をopen intervalより先へ置く。
        return lhs.lowerInclusive && !rhs.lowerInclusive;
    });
    mergeRangePieces(pieces);

    if (pieces.empty())
        return SolutionSet::empty(variables);
    if (pieces.size() == 1 && !pieces.front().lower && !pieces.front().upper)
        return SolutionSet::universal(variables);

    std::vector<SolutionBranch> branches;
    branches.reserve(pieces.size());
    for (RealRangePiece& piece : pieces) {
        if (piece.lower && piece.upper
            && *piece.lower == *piece.upper
            && piece.lowerInclusive && piece.upperInclusive) {
            branches.push_back(branch(variable, *piece.lower));
        }
        else {
            branches.push_back(realRangeBranch(variable, std::move(piece)));
        }
    }
    return SolutionSet::finite(variables, std::move(branches));
}

[[nodiscard]] int rationalSign(const Rational& value) noexcept {
    if (value < rational(0))
        return -1;
    if (value > rational(0))
        return 1;
    return 0;
}

[[nodiscard]] int compareRationalWithQuadraticRoot(
    const Rational& value,
    const Rational& a,
    const Rational& b,
    const Rational& c,
    bool lowerRoot) {
    const Rational vertex = -b / (rational(2) * a);
    const Rational q = a * value * value + b * value + c;
    if (q.isZero())
        return 0;

    const bool sameAsLeading = rationalSign(q) == rationalSign(a);
    if (lowerRoot) {
        if (value >= vertex)
            return 1;
        return sameAsLeading ? -1 : 1;
    }

    if (value <= vertex)
        return -1;
    return sameAsLeading ? 1 : -1;
}

[[nodiscard]] std::optional<std::vector<OrderedRealRoot>>
realRootsAfterRationalDeflation(
    symbolic::RationalPolynomial polynomial,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<std::pair<Rational, std::size_t>> rationalRoots;

    while (polynomial.degree() > 2) {
        const symbolic::RationalRootSearchResult search = symbolic::findRationalRoot(polynomial);
        if (!search.root)
            return std::nullopt;
        const Rational root = *search.root;
        std::size_t multiplicity = 0;
        while (auto quotient = symbolic::divideByLinearFactor(polynomial, root)) {
            polynomial = std::move(*quotient);
            ++multiplicity;
            if (polynomial.degree() == 0)
                break;
        }
        if (multiplicity == 0)
            return std::nullopt;
        rationalRoots.emplace_back(root, multiplicity);
    }

    if (polynomial.degree() == 1) {
        rationalRoots.emplace_back(
            -polynomial.coefficient(0) / polynomial.coefficient(1), 1);
        polynomial = symbolic::RationalPolynomial{{rational(1)}};
    }

    std::sort(rationalRoots.begin(), rationalRoots.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });

    std::vector<OrderedRealRoot> result;
    result.reserve(rationalRoots.size() + 2);

    if (polynomial.degree() != 2) {
        for (const auto& [root, multiplicity] : rationalRoots)
            result.push_back(OrderedRealRoot{Expr{Number{root}}, multiplicity});
        return result;
    }

    const Rational a = polynomial.coefficient(2);
    const Rational b = polynomial.coefficient(1);
    const Rational c = polynomial.coefficient(0);
    const Rational discriminant = b * b - rational(4) * a * c;
    if (discriminant < rational(0)) {
        for (const auto& [root, multiplicity] : rationalRoots)
            result.push_back(OrderedRealRoot{Expr{Number{root}}, multiplicity});
        return result;
    }
    if (discriminant.isZero()) {
        rationalRoots.emplace_back(-b / (rational(2) * a), 2);
        std::sort(rationalRoots.begin(), rationalRoots.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });
        for (const auto& [root, multiplicity] : rationalRoots)
            result.push_back(OrderedRealRoot{Expr{Number{root}}, multiplicity});
        return result;
    }

    const Expr rootDiscriminant = sqrtExpr(discriminant, builtins, mathematics, angles);
    if (rootDiscriminant.isNumber() && rootDiscriminant.asNumber().isReal()) {
        const Rational d = rootDiscriminant.asNumber().asReal().toRational();
        Rational first = (-b - d) / (rational(2) * a);
        Rational second = (-b + d) / (rational(2) * a);
        if (second < first)
            std::swap(first, second);
        rationalRoots.emplace_back(std::move(first), 1);
        rationalRoots.emplace_back(std::move(second), 1);
        std::sort(rationalRoots.begin(), rationalRoots.end(), [](const auto& lhs, const auto& rhs) {
            return lhs.first < rhs.first;
        });
        for (const auto& [root, multiplicity] : rationalRoots)
            result.push_back(OrderedRealRoot{Expr{Number{root}}, multiplicity});
        return result;
    }

    Expr lower = quadraticRoot(a, b, discriminant, false, builtins, mathematics, angles);
    Expr upper = quadraticRoot(a, b, discriminant, true, builtins, mathematics, angles);
    if (a < rational(0))
        std::swap(lower, upper);

    bool lowerInserted = false;
    bool upperInserted = false;
    for (const auto& [root, multiplicity] : rationalRoots) {
        if (!lowerInserted
            && compareRationalWithQuadraticRoot(root, a, b, c, true) > 0) {
            result.push_back(OrderedRealRoot{lower, 1});
            lowerInserted = true;
        }
        if (!upperInserted
            && compareRationalWithQuadraticRoot(root, a, b, c, false) > 0) {
            if (!lowerInserted) {
                result.push_back(OrderedRealRoot{lower, 1});
                lowerInserted = true;
            }
            result.push_back(OrderedRealRoot{upper, 1});
            upperInserted = true;
        }
        result.push_back(OrderedRealRoot{Expr{Number{root}}, multiplicity});
    }
    if (!lowerInserted)
        result.push_back(OrderedRealRoot{std::move(lower), 1});
    if (!upperInserted)
        result.push_back(OrderedRealRoot{std::move(upper), 1});
    return result;
}


[[nodiscard]] std::optional<std::vector<std::pair<Rational, std::size_t>>>
rationalRealRootsWithoutIrrationalResidual(symbolic::RationalPolynomial polynomial) {
    std::vector<std::pair<Rational, std::size_t>> roots;
    while (polynomial.degree() > 2) {
        const symbolic::RationalRootSearchResult search = symbolic::findRationalRoot(polynomial);
        if (!search.root)
            return std::nullopt;
        const Rational root = *search.root;
        std::size_t multiplicity = 0;
        while (auto quotient = symbolic::divideByLinearFactor(polynomial, root)) {
            polynomial = std::move(*quotient);
            ++multiplicity;
            if (polynomial.degree() == 0)
                break;
        }
        if (multiplicity == 0)
            return std::nullopt;
        roots.emplace_back(root, multiplicity);
    }

    if (polynomial.degree() == 1) {
        roots.emplace_back(-polynomial.coefficient(0) / polynomial.coefficient(1), 1);
    }
    else if (polynomial.degree() == 2) {
        const Rational a = polynomial.coefficient(2);
        const Rational b = polynomial.coefficient(1);
        const Rational c = polynomial.coefficient(0);
        const Rational discriminant = b * b - rational(4) * a * c;
        if (discriminant.isZero()) {
            roots.emplace_back(-b / (rational(2) * a), 2);
        }
        else if (discriminant > rational(0)) {
            const auto numeratorRoot = numeric::integerSqrt(discriminant.numerator());
            const auto denominatorRoot = numeric::integerSqrt(discriminant.denominator());
            if (!numeratorRoot.remainder.isZero() || !denominatorRoot.remainder.isZero())
                return std::nullopt;
            const Rational d{numeratorRoot.root, denominatorRoot.root};
            roots.emplace_back((-b - d) / (rational(2) * a), 1);
            roots.emplace_back((-b + d) / (rational(2) * a), 1);
        }
        // discriminant < 0 なら実根なし。符号表にはcritical pointを追加しない。
    }

    std::sort(roots.begin(), roots.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
    return roots;
}

struct RationalFunctionCriticalPoint final {
    Rational value;
    std::size_t numeratorMultiplicity = 0;
    std::size_t denominatorMultiplicity = 0;
};

[[nodiscard]] SolutionSet solveRationalFunctionInequality(
    const RationalFunctionPolynomial& function,
    const expression::Symbol& variable,
    RelationKind relation) {
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Real}};
    if (function.denominator.isZero())
        return SolutionSet::unresolved(variables);

    const auto denominatorRoots =
        rationalRealRootsWithoutIrrationalResidual(function.denominator);
    if (!denominatorRoots)
        return SolutionSet::unresolved(variables);

    if (function.numerator.isZero()) {
        if (!acceptsZero(relation))
            return SolutionSet::empty(variables);
        if (denominatorRoots->empty())
            return SolutionSet::universal(variables);

        std::vector<RealRangePiece> pieces;
        pieces.reserve(denominatorRoots->size() + 1);
        for (std::size_t i = 0; i <= denominatorRoots->size(); ++i) {
            RealRangePiece piece;
            if (i != 0)
                piece.lower = Expr{Number{(*denominatorRoots)[i - 1].first}};
            if (i != denominatorRoots->size())
                piece.upper = Expr{Number{(*denominatorRoots)[i].first}};
            pieces.push_back(std::move(piece));
        }
        std::vector<SolutionBranch> branches;
        branches.reserve(pieces.size());
        for (RealRangePiece& piece : pieces)
            branches.push_back(realRangeBranch(variable, std::move(piece)));
        return SolutionSet::finite(variables, std::move(branches));
    }

    const auto numeratorRoots =
        rationalRealRootsWithoutIrrationalResidual(function.numerator);
    if (!numeratorRoots)
        return SolutionSet::unresolved(variables);

    std::vector<RationalFunctionCriticalPoint> critical;
    auto addCritical = [&](const Rational& value, std::size_t multiplicity, bool denominator) {
        auto iterator = std::find_if(critical.begin(), critical.end(), [&](const auto& item) {
            return item.value == value;
        });
        if (iterator == critical.end()) {
            RationalFunctionCriticalPoint item;
            item.value = value;
            if (denominator)
                item.denominatorMultiplicity = multiplicity;
            else
                item.numeratorMultiplicity = multiplicity;
            critical.push_back(std::move(item));
            return;
        }
        if (denominator)
            iterator->denominatorMultiplicity += multiplicity;
        else
            iterator->numeratorMultiplicity += multiplicity;
    };
    for (const auto& [root, multiplicity] : *numeratorRoots)
        addCritical(root, multiplicity, false);
    for (const auto& [root, multiplicity] : *denominatorRoots)
        addCritical(root, multiplicity, true);
    std::sort(critical.begin(), critical.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.value < rhs.value;
    });

    const int rightSign =
        rationalSign(function.numerator.coefficient(function.numerator.degree()))
        * rationalSign(function.denominator.coefficient(function.denominator.degree()));
    std::vector<int> intervalSigns(critical.size() + 1, rightSign);
    for (std::size_t i = critical.size(); i > 0; --i) {
        intervalSigns[i - 1] = intervalSigns[i];
        const std::size_t crossingMultiplicity =
            critical[i - 1].numeratorMultiplicity
            + critical[i - 1].denominatorMultiplicity;
        if ((crossingMultiplicity & 1U) != 0U)
            intervalSigns[i - 1] = -intervalSigns[i - 1];
    }

    std::vector<RealRangePiece> pieces;
    for (std::size_t interval = 0; interval <= critical.size(); ++interval) {
        if (!acceptsSign(relation, intervalSigns[interval]))
            continue;
        RealRangePiece piece;
        if (interval != 0)
            piece.lower = Expr{Number{critical[interval - 1].value}};
        if (interval != critical.size())
            piece.upper = Expr{Number{critical[interval].value}};
        pieces.push_back(std::move(piece));
    }

    if (acceptsZero(relation)) {
        for (const auto& point : critical) {
            if (point.numeratorMultiplicity != 0 && point.denominatorMultiplicity == 0) {
                const Expr endpoint{Number{point.value}};
                pieces.push_back(RealRangePiece{endpoint, true, endpoint, true});
            }
        }
    }

    auto rank = [&](const std::optional<Expr>& endpoint, bool upper) -> std::size_t {
        if (!endpoint)
            return upper ? 2 * critical.size() + 1 : 0;
        for (std::size_t i = 0; i < critical.size(); ++i)
            if (*endpoint == Expr{Number{critical[i].value}})
                return 2 * i + 1;
        return 2 * critical.size() + 2;
    };
    std::stable_sort(pieces.begin(), pieces.end(), [&](const RealRangePiece& lhs, const RealRangePiece& rhs) {
        const std::size_t l = rank(lhs.lower, false);
        const std::size_t r = rank(rhs.lower, false);
        if (l != r)
            return l < r;
        return lhs.lowerInclusive && !rhs.lowerInclusive;
    });
    mergeRangePieces(pieces);

    if (pieces.empty())
        return SolutionSet::empty(variables);
    if (pieces.size() == 1 && !pieces.front().lower && !pieces.front().upper)
        return SolutionSet::universal(variables);

    std::vector<SolutionBranch> branches;
    branches.reserve(pieces.size());
    for (RealRangePiece& piece : pieces) {
        if (piece.lower && piece.upper && *piece.lower == *piece.upper
            && piece.lowerInclusive && piece.upperInclusive) {
            branches.push_back(branch(variable, *piece.lower));
        }
        else {
            branches.push_back(realRangeBranch(variable, std::move(piece)));
        }
    }
    return SolutionSet::finite(variables, std::move(branches));
}

[[nodiscard]] SolutionSet solveRationalPolynomialInequality(
    const symbolic::RationalPolynomial& polynomial,
    const expression::Symbol& variable,
    RelationKind relation,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Real}};
    if (polynomial.isZero())
        return acceptsZero(relation)
            ? SolutionSet::universal(variables)
            : SolutionSet::empty(variables);

    if (polynomial.degree() == 0) {
        const int sign = polynomial.coefficient(0) < rational(0) ? -1 : 1;
        return acceptsSign(relation, sign)
            ? SolutionSet::universal(variables)
            : SolutionSet::empty(variables);
    }

    const int leadingSign = polynomial.coefficient(polynomial.degree()) < rational(0) ? -1 : 1;
    std::vector<OrderedRealRoot> roots;

    if (polynomial.degree() == 1) {
        const Rational root = -polynomial.coefficient(0) / polynomial.coefficient(1);
        roots.push_back(OrderedRealRoot{Expr{Number{root}}, 1});
        return solutionFromSignChart(variable, std::move(roots), leadingSign, relation);
    }

    if (polynomial.degree() == 2) {
        const Rational a = polynomial.coefficient(2);
        const Rational b = polynomial.coefficient(1);
        const Rational c = polynomial.coefficient(0);
        const Rational discriminant = b * b - rational(4) * a * c;
        if (discriminant < rational(0))
            return acceptsSign(relation, leadingSign)
                ? SolutionSet::universal(variables)
                : SolutionSet::empty(variables);
        if (discriminant.isZero()) {
            roots.push_back(OrderedRealRoot{
                Expr{Number{(-b) / (rational(2) * a)}}, 2});
            return solutionFromSignChart(variable, std::move(roots), leadingSign, relation);
        }

        Expr minusRoot = quadraticRoot(
            a, b, discriminant, false, builtins, mathematics, angles);
        Expr plusRoot = quadraticRoot(
            a, b, discriminant, true, builtins, mathematics, angles);
        if (a < rational(0))
            std::swap(minusRoot, plusRoot);
        roots.push_back(OrderedRealRoot{std::move(minusRoot), 1});
        roots.push_back(OrderedRealRoot{std::move(plusRoot), 1});
        return solutionFromSignChart(variable, std::move(roots), leadingSign, relation);
    }

    // 高次はRational Root Theoremで線形因子をdeflateし、残りが高々二次なら
    // irrational quadratic rootsもexact radicalとして符号表へ入れる。一般の不可約三次以上
    // が残った場合は、実根なしと推測せずUnresolvedにする。
    const auto realRoots = realRootsAfterRationalDeflation(
        polynomial, builtins, mathematics, angles);
    if (!realRoots)
        return SolutionSet::unresolved(variables);
    return solutionFromSignChart(variable, *realRoots, leadingSign, relation);
}

[[nodiscard]] Expr equationZeroForm(
    const Expr& equation,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr zeroForm = equation;
    if (relationKindOf(equation, builtins)
        && equation.asCall().arguments.size() == 2) {
        zeroForm = Expr::call(
            builtins.symbol(BuiltinId::Subtract),
            {equation.asCall().arguments[0], equation.asCall().arguments[1]});
    }
    return simplify(std::move(zeroForm), builtins, mathematics, angles);
}

[[nodiscard]] bool sameVariable(
    const expression::Symbol& candidate,
    std::span<const expression::Symbol> variables) {
    return std::find(variables.begin(), variables.end(), candidate) != variables.end();
}

[[nodiscard]] bool containsAnySolverVariable(
    const Expr& expression,
    std::span<const expression::Symbol> variables) {
    for (const expression::Symbol& variable : variables)
        if (symbolic::containsSymbol(expression, variable))
            return true;
    return false;
}

void mergeAssumptions(
    mathematics::AssumptionSet& target,
    const mathematics::AssumptionSet& source) {
    for (const mathematics::Predicate& predicate : source.predicates())
        target.add(predicate);
}

[[nodiscard]] Expr multiplyFactors(
    std::vector<Expr> factors,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (factors.empty())
        return Expr{Number{BigInt{1}}};
    if (factors.size() == 1)
        return factors.front();
    return simplify(
        Expr::call(builtins.symbol(BuiltinId::Multiply), std::move(factors)),
        builtins, mathematics, angles);
}

struct SymbolicLinearEquation final {
    std::vector<Expr> coefficients;
    Expr rhs;
    mathematics::AssumptionSet domainConditions;
};

[[nodiscard]] std::optional<SymbolicLinearEquation> extractSymbolicLinearEquation(
    const Expr& equation,
    std::span<const expression::Symbol> variables,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr zeroForm = symbolic::expandExpression(
        equationZeroForm(equation, builtins, mathematics, angles),
        builtins, mathematics, angles);

    std::vector<Expr> terms;
    if (isHead(zeroForm, builtins, BuiltinId::Add))
        terms.assign(zeroForm.asCall().arguments.begin(), zeroForm.asCall().arguments.end());
    else
        terms.push_back(zeroForm);

    std::vector<Expr> coefficients(variables.size(), zeroExpr());
    Expr constant = zeroExpr();
    mathematics::AssumptionSet domainConditions;

    auto addTo = [&](Expr& destination, Expr value) {
        destination = simplify(
            Expr::call(builtins.symbol(BuiltinId::Add), {destination, std::move(value)}),
            builtins, mathematics, angles);
    };

    for (const Expr& term : terms) {
        if (!containsAnySolverVariable(term, variables)) {
            const auto domain = symbolic::scalarExpressionDomainConditions(term, builtins, mathematics);
            if (!domain)
                return std::nullopt;
            mergeAssumptions(domainConditions, *domain);
            addTo(constant, term);
            continue;
        }

        std::optional<std::size_t> variableIndex;
        std::vector<Expr> coefficientFactors;
        if (term.isSymbol()) {
            const auto iterator = std::find(variables.begin(), variables.end(), term.asSymbol());
            if (iterator == variables.end())
                return std::nullopt;
            variableIndex = static_cast<std::size_t>(iterator - variables.begin());
        }
        else if (isHead(term, builtins, BuiltinId::Multiply)) {
            for (const Expr& factor : term.asCall().arguments) {
                if (factor.isSymbol()) {
                    const auto iterator = std::find(variables.begin(), variables.end(), factor.asSymbol());
                    if (iterator != variables.end()) {
                        if (variableIndex)
                            return std::nullopt; // xy や x^2 は線形ではない。
                        variableIndex = static_cast<std::size_t>(iterator - variables.begin());
                        continue;
                    }
                }
                if (containsAnySolverVariable(factor, variables))
                    return std::nullopt;
                coefficientFactors.push_back(factor);
            }
            if (!variableIndex)
                return std::nullopt;
        }
        else
            return std::nullopt;

        Expr coefficient = multiplyFactors(
            std::move(coefficientFactors), builtins, mathematics, angles);
        const auto domain = symbolic::scalarExpressionDomainConditions(
            coefficient, builtins, mathematics);
        if (!domain)
            return std::nullopt;
        mergeAssumptions(domainConditions, *domain);
        addTo(coefficients[*variableIndex], std::move(coefficient));
    }

    Expr rhs = simplify(
        Expr::call(builtins.symbol(BuiltinId::Negate), {constant}),
        builtins, mathematics, angles);
    return SymbolicLinearEquation{
        std::move(coefficients), std::move(rhs), std::move(domainConditions)};
}

[[nodiscard]] Expr symbolicDeterminant(
    const std::vector<std::vector<Expr>>& matrix,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::size_t n = matrix.size();
    if (n == 0)
        return Expr{Number{BigInt{1}}};
    if (n == 1)
        return matrix.front().front();
    if (n == 2) {
        Expr first = Expr::call(
            builtins.symbol(BuiltinId::Multiply), {matrix[0][0], matrix[1][1]});
        Expr second = Expr::call(
            builtins.symbol(BuiltinId::Multiply), {matrix[0][1], matrix[1][0]});
        return simplify(
            Expr::call(builtins.symbol(BuiltinId::Subtract), {std::move(first), std::move(second)}),
            builtins, mathematics, angles);
    }

    Expr determinant = zeroExpr();
    for (std::size_t column = 0; column < n; ++column) {
        std::vector<std::vector<Expr>> minor;
        minor.reserve(n - 1);
        for (std::size_t row = 1; row < n; ++row) {
            std::vector<Expr> minorRow;
            minorRow.reserve(n - 1);
            for (std::size_t c = 0; c < n; ++c)
                if (c != column)
                    minorRow.push_back(matrix[row][c]);
            minor.push_back(std::move(minorRow));
        }
        Expr term = Expr::call(
            builtins.symbol(BuiltinId::Multiply),
            {matrix[0][column], symbolicDeterminant(minor, builtins, mathematics, angles)});
        if ((column & 1U) != 0)
            term = Expr::call(builtins.symbol(BuiltinId::Negate), {std::move(term)});
        determinant = Expr::call(
            builtins.symbol(BuiltinId::Add), {std::move(determinant), std::move(term)});
    }
    return simplify(std::move(determinant), builtins, mathematics, angles);
}

[[nodiscard]] SolutionSet solveSymbolicLinearSystem(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> variableSymbols,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<SolverVariable> variables;
    variables.reserve(variableSymbols.size());
    for (const expression::Symbol& variable : variableSymbols)
        variables.push_back(SolverVariable{variable, mathematics::NumericDomain::Complex});

    // Cramerの公式を基準実装にする。式爆発を避けるため、symbolic coefficient版は
    // まず4元まで。Rational係数系は下のGauss-Jordanが任意元を扱う。
    if (variableSymbols.empty() || equations.size() != variableSymbols.size()
        || variableSymbols.size() > 4)
        return SolutionSet::unresolved(variables);

    std::vector<std::vector<Expr>> matrix;
    std::vector<Expr> rhs;
    mathematics::AssumptionSet domainConditions;
    matrix.reserve(equations.size());
    rhs.reserve(equations.size());
    for (const Expr& equation : equations) {
        auto linear = extractSymbolicLinearEquation(
            equation, variableSymbols, builtins, mathematics, angles);
        if (!linear)
            return SolutionSet::unresolved(variables);
        mergeAssumptions(domainConditions, linear->domainConditions);
        matrix.push_back(std::move(linear->coefficients));
        rhs.push_back(std::move(linear->rhs));
    }

    const Expr determinant = symbolicDeterminant(matrix, builtins, mathematics, angles);
    const auto determinantZero = proveZero(determinant, builtins, mathematics);
    const auto determinantNonZero = proveNonZero(determinant, builtins, mathematics);
    if (determinantZero == mathematics::TruthValue::True)
        return SolutionSet::unresolved(variables).withAdditionalConditions(domainConditions);

    SolutionBranch unique;
    unique.bindings.reserve(variableSymbols.size());
    for (std::size_t column = 0; column < variableSymbols.size(); ++column) {
        auto replaced = matrix;
        for (std::size_t row = 0; row < replaced.size(); ++row)
            replaced[row][column] = rhs[row];
        Expr numerator = symbolicDeterminant(replaced, builtins, mathematics, angles);
        Expr value = simplify(
            Expr::call(builtins.symbol(BuiltinId::Divide), {std::move(numerator), determinant}),
            builtins, mathematics, angles);
        unique.bindings.push_back(SolutionBinding{variableSymbols[column], std::move(value)});
    }

    if (determinantNonZero == mathematics::TruthValue::True)
        return SolutionSet::finite(variables, {std::move(unique)})
            .withAdditionalConditions(domainConditions);

    SolutionSet result = SolutionSet::conditional(
        variables,
        {
            finiteCase(conditions({notEqualZero(determinant)}), {std::move(unique)}),
            SolutionCase{conditions({equalZero(determinant)}), SolutionSetKind::Unresolved, {}}
        });
    return result.withAdditionalConditions(domainConditions);
}


[[nodiscard]] std::optional<SolutionSet> solveDirectAlgebraicBinding(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    if (!isHead(relation, builtins, BuiltinId::Equal)
        || relation.asCall().arguments.size() != 2)
        return std::nullopt;

    const Expr& lhs = relation.asCall().arguments[0];
    const Expr& rhs = relation.asCall().arguments[1];
    const auto isVariable = [&](const Expr& value) {
        return value.isSymbol() && value.asSymbol() == variable;
    };

    const Expr* candidate = nullptr;
    if (isVariable(lhs) && !symbolic::containsSymbol(rhs, variable))
        candidate = &rhs;
    else if (isVariable(rhs) && !symbolic::containsSymbol(lhs, variable))
        candidate = &lhs;
    if (!candidate)
        return std::nullopt;

    // 7-6ではexact algebraic valueだけをdirect bindingする。Infinityや一般symbolic式まで
    // 「Complexの解」と推測せず，既存solverへfallbackする。
    if (!symbolic::exactAlgebraicValue(*candidate, builtins, mathematics))
        return std::nullopt;

    return SolutionSet::finite(
        {{variable, mathematics::NumericDomain::Complex}},
        {branch(variable, *candidate)});
}

} // namespace

std::optional<SolutionSet> solveDirectAlgebraicBindingRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    return solveDirectAlgebraicBinding(relation, variable, builtins, mathematics);
}

SolutionSet solveUnivariatePolynomialRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (const auto direct = solveDirectAlgebraicBinding(
            relation, variable, builtins, mathematics))
        return *direct;

    const auto relationKind = relationKindOf(relation, builtins);
    if (!relationKind || *relationKind == RelationKind::Equal)
        return solvePolynomialEquation(relation, variable, builtins, mathematics, angles);

    if (*relationKind == RelationKind::NotEqual) {
        // != は順序を必要としないのでComplex上でも意味を持つ。多項式のzero setを
        // exactに列挙できる場合、その有限集合を除いた補集合として保持する。
        const auto& args = relation.asCall().arguments;
        Expr equality = Expr::call(
            builtins.symbol(BuiltinId::Equal), {args[0], args[1]});
        SolutionSet zeros = solvePolynomialEquation(
            equality, variable, builtins, mathematics, angles);
        const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Complex}};
        if (zeros.kind() == SolutionSetKind::Empty)
            return SolutionSet::universal(variables);
        if (zeros.kind() == SolutionSetKind::Universal)
            return SolutionSet::empty(variables);
        if (zeros.kind() != SolutionSetKind::Finite)
            return SolutionSet::unresolved(variables);

        mathematics::AssumptionSet exclusions;
        for (const SolutionBranch& zero : zeros.branches()) {
            if (zero.bindings.size() != 1 || zero.bindings.front().variable != variable)
                return SolutionSet::unresolved(variables);
            exclusions.add(mathematics::relation(
                RelationKind::NotEqual, Expr{variable}, zero.bindings.front().value));
        }
        return SolutionSet::finite(
            variables,
            {SolutionBranch{
                {}, std::move(exclusions), std::nullopt,
                {SolverVariable{variable, mathematics::NumericDomain::Complex}}, std::nullopt}});
    }

    const Expr zeroForm = equationZeroForm(relation, builtins, mathematics, angles);
    if (const auto polynomial = symbolic::toRationalPolynomial(zeroForm, variable, builtins)) {
        return solveRationalPolynomialInequality(
            *polynomial, variable, *relationKind, builtins, mathematics, angles);
    }

    if (const auto rationalFunction =
        toRationalFunctionPolynomial(zeroForm, variable, builtins)) {
        return solveRationalFunctionInequality(
            *rationalFunction, variable, *relationKind);
    }

    return SolutionSet::unresolved(
        {{variable, mathematics::NumericDomain::Real}});
}

[[nodiscard]] symbolic::RationalPolynomial derivativePolynomial(
    const symbolic::RationalPolynomial& polynomial) {
    if (polynomial.degree() == 0)
        return symbolic::RationalPolynomial{};
    std::vector<Rational> coefficients(polynomial.degree(), rational(0));
    for (std::size_t exponent = 1; exponent <= polynomial.degree(); ++exponent)
        coefficients[exponent - 1] = polynomial.coefficient(exponent)
            * Rational{BigInt::fromUnsigned(exponent)};
    return symbolic::RationalPolynomial{std::move(coefficients)};
}

[[nodiscard]] symbolic::RationalPolynomial polynomialRemainder(
    const symbolic::RationalPolynomial& numerator,
    const symbolic::RationalPolynomial& denominator) {
    if (denominator.isZero())
        throw std::invalid_argument("Polynomial remainder requires nonzero denominator");
    std::vector<Rational> remainder = numerator.coefficients();
    const auto trim = [](std::vector<Rational>& coefficients) {
        while (coefficients.size() > 1 && coefficients.back().isZero())
            coefficients.pop_back();
        if (coefficients.empty())
            coefficients.push_back(rational(0));
    };
    trim(remainder);
    const Rational leading = denominator.coefficient(denominator.degree());
    while (!(remainder.size() == 1 && remainder.front().isZero())
        && remainder.size() - 1 >= denominator.degree()) {
        const std::size_t shift = remainder.size() - 1 - denominator.degree();
        const Rational factor = remainder.back() / leading;
        for (std::size_t i = 0; i <= denominator.degree(); ++i)
            remainder[i + shift] -= factor * denominator.coefficient(i);
        trim(remainder);
    }
    return symbolic::RationalPolynomial{std::move(remainder)};
}

[[nodiscard]] bool hasRepeatedPolynomialFactor(
    const symbolic::RationalPolynomial& polynomial) {
    if (polynomial.degree() < 2)
        return false;
    symbolic::RationalPolynomial lhs = polynomial;
    symbolic::RationalPolynomial rhs = derivativePolynomial(polynomial);
    while (!rhs.isZero()) {
        symbolic::RationalPolynomial remainder = polynomialRemainder(lhs, rhs);
        lhs = std::move(rhs);
        rhs = std::move(remainder);
    }
    return lhs.degree() > 0;
}

std::optional<SolutionSet> solveFactoredRealAlgebraicPolynomialEquation(
    const Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto relation = relationKindOf(equation, builtins);
    if (relation && *relation != RelationKind::Equal)
        return std::nullopt;

    const Expr zeroForm = equationZeroForm(equation, builtins, mathematics, angles);
    const auto polynomial = symbolic::toRationalPolynomial(zeroForm, variable, builtins);
    if (!polynomial || polynomial->degree() <= 2 || polynomial->isZero())
        return std::nullopt;

    if (const auto rationalRoots = rationalRealRootsWithoutIrrationalResidual(*polynomial)) {
        const std::vector<SolverVariable> variables{{
            variable, mathematics::NumericDomain::Real}};
        std::vector<SolutionBranch> branches;
        branches.reserve(rationalRoots->size());
        for (const auto& [root, multiplicity] : *rationalRoots)
            branches.push_back(branch(variable, Expr{Number{root}}, multiplicity));
        return branches.empty()
            ? SolutionSet::empty(variables)
            : SolutionSet::finite(variables, std::move(branches));
    }

    Expr factored = symbolic::factorExpression(zeroForm, builtins, mathematics, angles);
    if (!isHead(factored, builtins, BuiltinId::Multiply))
        return std::nullopt;

    std::vector<symbolic::RealAlgebraicNumber> roots;
    std::size_t nonconstantFactors = 0;
    for (const Expr& factor : factored.asCall().arguments) {
        const auto factorPolynomial = symbolic::toRationalPolynomial(factor, variable, builtins);
        if (!factorPolynomial) {
            if (!symbolic::containsSymbol(factor, variable))
                continue;
            return std::nullopt;
        }
        if (factorPolynomial->degree() == 0)
            continue;
        ++nonconstantFactors;
        const auto factorRoots = symbolic::RealAlgebraicNumber::isolateAll(
            factorPolynomial->coefficients());
        if (!factorRoots)
            return std::nullopt;
        roots.insert(roots.end(), factorRoots->begin(), factorRoots->end());
    }
    if (nonconstantFactors < 2)
        return std::nullopt;

    std::sort(roots.begin(), roots.end(),
        [](const symbolic::RealAlgebraicNumber& lhs,
           const symbolic::RealAlgebraicNumber& rhs) {
            const auto left = symbolic::AlgebraicNumber::fromRealRoot(lhs);
            const auto right = symbolic::AlgebraicNumber::fromRealRoot(rhs);
            const auto order = left.exactRealCompare(right);
            if (order)
                return *order == symbolic::AlgebraicOrder::Less;
            return lhs.isolatingInterval().upper <= rhs.isolatingInterval().lower;
        });
    roots.erase(std::unique(roots.begin(), roots.end(),
        [](const symbolic::RealAlgebraicNumber& lhs,
           const symbolic::RealAlgebraicNumber& rhs) {
            const auto left = symbolic::AlgebraicNumber::fromRealRoot(lhs);
            const auto right = symbolic::AlgebraicNumber::fromRealRoot(rhs);
            return left.exactEquals(right).value_or(false);
        }), roots.end());

    const std::vector<SolverVariable> variables{{
        variable, mathematics::NumericDomain::Real}};
    std::vector<SolutionBranch> branches;
    branches.reserve(roots.size());
    for (const auto& root : roots) {
        const auto canonical = symbolic::RealAlgebraicNumber::create(
            root.polynomial(), root.rootIndex());
        branches.push_back(branch(variable,
            symbolic::makeCanonicalRootExpression(
                canonical ? *canonical : root, builtins)));
    }
    return branches.empty()
        ? SolutionSet::empty(variables)
        : SolutionSet::finite(variables, std::move(branches));
}

std::optional<SolutionSet> solveRealAlgebraicPolynomialEquation(
    const Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto relation = relationKindOf(equation, builtins);
    if (relation && *relation != RelationKind::Equal)
        return std::nullopt;

    const Expr zeroForm = equationZeroForm(equation, builtins, mathematics, angles);
    const auto polynomial = symbolic::toRationalPolynomial(zeroForm, variable, builtins);
    if (!polynomial || polynomial->degree() <= 2)
        return std::nullopt;
    if (polynomial->isZero())
        return SolutionSet::universal({{variable, mathematics::NumericDomain::Real}});

    // 有理根だけで実根集合が閉じる場合は，Rootへ持ち上げずそのまま返す。
    // complexな残余因子は実解へ寄与しないため，x^4-1のような基本形のcanonical表示も保てる。
    if (const auto rationalRoots = rationalRealRootsWithoutIrrationalResidual(*polynomial)) {
        const std::vector<SolverVariable> variables{{
            variable, mathematics::NumericDomain::Real}};
        std::vector<SolutionBranch> branches;
        branches.reserve(rationalRoots->size());
        for (const auto& [root, multiplicity] : *rationalRoots)
            branches.push_back(branch(variable, Expr{Number{root}}, multiplicity));
        return branches.empty()
            ? SolutionSet::empty(variables)
            : SolutionSet::finite(variables, std::move(branches));
    }

    // x^4+b*x^2+cでy=x^2が異なる正の有理根を2個持つ場合は，
    // 4次全体を分離せず二つのx^2-y因子を直接分離する。
    if (polynomial->degree() == 4
        && polynomial->coefficient(1).isZero()
        && polynomial->coefficient(3).isZero()
        && polynomial->coefficient(4) == rational(1)) {
        const Rational b = polynomial->coefficient(2);
        const Rational c = polynomial->coefficient(0);
        const Rational discriminant = b * b - rational(4) * c;
        const auto numeratorRoot = numeric::integerSqrt(discriminant.numerator());
        const auto denominatorRoot = numeric::integerSqrt(discriminant.denominator());
        if (discriminant >= rational(0)
            && numeratorRoot.remainder.isZero()
            && denominatorRoot.remainder.isZero()) {
            const Rational sqrtDiscriminant{
                numeratorRoot.root, denominatorRoot.root};
            const Rational y1 = (-b - sqrtDiscriminant) / rational(2);
            const Rational y2 = (-b + sqrtDiscriminant) / rational(2);
            const auto isRationalSquare = [](const Rational& value) {
                if (value < rational(0)) return false;
                const auto n = numeric::integerSqrt(value.numerator());
                const auto d = numeric::integerSqrt(value.denominator());
                return n.remainder.isZero() && d.remainder.isZero();
            };
            if (y1 > rational(0) && y2 > rational(0)
                && y1 != y2 && !isRationalSquare(y1) && !isRationalSquare(y2)) {
                std::vector<symbolic::RealAlgebraicNumber> roots;
                for (const Rational& y : {y1, y2}) {
                    const std::array<Rational, 3> factor{-y, rational(0), rational(1)};
                    const auto factorRoots = symbolic::RealAlgebraicNumber::isolateAll(factor);
                    if (!factorRoots || factorRoots->size() != 2) {
                        roots.clear();
                        break;
                    }
                    roots.insert(roots.end(), factorRoots->begin(), factorRoots->end());
                }
                if (roots.size() == 4) {
                    std::sort(roots.begin(), roots.end(),
                        [](const auto& lhs, const auto& rhs) {
                            const auto left = symbolic::AlgebraicNumber::fromRealRoot(lhs);
                            const auto right = symbolic::AlgebraicNumber::fromRealRoot(rhs);
                            return left.exactRealCompare(right)
                                == symbolic::AlgebraicOrder::Less;
                        });
                    const std::vector<SolverVariable> variables{{
                        variable, mathematics::NumericDomain::Real}};
                    std::vector<SolutionBranch> branches;
                    for (const auto& root : roots)
                        branches.push_back(branch(variable,
                            symbolic::makeCanonicalRootExpression(root, builtins)));
                    return SolutionSet::finite(variables, std::move(branches));
                }
            }
        }
    }

    // reducible polynomialは全体を高次数のまま分離してからminimal polynomialへ戻さず，
    // Q上の因子ごとに実根を分離する。小次数因子のSturm計算を共有できるため，
    // (x^2-a)(x^2-b)のような基本形で高次数root分離の固定費を払わない。
    Expr factored = symbolic::factorExpression(
        zeroForm, builtins, mathematics, angles);
    if (isHead(factored, builtins, BuiltinId::Multiply)) {
        std::vector<symbolic::RealAlgebraicNumber> factorRoots;
        std::size_t nonconstantFactors = 0;
        bool usable = true;
        for (const Expr& factor : factored.asCall().arguments) {
            const auto factorPolynomial = symbolic::toRationalPolynomial(
                factor, variable, builtins);
            if (!factorPolynomial) {
                if (!symbolic::containsSymbol(factor, variable))
                    continue;
                usable = false;
                break;
            }
            if (factorPolynomial->degree() == 0)
                continue;
            ++nonconstantFactors;
            const auto roots = symbolic::RealAlgebraicNumber::isolateAll(
                factorPolynomial->coefficients());
            if (!roots) {
                usable = false;
                break;
            }
            factorRoots.insert(factorRoots.end(), roots->begin(), roots->end());
        }
        if (usable && nonconstantFactors >= 2) {
            std::sort(factorRoots.begin(), factorRoots.end(),
                [](const symbolic::RealAlgebraicNumber& lhs,
                   const symbolic::RealAlgebraicNumber& rhs) {
                    const auto left = symbolic::AlgebraicNumber::fromRealRoot(lhs);
                    const auto right = symbolic::AlgebraicNumber::fromRealRoot(rhs);
                    const auto order = left.exactRealCompare(right);
                    if (order)
                        return *order == symbolic::AlgebraicOrder::Less;
                    return lhs.isolatingInterval().upper <= rhs.isolatingInterval().lower;
                });
            factorRoots.erase(std::unique(factorRoots.begin(), factorRoots.end(),
                [](const symbolic::RealAlgebraicNumber& lhs,
                   const symbolic::RealAlgebraicNumber& rhs) {
                    const auto left = symbolic::AlgebraicNumber::fromRealRoot(lhs);
                    const auto right = symbolic::AlgebraicNumber::fromRealRoot(rhs);
                    return left.exactEquals(right).value_or(false);
                }), factorRoots.end());
            const std::vector<SolverVariable> variables{{
                variable, mathematics::NumericDomain::Real}};
            std::vector<SolutionBranch> branches;
            branches.reserve(factorRoots.size());
            for (const auto& root : factorRoots) {
                const auto canonical = symbolic::RealAlgebraicNumber::create(
                    root.polynomial(), root.rootIndex());
                const auto& value = canonical ? *canonical : root;
                if (value.degree() == 1) {
                    const auto coefficients = value.polynomial();
                    branches.push_back(branch(variable, Expr{Number{
                        -coefficients[0] / coefficients[1]}}));
                    continue;
                }
                branches.push_back(branch(variable,
                    symbolic::makeCanonicalRootExpression(value, builtins)));
            }
            return branches.empty()
                ? SolutionSet::empty(variables)
                : SolutionSet::finite(variables, std::move(branches));
        }
    }

    const auto roots = symbolic::RealAlgebraicNumber::isolateAll(polynomial->coefficients());
    if (!roots)
        return std::nullopt;
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Real}};
    if (roots->empty())
        return SolutionSet::empty(variables);

    std::vector<SolutionBranch> branches;
    branches.reserve(roots->size());
    for (const symbolic::RealAlgebraicNumber& root : *roots) {
        const auto canonical = symbolic::RealAlgebraicNumber::create(
            root.polynomial(), root.rootIndex());
        branches.push_back(branch(variable,
            symbolic::makeCanonicalRootExpression(canonical ? *canonical : root, builtins)));
    }
    return SolutionSet::finite(variables, std::move(branches));
}

std::optional<SolutionSet> solveRepeatedRealAlgebraicPolynomialEquation(
    const Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const auto relation = relationKindOf(equation, builtins);
    if (relation && *relation != RelationKind::Equal)
        return std::nullopt;
    const Expr zeroForm = equationZeroForm(equation, builtins, mathematics, angles);
    const auto polynomial = symbolic::toRationalPolynomial(zeroForm, variable, builtins);
    if (!polynomial || polynomial->degree() <= 2 || !hasRepeatedPolynomialFactor(*polynomial))
        return std::nullopt;
    return solveRealAlgebraicPolynomialEquation(
        equation, variable, builtins, mathematics, angles);
}

SolutionSet solvePolynomialEquation(
    const Expr& equation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Complex}};
    const Expr zeroForm = equationZeroForm(equation, builtins, mathematics, angles);
    const auto polynomial = symbolic::toRationalPolynomial(zeroForm, variable, builtins);
    if (!polynomial) {
        if (const auto rationalFunction =
            toRationalFunctionPolynomial(zeroForm, variable, builtins)) {
            // f(x)/g(x) == 0 は f(x)==0 かつ g(x)!=0。分母を消したことによる
            // extraneous rootをSolveConstraintsでexactに除外する。
            const Expr numeratorExpr = symbolic::polynomialToExpandedExpr(
                rationalFunction->numerator, variable, builtins);
            const Expr numeratorEquation = Expr::call(
                builtins.symbol(BuiltinId::Equal), {numeratorExpr, zeroExpr()});
            SolutionSet result = solvePolynomialEquation(
                numeratorEquation, variable, builtins, mathematics, angles);

            const Expr denominatorExpr = symbolic::polynomialToExpandedExpr(
                rationalFunction->denominator, variable, builtins);
            SolveConstraints denominatorConstraint;
            denominatorConstraint.assumptions.add(mathematics::relation(
                RelationKind::NotEqual, denominatorExpr, zeroExpr()));
            return applySolveConstraints(
                std::move(result), denominatorConstraint, builtins, mathematics, angles);
        }

        if (const auto conversion = symbolic::toExpressionPolynomialWithConditions(
            zeroForm, variable, builtins, mathematics, angles)) {
            SolutionSet result = solveExpressionPolynomial(
                conversion->polynomial, variable, builtins, mathematics, angles);
            return result.withAdditionalConditions(conversion->domainConditions);
        }
        return SolutionSet::unresolved(variables);
    }
    if (polynomial->isZero())
        return SolutionSet::universal(variables);

    switch (polynomial->degree()) {
    case 0:
        return SolutionSet::empty(variables);
    case 1: {
        const Rational root = -polynomial->coefficient(0) / polynomial->coefficient(1);
        return SolutionSet::finite(variables, {branch(variable, Expr{Number{root}})});
    }
    case 2:
        return SolutionSet::finite(
            variables,
            solveQuadratic(*polynomial, variable, builtins, mathematics, angles));
    default:
        break;
    }

    if (const auto cube = solvePerfectCubeBinomial(
        *polynomial, variable, builtins, mathematics, angles))
        return SolutionSet::finite(variables, *cube);

    if (isBinomial(*polynomial) && polynomial->degree() <= 256)
        return SolutionSet::finite(
            variables,
            solveBinomial(*polynomial, variable, builtins, mathematics, angles));

    if (auto branches = solveByRationalDeflation(
        *polynomial, variable, builtins, mathematics, angles))
        return branches->empty()
            ? SolutionSet::empty(variables)
            : SolutionSet::finite(variables, std::move(*branches));

    if (const auto roots = symbolic::ComplexAlgebraicNumber::isolateAll(polynomial->coefficients())) {
        if (roots->empty())
            return SolutionSet::empty(variables);
        const auto canonicalRoots = symbolic::ComplexAlgebraicNumber::canonicalizeAll(*roots);
        std::vector<SolutionBranch> branches;
        branches.reserve(roots->size());
        if (canonicalRoots) {
            for (const symbolic::ComplexAlgebraicNumber& root : *canonicalRoots)
                branches.push_back(branch(
                    variable, symbolic::makeCanonicalRootExpression(root, builtins)));
        }
        else {
            for (const symbolic::ComplexAlgebraicNumber& root : *roots) {
                const auto canonical = symbolic::ComplexAlgebraicNumber::create(
                    root.polynomial(), root.rootIndex());
                branches.push_back(branch(variable,
                    symbolic::makeCanonicalRootExpression(canonical ? *canonical : root, builtins)));
            }
        }
        return SolutionSet::finite(variables, std::move(branches));
    }
    return SolutionSet::unresolved(variables);
}


[[nodiscard]] SolutionSet applyAmbientDomainToSystemSolution(
    SolutionSet solution,
    mathematics::NumericDomain ambientDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (ambientDomain == mathematics::NumericDomain::Complex)
        return solution;
    SolveConstraints constraints;
    constraints.domain = ambientDomain;
    return applySolveConstraints(
        std::move(solution), constraints, builtins, mathematics, angles);
}

[[nodiscard]] SolutionSet solvePolynomialSubsystem(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> variables,
    mathematics::NumericDomain ambientDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    std::vector<SolverVariable> solverVariables;
    solverVariables.reserve(variables.size());
    for (const expression::Symbol& variable : variables)
        solverVariables.push_back(SolverVariable{variable, ambientDomain});

    if (variables.empty()) {
        for (const Expr& equation : equations) {
            const Expr zeroForm = equationZeroForm(equation, builtins, mathematics, angles);
            if (proveZero(zeroForm, builtins, mathematics) == mathematics::TruthValue::True)
                continue;
            if (proveNonZero(zeroForm, builtins, mathematics) == mathematics::TruthValue::True)
                return SolutionSet::empty({});
            return SolutionSet::unresolved({});
        }
        return SolutionSet::universal({});
    }

    if (variables.size() == 1) {
        if (equations.empty())
            return SolutionSet::universal(std::move(solverVariables));

        SolutionSet result = solveUnivariatePolynomialRelation(
            equations.front(), variables.front(), builtins, mathematics, angles);
        for (std::size_t i = 1; i < equations.size(); ++i) {
            const SolveConstraints relationConstraint = parseSolveConstraints(
                equations[i], variables, builtins, mathematics, angles);
            result = applySolveConstraints(
                std::move(result), relationConstraint, builtins, mathematics, angles);
        }
        return applyAmbientDomainToSystemSolution(
            std::move(result), ambientDomain, builtins, mathematics, angles);
    }

    SolutionSet result = [&] {
        if (auto polynomial = solvePolynomialSystem(
                equations, variables, ambientDomain, builtins, mathematics, angles))
            return std::move(*polynomial);
        return solveLinearPolynomialSystem(
            equations, variables, builtins, mathematics, angles);
    }();
    return applyAmbientDomainToSystemSolution(
        std::move(result), ambientDomain, builtins, mathematics, angles);
}

[[nodiscard]] std::optional<SolutionSet> solveByEliminantSpecialization(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> orderedSymbols,
    std::span<const expression::Symbol> outputSymbols,
    mathematics::NumericDomain ambientDomain,
    const expression::Symbol& lastVariable,
    const SolutionSet& lastSolutions,
    bool lastSolutionsCertifiedReal,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    constexpr std::size_t maximumSpecializationRoots = 64;
    constexpr std::size_t maximumSpecializationVariables = 4;
    if (orderedSymbols.size() > maximumSpecializationVariables
        || lastSolutions.branches().size() > maximumSpecializationRoots)
        return std::nullopt;

    std::vector<expression::Symbol> remainingSymbols{
        orderedSymbols.begin(), orderedSymbols.end() - 1};
    std::vector<SolverVariable> resultVariables;
    resultVariables.reserve(outputSymbols.size());
    for (const expression::Symbol& variable : outputSymbols)
        resultVariables.push_back(SolverVariable{variable, ambientDomain});

    std::vector<SolutionBranch> resultBranches;
    for (const SolutionBranch& rootBranch : lastSolutions.branches()) {
        if (!rootBranch.unconditional() || !rootBranch.freeVariables.empty()
            || rootBranch.bindings.size() != 1
            || !(rootBranch.bindings.front().variable == lastVariable))
            return std::nullopt;

        const Expr& root = rootBranch.bindings.front().value;
        std::vector<Expr> specializedEquations;
        specializedEquations.reserve(equations.size());
        for (const Expr& equation : equations) {
            Expr specialized = symbolic::substituteSymbol(equation, lastVariable, root);
            specializedEquations.push_back(
                simplify(std::move(specialized), builtins, mathematics, angles));
        }

        SolutionSet remaining = solvePolynomialSubsystem(
            specializedEquations, remainingSymbols, ambientDomain,
            builtins, mathematics, angles);
        if (remaining.kind() == SolutionSetKind::Unresolved
            || remaining.kind() == SolutionSetKind::Conditional)
            return std::nullopt;
        if (remaining.kind() == SolutionSetKind::Empty)
            continue;

        auto appendBranch = [&](const SolutionBranch* subBranch) {
            SolutionBranch combined;
            if (subBranch) {
                combined.bindings.assign(
                    subBranch->bindings.begin(), subBranch->bindings.end());
                combined.freeVariables.assign(
                    subBranch->freeVariables.begin(), subBranch->freeVariables.end());
                combined.conditions = subBranch->conditions;
                for (const mathematics::Predicate& predicate : remaining.conditions().predicates())
                    combined.conditions.add(predicate);
            }
            else {
                for (const expression::Symbol& variable : remainingSymbols)
                    combined.freeVariables.push_back(SolverVariable{variable, ambientDomain});
                combined.conditions = remaining.conditions();
            }
            combined.bindings.push_back(SolutionBinding{lastVariable, root});

            std::vector<SolutionBinding> orderedBindings;
            orderedBindings.reserve(combined.bindings.size());
            for (const expression::Symbol& variable : outputSymbols) {
                const auto found = std::find_if(
                    combined.bindings.begin(), combined.bindings.end(),
                    [&](const SolutionBinding& binding) { return binding.variable == variable; });
                if (found != combined.bindings.end())
                    orderedBindings.push_back(*found);
            }
            combined.bindings = std::move(orderedBindings);

            if (ambientDomain == mathematics::NumericDomain::Real
                && lastSolutionsCertifiedReal
                && (!subBranch || (subBranch->bindingsCertifiedDomain
                    && mathematics::isSubdomainOf(
                        *subBranch->bindingsCertifiedDomain,
                        mathematics::NumericDomain::Real))))
                combined.bindingsCertifiedDomain = mathematics::NumericDomain::Real;

            if (std::find(resultBranches.begin(), resultBranches.end(), combined)
                == resultBranches.end())
                resultBranches.push_back(std::move(combined));
        };

        if (remaining.kind() == SolutionSetKind::Universal)
            appendBranch(nullptr);
        else {
            for (const SolutionBranch& branch : remaining.branches())
                appendBranch(&branch);
        }
    }

    if (resultBranches.empty())
        return SolutionSet::empty(std::move(resultVariables));
    return SolutionSet::finite(std::move(resultVariables), std::move(resultBranches));
}


[[nodiscard]] std::optional<SolutionSet> solveSinglePolynomialProjection(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> variables,
    mathematics::NumericDomain ambientDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if ((ambientDomain != mathematics::NumericDomain::Complex
            && ambientDomain != mathematics::NumericDomain::Real)
        || equations.size() != 1 || variables.size() < 2
        || relationKindOf(equations.front(), builtins) != RelationKind::Equal)
        return std::nullopt;

    const Expr zeroForm = equationZeroForm(
        equations.front(), builtins, mathematics, angles);
    const auto multivariate = symbolic::toMultivariateRationalPolynomial(
        zeroForm, builtins, symbolic::PolynomialConversionOptions{64, 256});
    if (!multivariate || multivariate->totalDegree() <= 1)
        return std::nullopt;
    for (const expression::Symbol& symbol : multivariate->variables())
        if (!sameVariable(symbol, variables))
            return std::nullopt;

    const expression::Symbol* solvedVariable = nullptr;
    for (auto iterator = variables.rbegin(); iterator != variables.rend(); ++iterator) {
        const std::size_t degree = multivariate->degree(*iterator);
        if (degree >= 1 && degree <= 2) {
            solvedVariable = &*iterator;
            break;
        }
    }
    if (!solvedVariable)
        return std::nullopt;

    const auto conversion = symbolic::toExpressionPolynomialWithConditions(
        zeroForm, *solvedVariable, builtins, mathematics, angles,
        symbolic::PolynomialConversionOptions{2, 256});
    if (!conversion || conversion->polynomial.degree() > 2)
        return std::nullopt;

    // Realの正次元二次射影では，複素二次公式をそのまま返すと
    // sqrt[discriminant]の実在条件が自由parameterから落ちる。
    // 最高次係数が常に非零と証明できる場合だけ，判別式>=0をparameter domainとして
    // exactに保持する。1自由parameterなら既存の一変数不等式solverで区間へ正規化する。
    if (ambientDomain == mathematics::NumericDomain::Real
        && conversion->polynomial.degree() == 2) {
        const Expr a = conversion->polynomial.coefficient(2);
        if (proveNonZero(a, builtins, mathematics) != mathematics::TruthValue::True)
            return std::nullopt;

        const Expr b = conversion->polynomial.coefficient(1);
        const Expr c = conversion->polynomial.coefficient(0);
        const Expr discriminant = symbolicQuadraticDiscriminant(
            a, b, c, builtins, mathematics, angles);
        const bool centeredQuadratic = proveZero(b, builtins, mathematics)
            == mathematics::TruthValue::True;
        Expr realRootArgument = discriminant;
        if (centeredQuadratic) {
            Expr minusC = simplify(
                Expr::call(builtins.symbol(BuiltinId::Negate), {c}),
                builtins, mathematics, angles);
            realRootArgument = simplify(
                symbolic::expandExpression(
                    Expr::call(builtins.symbol(BuiltinId::Divide), {std::move(minusC), a}),
                    builtins, mathematics, angles, {256}),
                builtins, mathematics, angles);
        }

        std::vector<SolverVariable> resultVariables;
        resultVariables.reserve(variables.size());
        std::vector<SolverVariable> freeVariables;
        freeVariables.reserve(variables.size() - 1);
        for (const expression::Symbol& variable : variables) {
            resultVariables.push_back(SolverVariable{variable, ambientDomain});
            if (!(variable == *solvedVariable))
                freeVariables.push_back(SolverVariable{variable, mathematics::NumericDomain::Real});
        }

        mathematics::AssumptionSet realParameterAssumptions;
        for (const SolverVariable& parameter : freeVariables)
            realParameterAssumptions.add(mathematics::elementOf(
                Expr{parameter.symbol}, mathematics::NumericDomain::Real));
        const mathematics::KnowledgeContext realKnowledge{
            builtins, mathematics, realParameterAssumptions};
        const mathematics::Predicate nonnegative = mathematics::relation(
            RelationKind::GreaterEqual, realRootArgument, zeroExpr());
        const mathematics::TruthValue nonnegativeTruth = realKnowledge.prove(nonnegative);
        if (nonnegativeTruth == mathematics::TruthValue::False)
            return SolutionSet::empty(std::move(resultVariables));

        std::vector<mathematics::AssumptionSet> parameterRegions;
        if (nonnegativeTruth == mathematics::TruthValue::True) {
            parameterRegions.emplace_back();
        }
        else if (freeVariables.size() == 1) {
            Expr relation = Expr::call(
                builtins.symbol(BuiltinId::GreaterEqual), {realRootArgument, zeroExpr()});
            const SolutionSet range = solveUnivariatePolynomialRelation(
                relation, freeVariables.front().symbol, builtins, mathematics, angles);
            if (range.kind() == SolutionSetKind::Empty)
                return SolutionSet::empty(std::move(resultVariables));
            if (range.kind() == SolutionSetKind::Universal) {
                parameterRegions.emplace_back();
            }
            else if (range.kind() == SolutionSetKind::Finite) {
                for (const SolutionBranch& rangeBranch : range.branches()) {
                    // >= の一変数solverはfree-variable region branchを返す。
                    // point binding等へ変わる将来実装では，このbridgeで推測parameterizationをしない。
                    if (!rangeBranch.bindings.empty())
                        return std::nullopt;
                    parameterRegions.push_back(rangeBranch.conditions);
                }
            }
            else {
                parameterRegions.push_back(mathematics::AssumptionSet{std::vector<mathematics::Predicate>{nonnegative}});
            }
        }
        else {
            parameterRegions.push_back(mathematics::AssumptionSet{std::vector<mathematics::Predicate>{nonnegative}});
        }

        const auto discriminantZero = proveZero(realRootArgument, builtins, mathematics);
        const auto centeredRoot = [&](bool plus) {
            Expr root = simplify(
                Expr::call(builtins.symbol(BuiltinId::Sqrt), {realRootArgument}),
                builtins, mathematics, angles);
            if (plus)
                return root;
            return simplify(
                Expr::call(builtins.symbol(BuiltinId::Negate), {std::move(root)}),
                builtins, mathematics, angles);
        };
        std::vector<SolutionBranch> branches;
        const auto appendRoot = [&](Expr value, const mathematics::AssumptionSet& region,
                                    std::optional<std::size_t> multiplicity = std::nullopt) {
            SolutionBranch branch;
            branch.bindings.push_back(SolutionBinding{*solvedVariable, std::move(value)});
            branch.freeVariables = freeVariables;
            branch.multiplicity = multiplicity;
            branch.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
            for (const mathematics::Predicate& predicate : conversion->domainConditions.predicates())
                branch.conditions.add(predicate);
            for (const mathematics::Predicate& predicate : region.predicates())
                branch.conditions.add(predicate);
            if (std::find(branches.begin(), branches.end(), branch) == branches.end())
                branches.push_back(std::move(branch));
        };

        for (const mathematics::AssumptionSet& region : parameterRegions) {
            if (discriminantZero == mathematics::TruthValue::True) {
                appendRoot(
                    symbolicRepeatedQuadraticRoot(a, b, builtins, mathematics, angles),
                    region, std::size_t{2});
                continue;
            }
            appendRoot(centeredQuadratic
                ? centeredRoot(true)
                : symbolicQuadraticRoot(
                    a, b, discriminant, true, builtins, mathematics, angles), region);
            appendRoot(centeredQuadratic
                ? centeredRoot(false)
                : symbolicQuadraticRoot(
                    a, b, discriminant, false, builtins, mathematics, angles), region);
        }

        if (branches.empty())
            return SolutionSet::empty(std::move(resultVariables));
        return SolutionSet::finite(std::move(resultVariables), std::move(branches));
    }

    const SolutionSet projected = solveExpressionPolynomial(
        conversion->polynomial, *solvedVariable, builtins, mathematics, angles);
    if (projected.kind() == SolutionSetKind::Unresolved)
        return std::nullopt;

    std::vector<SolverVariable> resultVariables;
    resultVariables.reserve(variables.size());
    for (const expression::Symbol& variable : variables)
        resultVariables.push_back(SolverVariable{variable, ambientDomain});

    auto freeParameters = [&](bool includeSolved) {
        std::vector<SolverVariable> result;
        for (const expression::Symbol& variable : variables)
            if (includeSolved || !(variable == *solvedVariable))
                result.push_back(SolverVariable{variable, ambientDomain});
        return result;
    };

    std::vector<SolutionBranch> branches;
    auto appendProjectedBranch = [&](SolutionBranch branch,
                                     const mathematics::AssumptionSet& caseConditions) {
        for (const mathematics::Predicate& predicate : conversion->domainConditions.predicates())
            branch.conditions.add(predicate);
        for (const mathematics::Predicate& predicate : caseConditions.predicates())
            branch.conditions.add(predicate);
        branch.freeVariables = freeParameters(false);
        if (ambientDomain == mathematics::NumericDomain::Real
            && conversion->polynomial.degree() <= 1)
            branch.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
        if (std::find(branches.begin(), branches.end(), branch) == branches.end())
            branches.push_back(std::move(branch));
    };
    auto appendUniversalBranch = [&](const mathematics::AssumptionSet& caseConditions) {
        SolutionBranch branch;
        branch.freeVariables = freeParameters(true);
        for (const mathematics::Predicate& predicate : conversion->domainConditions.predicates())
            branch.conditions.add(predicate);
        for (const mathematics::Predicate& predicate : caseConditions.predicates())
            branch.conditions.add(predicate);
        if (std::find(branches.begin(), branches.end(), branch) == branches.end())
            branches.push_back(std::move(branch));
    };

    switch (projected.kind()) {
    case SolutionSetKind::Empty:
        return SolutionSet::empty(std::move(resultVariables));
    case SolutionSetKind::Universal:
        appendUniversalBranch(projected.conditions());
        break;
    case SolutionSetKind::Finite:
        for (const SolutionBranch& branch : projected.branches())
            appendProjectedBranch(branch, projected.conditions());
        break;
    case SolutionSetKind::Conditional:
        for (const SolutionCase& item : projected.cases()) {
            if (item.outcome == SolutionSetKind::Unresolved
                || item.outcome == SolutionSetKind::Conditional)
                return std::nullopt;
            if (item.outcome == SolutionSetKind::Empty)
                continue;
            if (item.outcome == SolutionSetKind::Universal) {
                appendUniversalBranch(item.conditions);
                continue;
            }
            for (const SolutionBranch& branch : item.branches)
                appendProjectedBranch(branch, item.conditions);
        }
        break;
    case SolutionSetKind::Unresolved:
        return std::nullopt;
    }

    if (branches.empty())
        return SolutionSet::empty(std::move(resultVariables));
    return SolutionSet::finite(std::move(resultVariables), std::move(branches));
}

[[nodiscard]] std::optional<SolutionSet> solveFactoredPolynomialSystem(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> variables,
    mathematics::NumericDomain ambientDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    constexpr std::size_t maximumFactorBranches = 16;
    for (std::size_t equationIndex = 0; equationIndex < equations.size(); ++equationIndex) {
        if (relationKindOf(equations[equationIndex], builtins) != RelationKind::Equal)
            continue;
        const Expr zeroForm = equationZeroForm(
            equations[equationIndex], builtins, mathematics, angles);
        Expr factored = zeroForm;
        if (!isHead(factored, builtins, BuiltinId::Multiply)) {
            const auto polynomial = symbolic::toMultivariateRationalPolynomial(
                zeroForm, builtins, symbolic::PolynomialConversionOptions{16, 64});
            if (polynomial && polynomial->totalDegree() <= 16 && polynomial->termCount() <= 64)
                factored = symbolic::factorExpression(
                    zeroForm, builtins, mathematics, angles);
        }
        if (!isHead(factored, builtins, BuiltinId::Multiply)
            || factored.asCall().arguments.size() < 2
            || factored.asCall().arguments.size() > maximumFactorBranches)
            continue;

        std::vector<Expr> factors;
        for (const Expr& factor : factored.asCall().arguments) {
            if (!containsAnySolverVariable(factor, variables)) {
                if (proveNonZero(factor, builtins, mathematics) == mathematics::TruthValue::True)
                    continue;
                return std::nullopt;
            }
            if (std::find(factors.begin(), factors.end(), factor) == factors.end())
                factors.push_back(factor);
        }
        if (factors.size() < 2)
            continue;

        std::vector<SolverVariable> resultVariables;
        resultVariables.reserve(variables.size());
        for (const expression::Symbol& variable : variables)
            resultVariables.push_back(SolverVariable{variable, ambientDomain});
        std::vector<SolutionBranch> unionBranches;

        for (const Expr& factor : factors) {
            std::vector<Expr> branchEquations;
            branchEquations.reserve(equations.size());
            for (std::size_t i = 0; i < equations.size(); ++i) {
                if (i == equationIndex)
                    continue;
                branchEquations.push_back(equations[i]);
            }
            branchEquations.push_back(Expr::call(
                builtins.symbol(BuiltinId::Equal), {factor, zeroExpr()}));

            SolutionSet branchSolution = solvePolynomialSubsystem(
                branchEquations, variables, ambientDomain,
                builtins, mathematics, angles);
            if (branchSolution.kind() == SolutionSetKind::Unresolved
                || branchSolution.kind() == SolutionSetKind::Conditional)
                return std::nullopt;
            if (branchSolution.kind() == SolutionSetKind::Universal)
                return SolutionSet::universal(std::move(resultVariables));
            if (branchSolution.kind() == SolutionSetKind::Empty)
                continue;
            for (const SolutionBranch& branch : branchSolution.branches())
                if (std::find(unionBranches.begin(), unionBranches.end(), branch)
                    == unionBranches.end())
                    unionBranches.push_back(branch);
        }

        if (unionBranches.empty())
            return SolutionSet::empty(std::move(resultVariables));
        return SolutionSet::finite(std::move(resultVariables), std::move(unionBranches));
    }
    return std::nullopt;
}

[[nodiscard]] std::optional<SolutionSet> solvePolynomialSystemInOrder(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> orderedSymbols,
    std::span<const expression::Symbol> outputSymbols,
    mathematics::NumericDomain ambientDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (equations.empty() || orderedSymbols.size() < 2)
        return std::nullopt;

    std::vector<SolverVariable> variables;
    variables.reserve(outputSymbols.size());
    for (const expression::Symbol& variable : outputSymbols)
        variables.push_back(SolverVariable{variable, ambientDomain});

    const symbolic::PolynomialRing ring{
        std::vector<expression::Symbol>{orderedSymbols.begin(), orderedSymbols.end()},
        symbolic::MonomialOrder::Lex};
    std::vector<symbolic::MultivariateRationalPolynomial> generators;
    generators.reserve(equations.size());
    std::vector<Expr> zeroForms;
    zeroForms.reserve(equations.size());
    bool nonlinear = false;

    for (const Expr& equation : equations) {
        const auto relation = relationKindOf(equation, builtins);
        if (!relation || *relation != RelationKind::Equal)
            return std::nullopt;
        Expr zeroForm = equationZeroForm(equation, builtins, mathematics, angles);
        const auto polynomial = symbolic::toMultivariateRationalPolynomial(
            zeroForm, ring, builtins,
            symbolic::PolynomialConversionOptions{4096, 100'000});
        if (!polynomial)
            return std::nullopt;
        nonlinear = nonlinear || polynomial->totalDegree() > 1;
        generators.push_back(*polynomial);
        zeroForms.push_back(std::move(zeroForm));
    }
    if (!nonlinear)
        return std::nullopt;

    symbolic::GroebnerComputation computation;
    try {
        computation = symbolic::groebnerBasis(generators, ring);
    }
    catch (const std::length_error&) {
        return SolutionSet::unresolved(std::move(variables));
    }

    if (computation.basis.size() == 1) {
        const auto leading = computation.basis.front().leadingTerm(ring);
        if (leading && leading->monomial.isOne())
            return SolutionSet::empty(std::move(variables));
    }

    const expression::Symbol& lastVariable = orderedSymbols.back();
    const symbolic::MultivariateRationalPolynomial* eliminant = nullptr;
    for (const auto& polynomial : computation.basis) {
        if (polynomial.degree(lastVariable) == 0)
            continue;
        bool onlyLast = true;
        for (const expression::Symbol& symbol : polynomial.variables())
            onlyLast = onlyLast && symbol == lastVariable;
        if (onlyLast && (!eliminant
                || polynomial.degree(lastVariable) < eliminant->degree(lastVariable)))
            eliminant = &polynomial;
    }
    if (!eliminant)
        return SolutionSet::unresolved(std::move(variables));

    const Expr eliminantExpr = symbolic::polynomialToExpandedExpr(*eliminant, builtins);
    const Expr eliminantEquation = Expr::call(
        builtins.symbol(BuiltinId::Equal), {eliminantExpr, zeroExpr()});
    bool lastSolutionsCertifiedReal = false;
    SolutionSet lastSolutions = [&] {
        if (ambientDomain == mathematics::NumericDomain::Real) {
            SolutionSet explicitRoots = solvePolynomialEquation(
                eliminantEquation, lastVariable, builtins, mathematics, angles);
            SolveConstraints realConstraint;
            realConstraint.domain = mathematics::NumericDomain::Real;
            explicitRoots = applySolveConstraints(
                std::move(explicitRoots), realConstraint, builtins, mathematics, angles);
            const bool explicitComplete = explicitRoots.kind() == SolutionSetKind::Finite
                && std::all_of(
                    explicitRoots.branches().begin(), explicitRoots.branches().end(),
                    [](const SolutionBranch& branch) {
                        return branch.unconditional() && branch.freeVariables.empty()
                            && branch.bindings.size() == 1;
                    });
            if (explicitComplete) {
                lastSolutionsCertifiedReal = true;
                return explicitRoots;
            }
            if (auto real = solveRealAlgebraicPolynomialEquation(
                    eliminantEquation, lastVariable, builtins, mathematics, angles)) {
                lastSolutionsCertifiedReal = true;
                return std::move(*real);
            }
        }
        return solvePolynomialEquation(
            eliminantEquation, lastVariable, builtins, mathematics, angles);
    }();
    if (lastSolutions.kind() != SolutionSetKind::Finite)
        return SolutionSet::unresolved(std::move(variables));

    // Shape-position relationをまず最後の変数の多項式として構築する。
    // Root式へ代入してからsimplifyで0判定すると，minimal polynomialの知識が
    // 必要になるため，公開候補を作る前にQ[t]上の剰余として一括検証する。
    std::vector<SolutionBinding> shapeBindings;
    shapeBindings.push_back(SolutionBinding{lastVariable, Expr{lastVariable}});
    for (std::size_t reverse = orderedSymbols.size() - 1; reverse-- > 0;) {
        const expression::Symbol& variable = orderedSymbols[reverse];
        const symbolic::MultivariateRationalPolynomial* relationPolynomial = nullptr;
        Rational variableCoefficient{BigInt{0}};

        for (const auto& polynomial : computation.basis) {
            if (polynomial.degree(variable) != 1)
                continue;
            bool containsEarlier = false;
            bool validVariableTerm = false;
            Rational coefficient{BigInt{0}};
            for (const symbolic::PolynomialTerm& term : polynomial.terms()) {
                for (std::size_t earlier = 0; earlier < reverse; ++earlier)
                    if (term.monomial.exponentOf(orderedSymbols[earlier]) != 0)
                        containsEarlier = true;
                const std::size_t exponent = term.monomial.exponentOf(variable);
                if (exponent == 0)
                    continue;
                if (exponent != 1 || term.monomial.totalDegree() != 1) {
                    containsEarlier = true;
                    break;
                }
                coefficient += term.coefficient;
                validVariableTerm = true;
            }
            if (!containsEarlier && validVariableTerm && !coefficient.isZero()) {
                relationPolynomial = &polynomial;
                variableCoefficient = coefficient;
                break;
            }
        }
        if (!relationPolynomial) {
            if (auto specialized = solveByEliminantSpecialization(
                    equations, orderedSymbols, outputSymbols, ambientDomain,
                    lastVariable, lastSolutions, lastSolutionsCertifiedReal,
                    builtins, mathematics, angles))
                return specialized;
            return SolutionSet::unresolved(std::move(variables));
        }

        std::vector<symbolic::PolynomialTerm> restTerms;
        for (const symbolic::PolynomialTerm& term : relationPolynomial->terms()) {
            if (term.monomial.exponentOf(variable) == 1
                && term.monomial.totalDegree() == 1)
                continue;
            restTerms.push_back(term);
        }
        Expr rest = symbolic::polynomialToExpandedExpr(
            symbolic::MultivariateRationalPolynomial{std::move(restTerms)}, builtins);
        for (const SolutionBinding& binding : shapeBindings)
            rest = symbolic::substituteSymbol(rest, binding.variable, binding.value);
        Expr value = simplify(
            mathematics::scaleExactExpression(
                -Rational{BigInt{1}} / variableCoefficient,
                std::move(rest), builtins),
            builtins, mathematics, angles);
        shapeBindings.push_back(SolutionBinding{variable, std::move(value)});
    }

    const symbolic::PolynomialRing eliminantRing{{lastVariable}, symbolic::MonomialOrder::Lex};
    const auto eliminantInOneVariable = symbolic::toMultivariateRationalPolynomial(
        eliminantExpr, eliminantRing, builtins,
        symbolic::PolynomialConversionOptions{4096, 100'000});
    if (!eliminantInOneVariable)
        return SolutionSet::unresolved(std::move(variables));
    const std::array<symbolic::MultivariateRationalPolynomial, 1> divisorBasis{
        *eliminantInOneVariable};

    for (const Expr& zeroForm : zeroForms) {
        Expr verified = zeroForm;
        for (const SolutionBinding& binding : shapeBindings) {
            if (binding.variable == lastVariable)
                continue;
            verified = symbolic::substituteSymbol(verified, binding.variable, binding.value);
        }
        verified = simplify(std::move(verified), builtins, mathematics, angles);
        const auto polynomial = symbolic::toMultivariateRationalPolynomial(
            verified, eliminantRing, builtins,
            symbolic::PolynomialConversionOptions{4096, 100'000});
        if (!polynomial)
            return SolutionSet::unresolved(std::move(variables));
        if (!symbolic::normalForm(*polynomial, divisorBasis, eliminantRing).isZero())
            return SolutionSet::unresolved(std::move(variables));
    }

    std::vector<SolutionBranch> resultBranches;
    resultBranches.reserve(lastSolutions.branches().size());
    for (const SolutionBranch& lastBranch : lastSolutions.branches()) {
        if (lastBranch.bindings.size() != 1
            || !(lastBranch.bindings.front().variable == lastVariable))
            return SolutionSet::unresolved(std::move(variables));

        const Expr& root = lastBranch.bindings.front().value;
        std::vector<SolutionBinding> bindings;
        bindings.reserve(outputSymbols.size());
        for (const expression::Symbol& variable : outputSymbols) {
            if (variable == lastVariable) {
                bindings.push_back(SolutionBinding{variable, root});
                continue;
            }
            const auto found = std::find_if(
                shapeBindings.begin(), shapeBindings.end(),
                [&](const SolutionBinding& binding) { return binding.variable == variable; });
            if (found == shapeBindings.end())
                return SolutionSet::unresolved(std::move(variables));
            Expr value = symbolic::substituteSymbol(found->value, lastVariable, root);
            value = simplify(std::move(value), builtins, mathematics, angles);
            bindings.push_back(SolutionBinding{variable, std::move(value)});
        }

        SolutionBranch branchResult;
        branchResult.bindings = std::move(bindings);
        branchResult.bindingsCertifiedDomain = lastSolutionsCertifiedReal
            ? mathematics::NumericDomain::Real
            : mathematics::NumericDomain::Complex;
        resultBranches.push_back(std::move(branchResult));
    }

    if (resultBranches.empty())
        return SolutionSet::empty(std::move(variables));
    return SolutionSet::finite(std::move(variables), std::move(resultBranches));
}



[[nodiscard]] std::optional<SolutionSet> solveByConstantLinearElimination(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> variables,
    mathematics::NumericDomain ambientDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (equations.size() < 2 || variables.size() < 2
        || (ambientDomain != mathematics::NumericDomain::Complex
            && ambientDomain != mathematics::NumericDomain::Real))
        return std::nullopt;

    constexpr std::size_t maximumEliminationVariables = 6;
    if (variables.size() > maximumEliminationVariables)
        return std::nullopt;

    for (std::size_t equationIndex = 0; equationIndex < equations.size(); ++equationIndex) {
        if (relationKindOf(equations[equationIndex], builtins) != RelationKind::Equal)
            continue;

        const Expr zeroForm = equationZeroForm(
            equations[equationIndex], builtins, mathematics, angles);
        const auto multivariate = symbolic::toMultivariateRationalPolynomial(
            zeroForm, builtins, symbolic::PolynomialConversionOptions{64, 512});
        if (!multivariate || multivariate->totalDegree() == 0)
            continue;

        for (auto variableIterator = variables.rbegin(); variableIterator != variables.rend(); ++variableIterator) {
            const expression::Symbol variable = *variableIterator;
            if (multivariate->degree(variable) != 1)
                continue;

            const auto conversion = symbolic::toExpressionPolynomialWithConditions(
                zeroForm, variable, builtins, mathematics, angles,
                symbolic::PolynomialConversionOptions{1, 512});
            if (!conversion || conversion->polynomial.degree() != 1
                || !conversion->domainConditions.empty())
                continue;

            const Expr coefficient = conversion->polynomial.coefficient(1);
            bool coefficientDependsOnSolverVariable = false;
            for (const expression::Symbol& candidate : variables)
                if (symbolic::containsSymbol(coefficient, candidate)) {
                    coefficientDependsOnSolverVariable = true;
                    break;
                }
            if (coefficientDependsOnSolverVariable
                || proveNonZero(coefficient, builtins, mathematics) != mathematics::TruthValue::True)
                continue;

            Expr negatedConstant = simplify(
                Expr::call(builtins.symbol(BuiltinId::Negate), {
                    conversion->polynomial.coefficient(0)}),
                builtins, mathematics, angles);
            Expr solvedValue = simplify(
                symbolic::expandExpression(
                    Expr::call(builtins.symbol(BuiltinId::Divide), {
                        std::move(negatedConstant), coefficient}),
                    builtins, mathematics, angles, {512}),
                builtins, mathematics, angles);

            std::vector<expression::Symbol> remainingVariables;
            remainingVariables.reserve(variables.size() - 1);
            for (const expression::Symbol& candidate : variables)
                if (!(candidate == variable))
                    remainingVariables.push_back(candidate);

            std::vector<Expr> reducedEquations;
            reducedEquations.reserve(equations.size() - 1);
            for (std::size_t i = 0; i < equations.size(); ++i) {
                if (i == equationIndex)
                    continue;
                Expr reduced = symbolic::substituteSymbol(equations[i], variable, solvedValue);
                reducedEquations.push_back(simplify(
                    std::move(reduced), builtins, mathematics, angles));
            }

            SolutionSet remaining = solvePolynomialSubsystem(
                reducedEquations, remainingVariables, ambientDomain,
                builtins, mathematics, angles);
            if (remaining.kind() == SolutionSetKind::Unresolved
                || remaining.kind() == SolutionSetKind::Conditional)
                continue;

            std::vector<SolverVariable> resultVariables;
            resultVariables.reserve(variables.size());
            for (const expression::Symbol& candidate : variables)
                resultVariables.push_back(SolverVariable{candidate, ambientDomain});

            if (remaining.kind() == SolutionSetKind::Empty)
                return SolutionSet::empty(std::move(resultVariables));

            auto combineBranch = [&](const SolutionBranch* subBranch) -> SolutionBranch {
                SolutionBranch combined;
                if (subBranch) {
                    combined.bindings.assign(
                        subBranch->bindings.begin(), subBranch->bindings.end());
                    combined.freeVariables.assign(
                        subBranch->freeVariables.begin(), subBranch->freeVariables.end());
                    combined.conditions = subBranch->conditions;
                    for (const mathematics::Predicate& predicate : remaining.conditions().predicates())
                        combined.conditions.add(predicate);
                }
                else {
                    for (const expression::Symbol& freeVariable : remainingVariables)
                        combined.freeVariables.push_back(SolverVariable{freeVariable, ambientDomain});
                    combined.conditions = remaining.conditions();
                }

                Expr specializedValue = solvedValue;
                if (subBranch)
                    for (const SolutionBinding& binding : subBranch->bindings)
                        specializedValue = symbolic::substituteSymbol(
                            specializedValue, binding.variable, binding.value);
                specializedValue = simplify(
                    std::move(specializedValue), builtins, mathematics, angles);
                combined.bindings.push_back(SolutionBinding{variable, std::move(specializedValue)});

                std::vector<SolutionBinding> orderedBindings;
                orderedBindings.reserve(combined.bindings.size());
                for (const expression::Symbol& outputVariable : variables) {
                    const auto found = std::find_if(
                        combined.bindings.begin(), combined.bindings.end(),
                        [&](const SolutionBinding& binding) {
                            return binding.variable == outputVariable;
                        });
                    if (found != combined.bindings.end())
                        orderedBindings.push_back(*found);
                }
                combined.bindings = std::move(orderedBindings);

                if (ambientDomain == mathematics::NumericDomain::Real
                    && (!subBranch || (subBranch->bindingsCertifiedDomain
                        && mathematics::isSubdomainOf(
                            *subBranch->bindingsCertifiedDomain,
                            mathematics::NumericDomain::Real))))
                    combined.bindingsCertifiedDomain = mathematics::NumericDomain::Real;
                return combined;
            };

            std::vector<SolutionBranch> branches;
            if (remaining.kind() == SolutionSetKind::Universal)
                branches.push_back(combineBranch(nullptr));
            else {
                branches.reserve(remaining.branches().size());
                for (const SolutionBranch& branch : remaining.branches())
                    branches.push_back(combineBranch(&branch));
            }

            return branches.empty()
                ? std::optional<SolutionSet>{SolutionSet::empty(std::move(resultVariables))}
                : std::optional<SolutionSet>{SolutionSet::finite(
                    std::move(resultVariables), std::move(branches))};
        }
    }

    return std::nullopt;
}

std::optional<SolutionSet> solvePolynomialSystem(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> variableSymbols,
    mathematics::NumericDomain ambientDomain,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (equations.empty() || variableSymbols.size() < 2)
        return std::nullopt;

    const auto primary = solvePolynomialSystemInOrder(
        equations, variableSymbols, variableSymbols, ambientDomain, builtins, mathematics, angles);
    if (!primary || primary->kind() != SolutionSetKind::Unresolved)
        return primary;

    // 2変数lex Groebner基底は，射影変数の選択だけでshape-position復元に失敗することがある。
    // 反対側の射影を1回だけ試し，解集合を変えずに定数係数の一次逆代入関係を探す。
    if (variableSymbols.size() == 2) {
        const std::array<expression::Symbol, 2> swapped{
            variableSymbols[1], variableSymbols[0]};
        const auto alternate = solvePolynomialSystemInOrder(
            equations, swapped, variableSymbols, ambientDomain, builtins, mathematics, angles);
        if (alternate && alternate->kind() != SolutionSetKind::Unresolved)
            return alternate;
    }

    // Gröbnerのzero-dimensional経路を先に尊重し，既存のcanonical Root/shape表示を
    // 奪わない。一般経路で完全解を構成できなかった場合だけ，定数係数の線形変数を
    // exactに消去して低次元systemへ落とす。その後に積=0のbranch分解と
    // 一方の変数を自由parameterにする低次projectionを試す。
    if (auto eliminated = solveByConstantLinearElimination(
            equations, variableSymbols, ambientDomain, builtins, mathematics, angles))
        return eliminated;
    if (auto factored = solveFactoredPolynomialSystem(
            equations, variableSymbols, ambientDomain, builtins, mathematics, angles))
        return factored;
    if (auto projected = solveSinglePolynomialProjection(
            equations, variableSymbols, ambientDomain, builtins, mathematics, angles))
        return projected;

    return primary;
}

SolutionSet solveLinearPolynomialSystem(
    std::span<const Expr> equations,
    std::span<const expression::Symbol> variableSymbols,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    bool hasOrderedRelation = false;
    for (const Expr& equation : equations) {
        const auto relation = relationKindOf(equation, builtins);
        if (relation && (*relation == RelationKind::Less
                || *relation == RelationKind::LessEqual
                || *relation == RelationKind::Greater
                || *relation == RelationKind::GreaterEqual)) {
            hasOrderedRelation = true;
            break;
        }
    }

    std::vector<SolverVariable> variables;
    variables.reserve(variableSymbols.size());
    for (const expression::Symbol& variable : variableSymbols)
        variables.push_back(SolverVariable{
            variable,
            hasOrderedRelation
                ? mathematics::NumericDomain::Real
                : mathematics::NumericDomain::Complex});

    // 多変数の半代数集合はまだ線形等式Gauss-Jordanとは別問題。
    // ordered relationを等式へ落として解くのは数学的に誤りなので、現段階では明示的に
    // Unresolvedとする。一変数の不等式はsolveUnivariatePolynomialRelationが扱う。
    if (hasOrderedRelation)
        return SolutionSet::unresolved(variables);
    if (variableSymbols.empty())
        return solveSymbolicLinearSystem(equations, variableSymbols, builtins, mathematics, angles);

    const std::size_t rows = equations.size();
    const std::size_t columns = variableSymbols.size();
    std::vector<std::vector<Rational>> matrix(
        rows, std::vector<Rational>(columns + 1, rational(0)));

    for (std::size_t row = 0; row < rows; ++row) {
        const Expr zeroForm = equationZeroForm(equations[row], builtins, mathematics, angles);
        const auto polynomial = symbolic::toMultivariateRationalPolynomial(zeroForm, builtins);
        if (!polynomial || polynomial->totalDegree() > 1)
            return solveSymbolicLinearSystem(equations, variableSymbols, builtins, mathematics, angles);

        for (const expression::Symbol& symbol : polynomial->variables())
            if (!sameVariable(symbol, variableSymbols))
                return solveSymbolicLinearSystem(equations, variableSymbols, builtins, mathematics, angles);

        Rational constant{BigInt{0}};
        for (const symbolic::PolynomialTerm& term : polynomial->terms()) {
            if (term.monomial.isOne()) {
                constant += term.coefficient;
                continue;
            }
            if (term.monomial.totalDegree() != 1 || term.monomial.factors().size() != 1)
                return solveSymbolicLinearSystem(equations, variableSymbols, builtins, mathematics, angles);
            const expression::Symbol& symbol = term.monomial.factors().front().variable;
            const auto iterator = std::find(variableSymbols.begin(), variableSymbols.end(), symbol);
            if (iterator == variableSymbols.end())
                return solveSymbolicLinearSystem(equations, variableSymbols, builtins, mathematics, angles);
            matrix[row][static_cast<std::size_t>(iterator - variableSymbols.begin())] += term.coefficient;
        }
        matrix[row][columns] = -constant;
    }

    // exact Rational Gauss-Jordan elimination。
    std::size_t pivotRow = 0;
    std::vector<std::optional<std::size_t>> pivotForColumn(columns);
    for (std::size_t column = 0; column < columns && pivotRow < rows; ++column) {
        std::size_t selected = pivotRow;
        while (selected < rows && matrix[selected][column].isZero())
            ++selected;
        if (selected == rows)
            continue;
        if (selected != pivotRow)
            std::swap(matrix[selected], matrix[pivotRow]);

        const Rational pivot = matrix[pivotRow][column];
        for (std::size_t j = column; j <= columns; ++j)
            matrix[pivotRow][j] /= pivot;

        for (std::size_t row = 0; row < rows; ++row) {
            if (row == pivotRow || matrix[row][column].isZero())
                continue;
            const Rational scale = matrix[row][column];
            for (std::size_t j = column; j <= columns; ++j)
                matrix[row][j] -= scale * matrix[pivotRow][j];
        }
        pivotForColumn[column] = pivotRow;
        ++pivotRow;
    }

    for (std::size_t row = 0; row < rows; ++row) {
        bool allZero = true;
        for (std::size_t column = 0; column < columns; ++column)
            allZero = allZero && matrix[row][column].isZero();
        if (allZero && !matrix[row][columns].isZero())
            return SolutionSet::empty(variables);
    }

    if (pivotRow == 0)
        return SolutionSet::universal(variables);

    SolutionBranch solution;
    std::vector<std::size_t> freeColumns;
    for (std::size_t column = 0; column < columns; ++column) {
        if (!pivotForColumn[column]) {
            freeColumns.push_back(column);
            solution.freeVariables.push_back(variables[column]);
        }
    }

    solution.bindings.reserve(columns - freeColumns.size());
    for (std::size_t column = 0; column < columns; ++column) {
        if (!pivotForColumn[column])
            continue;
        const std::size_t row = *pivotForColumn[column];
        Expr value{Number{matrix[row][columns]}};
        for (const std::size_t freeColumn : freeColumns) {
            if (matrix[row][freeColumn].isZero())
                continue;
            Expr term = mathematics::scaleExactExpression(
                -matrix[row][freeColumn], Expr{variableSymbols[freeColumn]}, builtins);
            value = Expr::call(builtins.symbol(BuiltinId::Add), {std::move(value), std::move(term)});
        }
        solution.bindings.push_back(SolutionBinding{
            variableSymbols[column],
            simplify(std::move(value), builtins, mathematics, angles)});
    }
    return SolutionSet::finite(variables, {std::move(solution)});
}

} // namespace mmcal::solver
