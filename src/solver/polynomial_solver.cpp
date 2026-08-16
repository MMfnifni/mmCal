// 多項式方程式solver
#include "polynomial_solver.hpp"

#include "mathematics/exact_algebra.hpp"
#include "mathematics/math_ids.hpp"
#include "mathematics/knowledge_context.hpp"
#include "mathematics/predicate.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "solver/solve_constraints.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/algebraic_number.hpp"
#include "symbolic/algebra_transforms.hpp"

#include <algorithm>
#include <cstdint>
#include <iterator>
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

[[nodiscard]] Expr algebraicRootExpr(
    const symbolic::RealAlgebraicNumber& algebraic,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Rational> coefficients(
        algebraic.polynomial().begin(), algebraic.polynomial().end());
    const std::size_t coefficientCount = coefficients.size();
    return Expr::call(builtins.symbol(BuiltinId::Root), {
        Expr::rationalArray({coefficientCount}, std::move(coefficients)),
        Expr{Number{BigInt::fromUnsigned(algebraic.rootIndex())}}
    });
}

[[nodiscard]] Expr algebraicRootExpr(
    const symbolic::ComplexAlgebraicNumber& algebraic,
    const evaluation::BuiltinRegistry& builtins) {
    std::vector<Rational> coefficients(
        algebraic.polynomial().begin(), algebraic.polynomial().end());
    const std::size_t coefficientCount = coefficients.size();
    return Expr::call(builtins.symbol(BuiltinId::Root), {
        Expr::rationalArray({coefficientCount}, std::move(coefficients)),
        Expr{Number{BigInt::fromUnsigned(algebraic.rootIndex())}},
        Expr{expression::Symbol{"Complex"}}
    });
}

[[nodiscard]] SolutionBranch branch(
    const expression::Symbol& variable,
    Expr value,
    std::optional<std::size_t> multiplicity = std::nullopt) {
    return SolutionBranch{{SolutionBinding{variable, std::move(value)}}, {}, multiplicity, {}};
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
    switch (definition->id) {
    case BuiltinId::Equal: return RelationKind::Equal;
    case BuiltinId::NotEqual: return RelationKind::NotEqual;
    case BuiltinId::Less: return RelationKind::Less;
    case BuiltinId::LessEqual: return RelationKind::LessEqual;
    case BuiltinId::Greater: return RelationKind::Greater;
    case BuiltinId::GreaterEqual: return RelationKind::GreaterEqual;
    default: return std::nullopt;
    }
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
        {SolverVariable{variable, mathematics::NumericDomain::Real}}};
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

} // namespace

SolutionSet solveUnivariatePolynomialRelation(
    const Expr& relation,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
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
                {SolverVariable{variable, mathematics::NumericDomain::Complex}}}});
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

    const auto roots = symbolic::RealAlgebraicNumber::isolateAll(polynomial->coefficients());
    if (!roots)
        return std::nullopt;
    const std::vector<SolverVariable> variables{{variable, mathematics::NumericDomain::Real}};
    if (roots->empty())
        return SolutionSet::empty(variables);

    std::vector<SolutionBranch> branches;
    branches.reserve(roots->size());
    for (const symbolic::RealAlgebraicNumber& root : *roots)
        branches.push_back(branch(variable, algebraicRootExpr(root, builtins)));
    return SolutionSet::finite(variables, std::move(branches));
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
        std::vector<SolutionBranch> branches;
        branches.reserve(roots->size());
        for (const symbolic::ComplexAlgebraicNumber& root : *roots)
            branches.push_back(branch(variable, algebraicRootExpr(root, builtins)));
        return SolutionSet::finite(variables, std::move(branches));
    }
    return SolutionSet::unresolved(variables);
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
