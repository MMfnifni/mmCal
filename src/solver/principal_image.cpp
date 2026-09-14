// principal branchの像判定をSolver間で共有する。
#include "principal_image.hpp"

#include "evaluation/builtin_registry.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/knowledge_context.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "solver_support.hpp"

#include <algorithm>
#include <utility>

namespace mmcal::solver {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::AssumptionSet;
using mathematics::FunctionBranchRule;
using mathematics::Predicate;
using mathematics::RelationKind;
using mathematics::TruthValue;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

using Region = std::vector<Predicate>;

[[nodiscard]] Expr piMultiple(
    const Rational& coefficient,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions) {
    const auto* pi = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!pi)
        return integerExpr(0);
    if (coefficient == Rational{BigInt{1}})
        return Expr{pi->symbol};
    if (coefficient == Rational{BigInt{-1}})
        return simplifyForSolve(
            Expr::call(builtins.symbol(BuiltinId::Negate), {Expr{pi->symbol}}),
            builtins, mathematics, angles, assumptions);
    return simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Multiply), {
            Expr{Number{coefficient}}, Expr{pi->symbol}}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] Expr angleValueFromTurns(
    const Rational& turns,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Degree:
        return Expr{Number{turns * Rational{BigInt{360}}}};
    case mathematics::AngleUnit::Gradian:
        return Expr{Number{turns * Rational{BigInt{400}}}};
    case mathematics::AngleUnit::Radian:
        return piMultiple(turns * Rational{BigInt{2}}, builtins, mathematics, angles, assumptions);
    }
    return integerExpr(0);
}

struct CartesianParts final {
    Expr real;
    Expr imaginary;
};

[[nodiscard]] std::optional<CartesianParts> exactCartesianParts(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions,
    std::size_t depth = 0) {
    if (depth > 64)
        return std::nullopt;
    if (expression.isNumber())
        return CartesianParts{
            Expr{Number{expression.asNumber().realPart()}},
            Expr{Number{expression.asNumber().imaginaryPart()}}};

    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    if (knowledge.prove(mathematics::elementOf(
            expression, mathematics::NumericDomain::Real)) == TruthValue::True)
        return CartesianParts{expression, integerExpr(0)};
    if (!expression.isCall())
        return std::nullopt;

    const auto& arguments = expression.asCall().arguments;
    const auto simplify = [&](Expr value) {
        return simplifyForSolve(
            std::move(value), builtins, mathematics, angles, assumptions);
    };
    const auto addParts = [&](CartesianParts lhs, CartesianParts rhs) {
        return CartesianParts{
            simplify(Expr::call(builtins.symbol(BuiltinId::Add), {
                std::move(lhs.real), std::move(rhs.real)})),
            simplify(Expr::call(builtins.symbol(BuiltinId::Add), {
                std::move(lhs.imaginary), std::move(rhs.imaginary)}))};
    };

    if (builtins.isCallTo(expression, BuiltinId::Negate) && arguments.size() == 1) {
        auto value = exactCartesianParts(
            arguments[0], builtins, mathematics, angles, assumptions, depth + 1);
        if (!value)
            return std::nullopt;
        return CartesianParts{
            simplify(Expr::call(
                builtins.symbol(BuiltinId::Negate), {std::move(value->real)})),
            simplify(Expr::call(
                builtins.symbol(BuiltinId::Negate), {std::move(value->imaginary)}))};
    }
    if (builtins.isCallTo(expression, BuiltinId::Add)) {
        CartesianParts result{integerExpr(0), integerExpr(0)};
        for (const Expr& argument : arguments) {
            auto value = exactCartesianParts(
                argument, builtins, mathematics, angles, assumptions, depth + 1);
            if (!value)
                return std::nullopt;
            result = addParts(std::move(result), std::move(*value));
        }
        return result;
    }
    if (builtins.isCallTo(expression, BuiltinId::Subtract) && arguments.size() == 2) {
        auto lhs = exactCartesianParts(
            arguments[0], builtins, mathematics, angles, assumptions, depth + 1);
        auto rhs = exactCartesianParts(
            arguments[1], builtins, mathematics, angles, assumptions, depth + 1);
        if (!lhs || !rhs)
            return std::nullopt;
        rhs->real = simplify(Expr::call(
            builtins.symbol(BuiltinId::Negate), {std::move(rhs->real)}));
        rhs->imaginary = simplify(Expr::call(
            builtins.symbol(BuiltinId::Negate), {std::move(rhs->imaginary)}));
        return addParts(std::move(*lhs), std::move(*rhs));
    }
    if (builtins.isCallTo(expression, BuiltinId::Divide) && arguments.size() == 2) {
        auto lhs = exactCartesianParts(
            arguments[0], builtins, mathematics, angles, assumptions, depth + 1);
        auto rhs = exactCartesianParts(
            arguments[1], builtins, mathematics, angles, assumptions, depth + 1);
        if (!lhs || !rhs)
            return std::nullopt;
        Expr denominator = simplify(Expr::call(builtins.symbol(BuiltinId::Add), {
            simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {rhs->real, rhs->real})),
            simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
                rhs->imaginary, rhs->imaginary}))}));
        Expr ac = simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
            lhs->real, rhs->real}));
        Expr bd = simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
            lhs->imaginary, rhs->imaginary}));
        Expr bc = simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
            lhs->imaginary, rhs->real}));
        Expr ad = simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
            lhs->real, rhs->imaginary}));
        return CartesianParts{
            simplify(Expr::call(builtins.symbol(BuiltinId::Divide), {
                simplify(Expr::call(builtins.symbol(BuiltinId::Add), {
                    std::move(ac), std::move(bd)})), denominator})),
            simplify(Expr::call(builtins.symbol(BuiltinId::Divide), {
                simplify(Expr::call(builtins.symbol(BuiltinId::Subtract), {
                    std::move(bc), std::move(ad)})), std::move(denominator)}))};
    }
    if (builtins.isCallTo(expression, BuiltinId::Multiply)) {
        CartesianParts result{integerExpr(1), integerExpr(0)};
        for (const Expr& argument : arguments) {
            auto rhs = exactCartesianParts(
                argument, builtins, mathematics, angles, assumptions, depth + 1);
            if (!rhs)
                return std::nullopt;
            Expr ac = simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
                result.real, rhs->real}));
            Expr bd = simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
                result.imaginary, rhs->imaginary}));
            Expr ad = simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
                result.real, rhs->imaginary}));
            Expr bc = simplify(Expr::call(builtins.symbol(BuiltinId::Multiply), {
                result.imaginary, rhs->real}));
            result.real = simplify(Expr::call(
                builtins.symbol(BuiltinId::Subtract), {std::move(ac), std::move(bd)}));
            result.imaginary = simplify(Expr::call(
                builtins.symbol(BuiltinId::Add), {std::move(ad), std::move(bc)}));
        }
        return result;
    }
    return std::nullopt;
}

[[nodiscard]] Expr coordinate(
    const Expr& value,
    BuiltinId component,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions) {
    return simplifyForSolve(
        Expr::call(builtins.symbol(component), {value}),
        builtins, mathematics, angles, assumptions);
}

[[nodiscard]] TruthValue relationTruthFromOrder(
    RelationKind relation,
    CertifiedOrder order) {
    if (order == CertifiedOrder::Unknown)
        return TruthValue::Unknown;
    switch (relation) {
    case RelationKind::Less:
        return order == CertifiedOrder::Less ? TruthValue::True : TruthValue::False;
    case RelationKind::LessEqual:
        return order == CertifiedOrder::Less || order == CertifiedOrder::Equal
            ? TruthValue::True : TruthValue::False;
    case RelationKind::Greater:
        return order == CertifiedOrder::Greater ? TruthValue::True : TruthValue::False;
    case RelationKind::GreaterEqual:
        return order == CertifiedOrder::Greater || order == CertifiedOrder::Equal
            ? TruthValue::True : TruthValue::False;
    case RelationKind::Equal:
        return order == CertifiedOrder::Equal ? TruthValue::True : TruthValue::False;
    case RelationKind::NotEqual:
        return order == CertifiedOrder::Equal ? TruthValue::False : TruthValue::True;
    }
    return TruthValue::Unknown;
}

[[nodiscard]] Predicate simplifyPredicate(
    const Predicate& predicate,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions) {
    if (const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate))
        return mathematics::relation(
            relation->relation,
            simplifyForSolve(relation->lhs, builtins, mathematics, angles, assumptions),
            simplifyForSolve(relation->rhs, builtins, mathematics, angles, assumptions));
    const auto& domain = std::get<mathematics::DomainPredicate>(predicate);
    return mathematics::elementOf(
        simplifyForSolve(domain.expression, builtins, mathematics, angles, assumptions),
        domain.domain);
}

[[nodiscard]] TruthValue provePredicate(
    const Predicate& predicate,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions) {
    const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
    TruthValue truth = knowledge.prove(predicate);
    if (truth != TruthValue::Unknown)
        return truth;

    const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate);
    if (!relation)
        return TruthValue::Unknown;

    const Expr difference = simplifyForSolve(
        Expr::call(builtins.symbol(BuiltinId::Subtract), {relation->lhs, relation->rhs}),
        builtins, mathematics, angles, assumptions);
    if (difference.isNumber() && difference.asNumber().isReal()
        && difference.asNumber().isZero()) {
        switch (relation->relation) {
        case RelationKind::Less:
        case RelationKind::Greater:
        case RelationKind::NotEqual:
            return TruthValue::False;
        case RelationKind::LessEqual:
        case RelationKind::GreaterEqual:
        case RelationKind::Equal:
            return TruthValue::True;
        }
    }

    return relationTruthFromOrder(
        relation->relation,
        certifiedConstantOrder(
            relation->lhs, relation->rhs, builtins, mathematics, angles));
}

[[nodiscard]] std::optional<AssumptionSet> evaluateRegion(
    const Region& region,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions) {
    AssumptionSet local = assumptions;
    AssumptionSet remaining;
    for (const Predicate& raw : region) {
        const Predicate predicate = simplifyPredicate(
            raw, builtins, mathematics, angles, local);
        const TruthValue truth = provePredicate(
            predicate, builtins, mathematics, angles, local);
        if (truth == TruthValue::False)
            return std::nullopt;
        if (truth == TruthValue::Unknown) {
            remaining.add(predicate);
            local.add(predicate);
        }
    }
    return remaining;
}

void addComplexDomain(Region& region, const Expr& value) {
    region.push_back(mathematics::elementOf(
        value, mathematics::NumericDomain::Complex));
}

void addRelation(Region& region, RelationKind kind, Expr lhs, Expr rhs) {
    region.push_back(mathematics::relation(kind, std::move(lhs), std::move(rhs)));
}

[[nodiscard]] std::vector<Region> principalComplexRegions(
    FunctionBranchRule rule,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions) {
    const Expr zero = integerExpr(0);
    const auto cartesian = exactCartesianParts(
        value, builtins, mathematics, angles, assumptions);
    const Expr re = cartesian ? cartesian->real : coordinate(
        value, BuiltinId::Re, builtins, mathematics, angles, assumptions);
    const Expr im = cartesian ? cartesian->imaginary : coordinate(
        value, BuiltinId::Im, builtins, mathematics, angles, assumptions);

    const Expr negativePi = piMultiple(
        Rational{BigInt{-1}}, builtins, mathematics, angles, assumptions);
    const Expr pi = piMultiple(
        Rational{BigInt{1}}, builtins, mathematics, angles, assumptions);
    const Expr negativeHalfPi = piMultiple(
        Rational{BigInt{-1}, BigInt{2}}, builtins, mathematics, angles, assumptions);
    const Expr halfPi = piMultiple(
        Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles, assumptions);
    const Expr negativeQuarterTurn = angleValueFromTurns(
        Rational{BigInt{-1}, BigInt{4}}, builtins, mathematics, angles, assumptions);
    const Expr quarterTurn = angleValueFromTurns(
        Rational{BigInt{1}, BigInt{4}}, builtins, mathematics, angles, assumptions);
    const Expr halfTurn = angleValueFromTurns(
        Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles, assumptions);

    std::vector<Region> regions;
    switch (rule) {
    case FunctionBranchRule::PrincipalSquareRoot: {
        // sqrtの像は右半平面と正の虚軸。実値と証明済みなら単一条件へ縮約する。
        const mathematics::KnowledgeContext knowledge{builtins, mathematics, assumptions};
        if (knowledge.prove(mathematics::elementOf(
                value, mathematics::NumericDomain::Real)) == TruthValue::True) {
            Region realRegion;
            addRelation(realRegion, RelationKind::GreaterEqual, value, zero);
            regions.push_back(std::move(realRegion));
            break;
        }
        Region interior;
        addRelation(interior, RelationKind::Greater, re, zero);
        regions.push_back(std::move(interior));
        Region boundary;
        addRelation(boundary, RelationKind::Equal, re, zero);
        addRelation(boundary, RelationKind::GreaterEqual, im, zero);
        regions.push_back(std::move(boundary));
        break;
    }
    case FunctionBranchRule::PrincipalLogarithm: {
        Region region;
        addComplexDomain(region, value);
        addRelation(region, RelationKind::Greater, im, negativePi);
        addRelation(region, RelationKind::LessEqual, im, pi);
        regions.push_back(std::move(region));
        break;
    }
    case FunctionBranchRule::PrincipalArcSine: {
        Region interior;
        addComplexDomain(interior, value);
        addRelation(interior, RelationKind::Greater, re, negativeQuarterTurn);
        addRelation(interior, RelationKind::Less, re, quarterTurn);
        regions.push_back(std::move(interior));

        Region left;
        addComplexDomain(left, value);
        addRelation(left, RelationKind::Equal, re, negativeQuarterTurn);
        addRelation(left, RelationKind::GreaterEqual, im, zero);
        regions.push_back(std::move(left));

        Region right;
        addComplexDomain(right, value);
        addRelation(right, RelationKind::Equal, re, quarterTurn);
        addRelation(right, RelationKind::LessEqual, im, zero);
        regions.push_back(std::move(right));
        break;
    }
    case FunctionBranchRule::PrincipalArcCosine: {
        Region interior;
        addComplexDomain(interior, value);
        addRelation(interior, RelationKind::Greater, re, zero);
        addRelation(interior, RelationKind::Less, re, halfTurn);
        regions.push_back(std::move(interior));

        Region left;
        addComplexDomain(left, value);
        addRelation(left, RelationKind::Equal, re, zero);
        addRelation(left, RelationKind::GreaterEqual, im, zero);
        regions.push_back(std::move(left));

        Region right;
        addComplexDomain(right, value);
        addRelation(right, RelationKind::Equal, re, halfTurn);
        addRelation(right, RelationKind::LessEqual, im, zero);
        regions.push_back(std::move(right));
        break;
    }
    case FunctionBranchRule::PrincipalArcTangent: {
        Region interior;
        addComplexDomain(interior, value);
        addRelation(interior, RelationKind::Greater, re, negativeQuarterTurn);
        addRelation(interior, RelationKind::Less, re, quarterTurn);
        regions.push_back(std::move(interior));

        Region left;
        addComplexDomain(left, value);
        addRelation(left, RelationKind::Equal, re, negativeQuarterTurn);
        addRelation(left, RelationKind::Less, im, zero);
        regions.push_back(std::move(left));

        Region right;
        addComplexDomain(right, value);
        addRelation(right, RelationKind::Equal, re, quarterTurn);
        addRelation(right, RelationKind::Greater, im, zero);
        regions.push_back(std::move(right));
        break;
    }
    case FunctionBranchRule::PrincipalAreaHyperbolicSine: {
        Region interior;
        addComplexDomain(interior, value);
        addRelation(interior, RelationKind::Greater, im, negativeHalfPi);
        addRelation(interior, RelationKind::Less, im, halfPi);
        regions.push_back(std::move(interior));

        Region lower;
        addComplexDomain(lower, value);
        addRelation(lower, RelationKind::Equal, im, negativeHalfPi);
        addRelation(lower, RelationKind::LessEqual, re, zero);
        regions.push_back(std::move(lower));

        Region upper;
        addComplexDomain(upper, value);
        addRelation(upper, RelationKind::Equal, im, halfPi);
        addRelation(upper, RelationKind::GreaterEqual, re, zero);
        regions.push_back(std::move(upper));
        break;
    }
    case FunctionBranchRule::PrincipalAreaHyperbolicCosine: {
        Region interior;
        addComplexDomain(interior, value);
        addRelation(interior, RelationKind::Greater, re, zero);
        addRelation(interior, RelationKind::Greater, im, negativePi);
        addRelation(interior, RelationKind::LessEqual, im, pi);
        regions.push_back(std::move(interior));

        Region boundary;
        addComplexDomain(boundary, value);
        addRelation(boundary, RelationKind::Equal, re, zero);
        addRelation(boundary, RelationKind::GreaterEqual, im, zero);
        addRelation(boundary, RelationKind::LessEqual, im, pi);
        regions.push_back(std::move(boundary));
        break;
    }
    case FunctionBranchRule::PrincipalAreaHyperbolicTangent: {
        Region interior;
        addComplexDomain(interior, value);
        addRelation(interior, RelationKind::Greater, im, negativeHalfPi);
        addRelation(interior, RelationKind::Less, im, halfPi);
        regions.push_back(std::move(interior));

        Region lower;
        addComplexDomain(lower, value);
        addRelation(lower, RelationKind::Equal, im, negativeHalfPi);
        addRelation(lower, RelationKind::Greater, re, zero);
        regions.push_back(std::move(lower));

        Region upper;
        addComplexDomain(upper, value);
        addRelation(upper, RelationKind::Equal, im, halfPi);
        addRelation(upper, RelationKind::Less, re, zero);
        regions.push_back(std::move(upper));
        break;
    }
    default:
        break;
    }
    return regions;
}


[[nodiscard]] std::vector<Region> principalRealInputRegions(
    FunctionBranchRule rule,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions) {
    const Expr zero = integerExpr(0);
    const auto cartesian = exactCartesianParts(
        value, builtins, mathematics, angles, assumptions);
    const Expr re = cartesian ? cartesian->real : coordinate(
        value, BuiltinId::Re, builtins, mathematics, angles, assumptions);
    const Expr im = cartesian ? cartesian->imaginary : coordinate(
        value, BuiltinId::Im, builtins, mathematics, angles, assumptions);

    const Expr pi = piMultiple(
        Rational{BigInt{1}}, builtins, mathematics, angles, assumptions);
    const Expr negativeHalfPi = piMultiple(
        Rational{BigInt{-1}, BigInt{2}}, builtins, mathematics, angles, assumptions);
    const Expr halfPi = piMultiple(
        Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles, assumptions);
    const Expr negativeQuarterTurn = angleValueFromTurns(
        Rational{BigInt{-1}, BigInt{4}}, builtins, mathematics, angles, assumptions);
    const Expr quarterTurn = angleValueFromTurns(
        Rational{BigInt{1}, BigInt{4}}, builtins, mathematics, angles, assumptions);
    const Expr halfTurn = angleValueFromTurns(
        Rational{BigInt{1}, BigInt{2}}, builtins, mathematics, angles, assumptions);

    std::vector<Region> regions;
    switch (rule) {
    case FunctionBranchRule::PrincipalSquareRoot: {
        // x>=0 は非負実軸，x<0 は正の虚軸へ写る。
        Region realAxis;
        realAxis.push_back(mathematics::elementOf(
            value, mathematics::NumericDomain::Real));
        addRelation(realAxis, RelationKind::GreaterEqual, value, zero);
        regions.push_back(std::move(realAxis));

        Region imaginaryAxis;
        addComplexDomain(imaginaryAxis, value);
        addRelation(imaginaryAxis, RelationKind::Equal, re, zero);
        addRelation(imaginaryAxis, RelationKind::Greater, im, zero);
        regions.push_back(std::move(imaginaryAxis));
        break;
    }
    case FunctionBranchRule::PrincipalLogarithm: {
        // principal Log(real nonzero) は実軸（正入力）または Im=+Pi（負入力）。
        Region positiveInputImage;
        positiveInputImage.push_back(mathematics::elementOf(
            value, mathematics::NumericDomain::Real));
        regions.push_back(std::move(positiveInputImage));

        Region negativeInputImage;
        addComplexDomain(negativeInputImage, value);
        addRelation(negativeInputImage, RelationKind::Equal, im, pi);
        regions.push_back(std::move(negativeInputImage));
        break;
    }
    case FunctionBranchRule::PrincipalArcSine: {
        // real x in [-1,1] は実区間，|x|>1 は両端から伸びる縦半直線へ写る。
        Region realSegment;
        realSegment.push_back(mathematics::elementOf(
            value, mathematics::NumericDomain::Real));
        addRelation(realSegment, RelationKind::GreaterEqual, value, negativeQuarterTurn);
        addRelation(realSegment, RelationKind::LessEqual, value, quarterTurn);
        regions.push_back(std::move(realSegment));

        Region leftRay;
        addComplexDomain(leftRay, value);
        addRelation(leftRay, RelationKind::Equal, re, negativeQuarterTurn);
        addRelation(leftRay, RelationKind::Greater, im, zero);
        regions.push_back(std::move(leftRay));

        Region rightRay;
        addComplexDomain(rightRay, value);
        addRelation(rightRay, RelationKind::Equal, re, quarterTurn);
        addRelation(rightRay, RelationKind::Less, im, zero);
        regions.push_back(std::move(rightRay));
        break;
    }
    case FunctionBranchRule::PrincipalArcCosine: {
        Region realSegment;
        realSegment.push_back(mathematics::elementOf(
            value, mathematics::NumericDomain::Real));
        addRelation(realSegment, RelationKind::GreaterEqual, value, zero);
        addRelation(realSegment, RelationKind::LessEqual, value, halfTurn);
        regions.push_back(std::move(realSegment));

        Region zeroRay;
        addComplexDomain(zeroRay, value);
        addRelation(zeroRay, RelationKind::Equal, re, zero);
        addRelation(zeroRay, RelationKind::Greater, im, zero);
        regions.push_back(std::move(zeroRay));

        Region piRay;
        addComplexDomain(piRay, value);
        addRelation(piRay, RelationKind::Equal, re, halfTurn);
        addRelation(piRay, RelationKind::Less, im, zero);
        regions.push_back(std::move(piRay));
        break;
    }
    case FunctionBranchRule::PrincipalArcTangent: {
        Region realSegment;
        realSegment.push_back(mathematics::elementOf(
            value, mathematics::NumericDomain::Real));
        addRelation(realSegment, RelationKind::Greater, value, negativeQuarterTurn);
        addRelation(realSegment, RelationKind::Less, value, quarterTurn);
        regions.push_back(std::move(realSegment));
        break;
    }
    case FunctionBranchRule::PrincipalAreaHyperbolicSine: {
        Region realAxis;
        realAxis.push_back(mathematics::elementOf(
            value, mathematics::NumericDomain::Real));
        regions.push_back(std::move(realAxis));
        break;
    }
    case FunctionBranchRule::PrincipalAreaHyperbolicCosine: {
        // x>=1: 非負実軸，-1<=x<=1: 正の虚軸，x<-1: Im=Pi の上側境界。
        Region realRay;
        realRay.push_back(mathematics::elementOf(
            value, mathematics::NumericDomain::Real));
        addRelation(realRay, RelationKind::GreaterEqual, value, zero);
        regions.push_back(std::move(realRay));

        Region imaginarySegment;
        addComplexDomain(imaginarySegment, value);
        addRelation(imaginarySegment, RelationKind::Equal, re, zero);
        addRelation(imaginarySegment, RelationKind::Greater, im, zero);
        addRelation(imaginarySegment, RelationKind::Less, im, pi);
        regions.push_back(std::move(imaginarySegment));

        Region upperBoundary;
        addComplexDomain(upperBoundary, value);
        addRelation(upperBoundary, RelationKind::Equal, im, pi);
        addRelation(upperBoundary, RelationKind::GreaterEqual, re, zero);
        regions.push_back(std::move(upperBoundary));
        break;
    }
    case FunctionBranchRule::PrincipalAreaHyperbolicTangent: {
        // |x|<1 は実軸，x>1 / x<-1 はそれぞれ下/上の半開境界へ写る。
        Region realAxis;
        realAxis.push_back(mathematics::elementOf(
            value, mathematics::NumericDomain::Real));
        regions.push_back(std::move(realAxis));

        Region lower;
        addComplexDomain(lower, value);
        addRelation(lower, RelationKind::Equal, im, negativeHalfPi);
        addRelation(lower, RelationKind::Greater, re, zero);
        regions.push_back(std::move(lower));

        Region upper;
        addComplexDomain(upper, value);
        addRelation(upper, RelationKind::Equal, im, halfPi);
        addRelation(upper, RelationKind::Less, re, zero);
        regions.push_back(std::move(upper));
        break;
    }
    default:
        break;
    }
    return regions;
}

[[nodiscard]] std::vector<Region> principalRegions(
    FunctionBranchRule rule,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions,
    PrincipalImageInputDomain inputDomain) {
    if (inputDomain == PrincipalImageInputDomain::Real)
        return principalRealInputRegions(
            rule, value, builtins, mathematics, angles, assumptions);
    return principalComplexRegions(
        rule, value, builtins, mathematics, angles, assumptions);
}

[[nodiscard]] bool supportedRule(FunctionBranchRule rule) noexcept {
    switch (rule) {
    case FunctionBranchRule::PrincipalSquareRoot:
    case FunctionBranchRule::PrincipalLogarithm:
    case FunctionBranchRule::PrincipalArcSine:
    case FunctionBranchRule::PrincipalArcCosine:
    case FunctionBranchRule::PrincipalArcTangent:
    case FunctionBranchRule::PrincipalAreaHyperbolicSine:
    case FunctionBranchRule::PrincipalAreaHyperbolicCosine:
    case FunctionBranchRule::PrincipalAreaHyperbolicTangent:
        return true;
    default:
        return false;
    }
}

} // namespace

std::optional<PrincipalImageAnalysis> analyzePrincipalImage(
    FunctionBranchRule branchRule,
    const Expr& value,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles,
    const AssumptionSet& assumptions,
    PrincipalImageInputDomain inputDomain) {
    if (!supportedRule(branchRule))
        return std::nullopt;

    // principal sqrtの出力を再度sqrtの像判定へ掛ける場合は，定義そのものから像内。
    if (inputDomain == PrincipalImageInputDomain::Complex
        && branchRule == FunctionBranchRule::PrincipalSquareRoot
        && builtins.isCallTo(value, BuiltinId::Sqrt))
        return PrincipalImageAnalysis{{AssumptionSet{}}};

    const auto regions = principalRegions(
        branchRule, value, builtins, mathematics, angles, assumptions, inputDomain);
    PrincipalImageAnalysis result;
    for (const Region& region : regions) {
        auto remaining = evaluateRegion(
            region, builtins, mathematics, angles, assumptions);
        if (!remaining)
            continue;
        if (remaining->empty())
            return PrincipalImageAnalysis{{AssumptionSet{}}};
        if (std::find(result.alternatives.begin(), result.alternatives.end(), *remaining)
            == result.alternatives.end())
            result.alternatives.push_back(std::move(*remaining));
    }
    return result;
}

} // namespace mmcal::solver
