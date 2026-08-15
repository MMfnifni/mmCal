// 式のcanonical順序
#include "expression_ordering.hpp"

#include "mathematics/predicate.hpp"
#include "solver/solution_set.hpp"

#include <cstdint>
#include <string>
#include <string_view>

namespace mmcal::simplification {
namespace {

using expression::Expr;


[[nodiscard]] std::string legacyOrderKey(const Expr& expression) {
    switch (expression.kind()) {
    case expression::ExprKind::Number:
        return "0:" + expression.asNumber().toString();
    case expression::ExprKind::DecimalApproximation:
        return "1:" + std::string{expression.asDecimalApproximation().text()};
    case expression::ExprKind::ComplexDecimalApproximation:
        return "2:" + std::string{expression.asComplexDecimalApproximation().text()};
    case expression::ExprKind::Boolean:
        return expression.asBoolean() ? "3:1" : "3:0";
    case expression::ExprKind::String:
        return "4:" + expression.asString();
    case expression::ExprKind::Symbol:
        return "5:" + expression.asSymbol().name();
    case expression::ExprKind::Array: {
        std::string key{"6:"};
        for (const std::size_t dimension : expression.asArray().shape)
            key += std::to_string(dimension) + ',';
        key.push_back(':');
        for (std::size_t i = 0; i < expression.asArray().size(); ++i)
            key += legacyOrderKey(expression.asArray().element(i)) + ';';
        return key;
    }
    case expression::ExprKind::List: {
        std::string key{"9:"};
        for (const Expr& element : expression.asList().elements)
            key += legacyOrderKey(element) + ';';
        return key;
    }
    case expression::ExprKind::Call: {
        std::string key = "7:" + expression.asCall().head.name() + '[';
        for (const Expr& argument : expression.asCall().arguments)
            key += legacyOrderKey(argument) + ';';
        key.push_back(']');
        return key;
    }
    case expression::ExprKind::SolutionSet:
        // V1.5.0までのcanonical表示順を一次キーとして維持する。
        // 実際のSolutionSet構造はcollision-freeな二次キーで完全順序化する。
        return "8:SolutionSet";
    }
    return {};
}

void appendField(std::string& output, std::string_view value) {
    output += std::to_string(value.size());
    output.push_back(':');
    output.append(value);
    output.push_back(';');
}

void appendUnsigned(std::string& output, std::size_t value) {
    appendField(output, std::to_string(value));
}

void appendRational(std::string& output, const numeric::Rational& value) {
    appendField(output, value.toString());
}

void appendExpr(std::string& output, const Expr& expression);

void appendPredicate(std::string& output, const mathematics::Predicate& predicate);

void appendAssumptions(
    std::string& output,
    const mathematics::AssumptionSet& assumptions) {
    appendUnsigned(output, assumptions.size());
    for (const auto& predicate : assumptions.predicates())
        appendPredicate(output, predicate);
}

void appendPredicate(std::string& output, const mathematics::Predicate& predicate) {
    if (const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate)) {
        appendField(output, "R");
        appendUnsigned(output, static_cast<std::size_t>(relation->relation));
        appendExpr(output, relation->lhs);
        appendExpr(output, relation->rhs);
        return;
    }

    const auto& domain = std::get<mathematics::DomainPredicate>(predicate);
    appendField(output, "D");
    appendUnsigned(output, static_cast<std::size_t>(domain.domain));
    appendExpr(output, domain.expression);
}

void appendSolverVariable(std::string& output, const solver::SolverVariable& variable) {
    appendField(output, variable.symbol.name());
    appendUnsigned(output, static_cast<std::size_t>(variable.domain));
}

void appendBranch(std::string& output, const solver::SolutionBranch& branch) {
    appendUnsigned(output, branch.bindings.size());
    for (const auto& binding : branch.bindings) {
        appendField(output, binding.variable.name());
        appendExpr(output, binding.value);
    }
    appendAssumptions(output, branch.conditions);
    appendField(output, branch.multiplicity
        ? std::to_string(*branch.multiplicity) : std::string{"none"});
    appendUnsigned(output, branch.freeVariables.size());
    for (const auto& variable : branch.freeVariables)
        appendSolverVariable(output, variable);
}

void appendSolutionSet(std::string& output, const solver::SolutionSet& solutions) {
    appendUnsigned(output, static_cast<std::size_t>(solutions.kind()));
    appendUnsigned(output, solutions.variables().size());
    for (const auto& variable : solutions.variables())
        appendSolverVariable(output, variable);

    appendUnsigned(output, solutions.branches().size());
    for (const auto& branch : solutions.branches())
        appendBranch(output, branch);

    appendUnsigned(output, solutions.cases().size());
    for (const auto& solutionCase : solutions.cases()) {
        appendAssumptions(output, solutionCase.conditions);
        appendUnsigned(output, static_cast<std::size_t>(solutionCase.outcome));
        appendUnsigned(output, solutionCase.branches.size());
        for (const auto& branch : solutionCase.branches)
            appendBranch(output, branch);
    }

    appendAssumptions(output, solutions.conditions());
}

void appendDecimal(std::string& output, const numeric::DecimalApproximation& value) {
    appendField(output, value.text());
    appendUnsigned(output, value.fractionalDigits());
    appendUnsigned(output, value.requestedFractionalDigits());
    appendField(output, value.isRounded() ? "1" : "0");
    appendField(output, value.origin() == numeric::ApproximationOrigin::ExactValue ? "0" : "1");
    appendRational(output, value.displayedValue());
    appendRational(output, value.certifiedLower());
    appendRational(output, value.certifiedUpper());
}

void appendExpr(std::string& output, const Expr& expression) {
    output.push_back(static_cast<char>('0' + static_cast<int>(expression.kind())));
    output.push_back('{');

    switch (expression.kind()) {
    case expression::ExprKind::Number:
        appendField(output, expression.asNumber().toString());
        break;
    case expression::ExprKind::DecimalApproximation:
        appendDecimal(output, expression.asDecimalApproximation());
        break;
    case expression::ExprKind::ComplexDecimalApproximation: {
        const auto& value = expression.asComplexDecimalApproximation();
        appendDecimal(output, value.real());
        appendDecimal(output, value.imaginary());
        appendField(output, value.realExactlyZero() ? "1" : "0");
        appendField(output, value.imaginaryExactlyZero() ? "1" : "0");
        break;
    }
    case expression::ExprKind::Boolean:
        appendField(output, expression.asBoolean() ? "1" : "0");
        break;
    case expression::ExprKind::String:
        appendField(output, expression.asString());
        break;
    case expression::ExprKind::Symbol:
        appendField(output, expression.asSymbol().name());
        break;
    case expression::ExprKind::Array:
        appendUnsigned(output, expression.asArray().shape.size());
        for (const std::size_t dimension : expression.asArray().shape)
            appendUnsigned(output, dimension);
        appendUnsigned(output, expression.asArray().size());
        for (std::size_t i = 0; i < expression.asArray().size(); ++i)
            appendExpr(output, expression.asArray().element(i));
        break;
    case expression::ExprKind::List:
        appendUnsigned(output, expression.asList().elements.size());
        for (const Expr& element : expression.asList().elements)
            appendExpr(output, element);
        break;
    case expression::ExprKind::Call:
        appendField(output, expression.asCall().head.name());
        appendUnsigned(output, expression.asCall().arguments.size());
        for (const Expr& argument : expression.asCall().arguments)
            appendExpr(output, argument);
        break;
    case expression::ExprKind::SolutionSet:
        appendSolutionSet(output, expression.asSolutionSet());
        break;
    }

    output.push_back('}');
}

} // namespace

std::string expressionOrderKey(const Expr& expression) {
    static constexpr char hex[] = "0123456789ABCDEF";
    const std::string primary = legacyOrderKey(expression);

    // primaryを固定幅hexへ写すと元のbyte列の辞書順を保ったままdelimiterを
    // 安全に挿入できる。従来順を維持し、同一primaryの衝突だけ完全構造でtie-breakする。
    std::string output;
    output.reserve(primary.size() * 2 + 65);
    for (const unsigned char byte : primary) {
        output.push_back(hex[byte >> 4]);
        output.push_back(hex[byte & 0x0f]);
    }
    output.push_back('/'); // hex digitより小さいため、primaryのprefix順も維持する。
    appendExpr(output, expression);
    return output;
}

bool ExpressionLess::operator()(const Expr& lhs, const Expr& rhs) const {
    if (lhs == rhs)
        return false;
    return expressionOrderKey(lhs) < expressionOrderKey(rhs);
}

} // namespace mmcal::simplification
