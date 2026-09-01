// 再parse可能な標準数式表示
#include "expr_formatter.hpp"

#include "builtins/names.hpp"
#include "expression/array_utils.hpp"
#include "numeric/integer_algorithms.hpp"
#include "solver/solution_set.hpp"
#include "mathematics/predicate.hpp"

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <variant>

namespace mmcal::formatting {
namespace {

using expression::CallExpr;
using expression::Expr;

constexpr int precedenceLowest = 0;
constexpr int precedenceAssignment = 5;
constexpr int precedenceComparison = 8;
constexpr int precedenceAdditive = 10;
constexpr int precedenceMultiplicative = 20;
constexpr int precedenceUnary = 30;
constexpr int precedencePower = 40;
constexpr int precedencePostfix = 50;

void appendExpr(std::string& output, const Expr& expression, unsigned radix, int parentPrecedence);

[[nodiscard]] int numberPrecedence(const numeric::Number& number) noexcept {
    // Numberは式木上はAtomだが、線形表示は "-8" や "1/3" のように
    // 演算子を含み得る。親の優先順位を無視すると、
    //     (-8)^(1/3) -> -8 ^ 1/3
    // のように再parse時の意味が変わるため、表示上の優先順位を別途与える。
    if (number.isComplex()) {
        // a+bI は加法式だが、純虚数 aI は暗黙乗算相当。Iを通常の積の中で
        // 不要に (I) とせず、Powerの底では (2I)^2 のように必要な括弧を付ける。
        return number.realPart().isZero()
            ? precedenceMultiplicative
            : precedenceAdditive;
    }

    const auto& real = number.asReal();
    if (!real.isInteger())
        return precedenceMultiplicative;
    if (real.isNegative())
        return precedenceUnary;
    return precedencePostfix + 1;
}

void appendNumber(
    std::string& output,
    const numeric::Number& number,
    unsigned radix,
    int parentPrecedence) {
    const int precedence = numberPrecedence(number);
    const bool parenthesize = precedence < parentPrecedence;
    if (parenthesize)
        output.push_back('(');
    output += number.toString(radix);
    if (parenthesize)
        output.push_back(')');
}

void appendEscapedString(std::string& output, std::string_view value) {
    output.push_back('"');

    for (const char character : value) {
        switch (character) {
        case '\\':
            output += "\\\\";
            break;
        case '"':
            output += "\\\"";
            break;
        case '\n':
            output += "\\n";
            break;
        case '\r':
            output += "\\r";
            break;
        case '\t':
            output += "\\t";
            break;
        default:
            output.push_back(character);
            break;
        }
    }

    output.push_back('"');
}


[[nodiscard]] bool appendRelationAgainstZeroInSolvedForm(
    std::string& output,
    const mathematics::RelationPredicate& relation,
    unsigned radix) {
    if (!relation.rhs.isNumber() || !relation.rhs.asNumber().isZero() || !relation.lhs.isCall())
        return false;

    const auto appendSeparator = [&]() {
        switch (relation.relation) {
        case mathematics::RelationKind::Equal: output += " == "; break;
        case mathematics::RelationKind::NotEqual: output += " != "; break;
        case mathematics::RelationKind::Less: output += " < "; break;
        case mathematics::RelationKind::LessEqual: output += " <= "; break;
        case mathematics::RelationKind::Greater: output += " > "; break;
        case mathematics::RelationKind::GreaterEqual: output += " >= "; break;
        }
    };

    const auto& call = relation.lhs.asCall();
    const std::string_view head = call.head.view();

    // x-c ? 0 は表示上 x ? c へ戻す。Solver内部ではzero-formの方が扱いやすいが、
    // ユーザーへ x-2 < 0 と見せる理由はない。これは表示だけの同値変形である。
    if (head == builtins::names::subtract && call.arguments.size() == 2
        && call.arguments[1].isNumber()) {
        appendExpr(output, call.arguments[0], radix, precedenceComparison);
        appendSeparator();
        appendExpr(output, call.arguments[1], radix, precedenceComparison + 1);
        return true;
    }

    if (head != builtins::names::add || call.arguments.size() < 2)
        return false;

    std::optional<std::size_t> numericIndex;
    for (std::size_t i = 0; i < call.arguments.size(); ++i) {
        if (!call.arguments[i].isNumber())
            continue;
        if (numericIndex)
            return false;
        numericIndex = i;
    }
    if (!numericIndex || call.arguments[*numericIndex].asNumber().isZero())
        return false;

    std::vector<Expr> remaining;
    remaining.reserve(call.arguments.size() - 1);
    for (std::size_t i = 0; i < call.arguments.size(); ++i)
        if (i != *numericIndex)
            remaining.push_back(call.arguments[i]);
    if (remaining.empty())
        return false;

    const Expr lhs = remaining.size() == 1
        ? remaining.front()
        : Expr::call(call.head, std::move(remaining));
    const Expr rhs{-call.arguments[*numericIndex].asNumber()};

    appendExpr(output, lhs, radix, precedenceComparison);
    appendSeparator();
    appendExpr(output, rhs, radix, precedenceComparison + 1);
    return true;
}

void appendPredicate(
    std::string& output,
    const mathematics::Predicate& predicate,
    unsigned radix) {
    if (const auto* relation = std::get_if<mathematics::RelationPredicate>(&predicate)) {
        if (appendRelationAgainstZeroInSolvedForm(output, *relation, radix))
            return;
        appendExpr(output, relation->lhs, radix, precedenceComparison);
        switch (relation->relation) {
        case mathematics::RelationKind::Equal: output += " == "; break;
        case mathematics::RelationKind::NotEqual: output += " != "; break;
        case mathematics::RelationKind::Less: output += " < "; break;
        case mathematics::RelationKind::LessEqual: output += " <= "; break;
        case mathematics::RelationKind::Greater: output += " > "; break;
        case mathematics::RelationKind::GreaterEqual: output += " >= "; break;
        }
        appendExpr(output, relation->rhs, radix, precedenceComparison + 1);
        return;
    }

    const auto& domain = std::get<mathematics::DomainPredicate>(predicate);
    appendExpr(output, domain.expression, radix, precedenceLowest);
    output += " in ";
    switch (domain.domain) {
    case mathematics::NumericDomain::Integer: output += "Integer"; break;
    case mathematics::NumericDomain::Rational: output += "Rational"; break;
    case mathematics::NumericDomain::Real: output += "Real"; break;
    case mathematics::NumericDomain::Complex: output += "Complex"; break;
    case mathematics::NumericDomain::Unknown: output += "Unknown"; break;
    }
}

void appendConditions(
    std::string& output,
    const mathematics::AssumptionSet& conditions,
    unsigned radix) {
    for (std::size_t i = 0; i < conditions.size(); ++i) {
        if (i != 0) output += " && ";
        appendPredicate(output, conditions.predicates()[i], radix);
    }
}

void appendSolutionBranches(
    std::string& output,
    std::span<const solver::SolutionBranch> branches,
    unsigned radix) {
    output.push_back('{');
    for (std::size_t i = 0; i < branches.size(); ++i) {
        if (i != 0) output += ", ";
        const auto& branch = branches[i];
        if (branch.bindings.size() > 1) output.push_back('{');
        for (std::size_t j = 0; j < branch.bindings.size(); ++j) {
            if (j != 0) output += ", ";
            output += branch.bindings[j].variable.name();
            output += " == ";
            appendExpr(output, branch.bindings[j].value, radix, precedenceComparison + 1);
        }
        if (branch.bindings.size() > 1) output.push_back('}');
        if (!branch.freeVariables.empty()) {
            // 不等式解のようにbindingを持たず、変数自身が自由parameterであるbranchは
            // "{ where x in Real ...}" ではなく "{x in Real if ...}" と表示する。
            const bool regionBranch = branch.bindings.empty();
            if (!regionBranch)
                output += " where ";
            for (std::size_t k = 0; k < branch.freeVariables.size(); ++k) {
                if (k != 0) output += ", ";
                output += branch.freeVariables[k].symbol.name();
                output += " in ";
                switch (branch.freeVariables[k].domain) {
                case mathematics::NumericDomain::Integer: output += "Integer"; break;
                case mathematics::NumericDomain::Rational: output += "Rational"; break;
                case mathematics::NumericDomain::Real: output += "Real"; break;
                case mathematics::NumericDomain::Complex: output += "Complex"; break;
                case mathematics::NumericDomain::Unknown: output += "Unknown"; break;
                }
            }
        }
        if (!branch.conditions.empty()) {
            output += " if ";
            appendConditions(output, branch.conditions, radix);
        }
        if (branch.multiplicity && *branch.multiplicity > 1)
            output += " (multiplicity " + std::to_string(*branch.multiplicity) + ")";
    }
    output.push_back('}');
}

void appendSolutionSet(
    std::string& output,
    const solver::SolutionSet& solutions,
    unsigned radix) {
    using solver::SolutionSetKind;
    if (solutions.kind() == SolutionSetKind::Empty) {
        output += "{}";
        if (!solutions.conditions().empty()) {
            output += " if ";
            appendConditions(output, solutions.conditions(), radix);
        }
        return;
    }
    if (solutions.kind() == SolutionSetKind::Universal) {
        output += "All";
        if (!solutions.conditions().empty()) {
            output += " if ";
            appendConditions(output, solutions.conditions(), radix);
        }
        return;
    }
    if (solutions.kind() == SolutionSetKind::Conditional) {
        output += "cases[";
        for (std::size_t i = 0; i < solutions.cases().size(); ++i) {
            if (i != 0) output += "; ";
            const auto& solutionCase = solutions.cases()[i];
            switch (solutionCase.outcome) {
            case SolutionSetKind::Empty:
                output += "{}";
                break;
            case SolutionSetKind::Finite:
                appendSolutionBranches(output, solutionCase.branches, radix);
                break;
            case SolutionSetKind::Universal:
                output += "All";
                break;
            case SolutionSetKind::Unresolved:
                output += "Unresolved";
                break;
            case SolutionSetKind::Conditional:
                output += "Unresolved";
                break;
            }
            if (!solutionCase.conditions.empty()) {
                output += " if ";
                appendConditions(output, solutionCase.conditions, radix);
            }
        }
        output += "]";
        if (!solutions.conditions().empty()) {
            output += " if ";
            appendConditions(output, solutions.conditions(), radix);
        }
        return;
    }
    if (solutions.kind() == SolutionSetKind::Unresolved) {
        output += "UnresolvedSolutionSet[";
        for (std::size_t i = 0; i < solutions.variables().size(); ++i) {
            if (i != 0) output += ", ";
            output += solutions.variables()[i].symbol.name();
        }
        output += "]";
        if (!solutions.conditions().empty()) {
            output += " if ";
            appendConditions(output, solutions.conditions(), radix);
        }
        return;
    }

    appendSolutionBranches(output, solutions.branches(), radix);
    if (!solutions.conditions().empty()) {
        output += " if ";
        appendConditions(output, solutions.conditions(), radix);
    }
}

[[nodiscard]] std::size_t trailingBlockSize(
    std::span<const std::size_t> shape,
    std::size_t dimension) {
    std::size_t blockSize = 1;

    for (std::size_t i = dimension + 1; i < shape.size(); ++i)
        blockSize *= shape[i];

    return blockSize;
}

void appendArrayDimension(
    std::string& output,
    const expression::ArrayExpr& array,
    std::size_t dimension,
    std::size_t offset,
    unsigned radix) {
    output.push_back('{');

    const std::size_t count = array.shape[dimension];
    const std::size_t blockSize = trailingBlockSize(array.shape, dimension);

    for (std::size_t i = 0; i < count; ++i) {
        if (i != 0)
            output += ", ";

        const std::size_t elementOffset = offset + i * blockSize;

        if (dimension + 1 == array.shape.size()) {
            switch (array.storedKindAt(elementOffset)) {
            case expression::ArrayStorageKind::Integer:
                output += array.integerAt(elementOffset).toString(radix);
                break;
            case expression::ArrayStorageKind::Rational:
                output += array.rationalAt(elementOffset).toString(radix);
                break;
            case expression::ArrayStorageKind::Number:
                appendNumber(output, array.numberAt(elementOffset), radix, precedenceLowest);
                break;
            case expression::ArrayStorageKind::DecimalApproximation:
                output += array.decimalAt(elementOffset).text();
                break;
            case expression::ArrayStorageKind::ComplexDecimalApproximation:
                output += array.complexDecimalAt(elementOffset).text();
                break;
            case expression::ArrayStorageKind::Generic:
                appendExpr(output, array.expressionAt(elementOffset), radix, precedenceLowest);
                break;
            }
        }
        else
            appendArrayDimension(output, array, dimension + 1, elementOffset, radix);
    }

    output.push_back('}');
}


void appendShapeLiteral(std::string& output, std::span<const std::size_t> shape) {
    output.push_back('{');
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (i != 0)
            output += ", ";
        output += std::to_string(shape[i]);
    }
    output.push_back('}');
}

void appendBinaryCall(
    std::string& output,
    const CallExpr& call,
    std::string_view separator,
    int precedence,
    unsigned radix,
    int parentPrecedence,
    bool protectRightOperand) {
    if (call.arguments.size() != 2)
        return;

    const bool parenthesize = precedence < parentPrecedence;
    if (parenthesize)
        output.push_back('(');

    appendExpr(output, call.arguments[0], radix, precedence);
    output += separator;
    appendExpr(
        output,
        call.arguments[1],
        radix,
        protectRightOperand ? precedence + 1 : precedence);

    if (parenthesize)
        output.push_back(')');
}

[[nodiscard]] std::optional<Expr> positiveMagnitudeOfNegative(const Expr& expression) {
    if (expression.isNumber() && expression.asNumber().isReal()
        && expression.asNumber().asReal().isNegative())
        return Expr{-expression.asNumber()};

    if (!expression.isCall())
        return std::nullopt;

    const auto& call = expression.asCall();
    const std::string_view head = call.head.view();
    if (head == builtins::names::negate && call.arguments.size() == 1)
        return call.arguments.front();

    if (head == builtins::names::multiply && !call.arguments.empty()) {
        const Expr& first = call.arguments.front();
        if (!first.isNumber() || !first.asNumber().isReal()
            || !first.asNumber().asReal().isNegative())
            return std::nullopt;

        std::vector<Expr> factors{call.arguments.begin(), call.arguments.end()};
        factors.front() = Expr{-first.asNumber()};
        if (factors.front().isNumber()
            && factors.front().asNumber().isReal()
            && factors.front().asNumber().asReal().isInteger()
            && factors.front().asNumber().asReal().asInteger() == numeric::BigInt{1})
            factors.erase(factors.begin());
        if (factors.size() == 1)
            return factors.front();
        return Expr::call(call.head, std::move(factors));
    }

    if (head == builtins::names::divide && call.arguments.size() == 2) {
        if (auto numerator = positiveMagnitudeOfNegative(call.arguments.front()))
            return Expr::call(call.head, {std::move(*numerator), call.arguments[1]});
    }

    return std::nullopt;
}

struct SignedTerm final {
    bool negative = false;
    Expr magnitude;
};

struct DisplayMonomial final {
    std::optional<symbols::SymbolId> variable;
    std::uint64_t degree = 0;
};

[[nodiscard]] std::optional<DisplayMonomial> displayMonomial(const Expr& expression) {
    if (expression.isNumber())
        return DisplayMonomial{};

    if (expression.isSymbol())
        return DisplayMonomial{expression.asSymbol().id(), 1};

    if (!expression.isCall())
        return std::nullopt;

    const auto& call = expression.asCall();
    const std::string_view head = call.head.view();
    if (head == builtins::names::power && call.arguments.size() == 2
        && call.arguments[0].isSymbol() && call.arguments[1].isNumber()) {
        const auto& exponent = call.arguments[1].asNumber();
        if (!exponent.isReal() || !exponent.asReal().isInteger()
            || exponent.asReal().isNegative())
            return std::nullopt;
        const auto degree = numeric::tryToUint64(exponent.asReal().asInteger());
        if (!degree)
            return std::nullopt;
        return DisplayMonomial{call.arguments[0].asSymbol().id(), *degree};
    }

    if (head != builtins::names::multiply || call.arguments.empty())
        return std::nullopt;

    DisplayMonomial result;
    for (const Expr& factor : call.arguments) {
        if (factor.isNumber())
            continue;
        const auto monomial = displayMonomial(factor);
        if (!monomial || !monomial->variable)
            return std::nullopt;
        if (result.variable && *result.variable != *monomial->variable)
            return std::nullopt;
        result.variable = monomial->variable;
        if (result.degree > std::numeric_limits<std::uint64_t>::max() - monomial->degree)
            return std::nullopt;
        result.degree += monomial->degree;
    }
    return result;
}

void orderUnivariatePolynomialTermsForDisplay(std::vector<SignedTerm>& terms) {
    if (terms.size() < 2)
        return;

    std::vector<std::uint64_t> degrees;
    degrees.reserve(terms.size());
    std::optional<symbols::SymbolId> variable;
    bool hasVariableTerm = false;
    for (const SignedTerm& term : terms) {
        const auto monomial = displayMonomial(term.magnitude);
        if (!monomial)
            return;
        if (monomial->variable) {
            if (variable && *variable != *monomial->variable)
                return;
            variable = monomial->variable;
            hasVariableTerm = true;
        }
        degrees.push_back(monomial->degree);
    }
    if (!hasVariableTerm)
        return;

    std::vector<std::size_t> order(terms.size());
    for (std::size_t i = 0; i < order.size(); ++i)
        order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](std::size_t lhs, std::size_t rhs) {
        return degrees[lhs] > degrees[rhs];
    });

    std::vector<SignedTerm> sorted;
    sorted.reserve(terms.size());
    for (const std::size_t index : order)
        sorted.push_back(std::move(terms[index]));
    terms = std::move(sorted);
}

void collectSignedTerms(
    const Expr& expression,
    bool negative,
    std::vector<SignedTerm>& terms) {
    if (expression.isCall()) {
        const auto& call = expression.asCall();
        const std::string_view head = call.head.view();
        if (head == builtins::names::add) {
            for (const Expr& argument : call.arguments)
                collectSignedTerms(argument, negative, terms);
            return;
        }
        if (head == builtins::names::subtract && call.arguments.size() == 2) {
            collectSignedTerms(call.arguments[0], negative, terms);
            collectSignedTerms(call.arguments[1], !negative, terms);
            return;
        }
    }

    if (auto magnitude = positiveMagnitudeOfNegative(expression)) {
        collectSignedTerms(*magnitude, !negative, terms);
        return;
    }

    terms.push_back(SignedTerm{negative, expression});
}

void appendSignedTerms(
    std::string& output,
    std::span<const SignedTerm> terms,
    unsigned radix,
    int parentPrecedence) {
    const bool parenthesize = precedenceAdditive < parentPrecedence;
    if (parenthesize)
        output.push_back('(');

    for (std::size_t i = 0; i < terms.size(); ++i) {
        const SignedTerm& term = terms[i];
        if (i == 0) {
            if (term.negative)
                output.push_back('-');
        }
        else
            output.push_back(term.negative ? '-' : '+');

        // '+'/'-' の右側は加法優先順位より一段強く表示する。
        // これにより複素数Atomの a+bI などだけに必要な括弧が付き、
        // 積や冪には不要な括弧を増やさない。
        const int termPrecedence = i == 0 && !term.negative
            ? precedenceAdditive
            : precedenceAdditive + 1;
        appendExpr(output, term.magnitude, radix, termPrecedence);
    }

    if (parenthesize)
        output.push_back(')');
}

void appendAddCall(
    std::string& output,
    const CallExpr& call,
    unsigned radix,
    int parentPrecedence) {
    std::vector<SignedTerm> terms;
    terms.reserve(call.arguments.size());
    for (const Expr& argument : call.arguments)
        collectSignedTerms(argument, false, terms);
    orderUnivariatePolynomialTermsForDisplay(terms);
    appendSignedTerms(output, terms, radix, parentPrecedence);
}

void appendMultiplyCall(
    std::string& output,
    const CallExpr& call,
    unsigned radix,
    int parentPrecedence) {
    const bool parenthesize = precedenceMultiplicative < parentPrecedence;
    if (parenthesize)
        output.push_back('(');

    std::string previousText;
    for (std::size_t i = 0; i < call.arguments.size(); ++i) {
        std::string currentText;
        appendExpr(currentText, call.arguments[i], radix, precedenceMultiplicative);

        if (i != 0 && !previousText.empty() && !currentText.empty()) {
            const unsigned char previousLast = static_cast<unsigned char>(previousText.back());
            const unsigned char currentFirst = static_cast<unsigned char>(currentText.front());
            const bool previousIdentifier = std::isalpha(previousLast) || previousText.back() == '_';
            const bool previousDigit = std::isdigit(previousLast) != 0;
            const bool currentIdentifier = std::isalpha(currentFirst) || currentText.front() == '_';
            const bool currentDigit = std::isdigit(currentFirst) != 0;
            const bool exponentMarkerCollision = previousDigit && currentIdentifier
                && (currentText.front() == 'e' || currentText.front() == 'E');
            const bool standaloneZeroSuffix = previousText.back() == '0'
                && (previousText.size() == 1
                    || (!std::isdigit(static_cast<unsigned char>(previousText[previousText.size() - 2]))
                        && previousText[previousText.size() - 2] != '.'));
            const bool radixPrefixCollision = standaloneZeroSuffix && currentIdentifier
                && (currentText.front() == 'b' || currentText.front() == 'B'
                    || currentText.front() == 'o' || currentText.front() == 'O'
                    || currentText.front() == 'x' || currentText.front() == 'X');
            const bool arrayBoundary = previousText.back() == '}' || currentText.front() == '{';
            const bool infinityBoundary = previousText == "Infinity" || currentText == "Infinity";

            // 数値literalとの字句衝突、配列境界、Infinityとの積は '*' を明示する。
            // Infinityは一般symbolと異なり拡張実数sentinelなので、0*Infinity等を暗黙積で曖昧にしない。
            if (exponentMarkerCollision || radixPrefixCollision || arrayBoundary || infinityBoundary)
                output.push_back('*');
            else if (previousIdentifier && currentIdentifier)
                output.push_back(' ');
            else if ((previousIdentifier && (currentDigit || currentText.front() == '('))
                || (previousDigit && currentDigit)
                || currentText.front() == '+' || currentText.front() == '-')
                output.push_back('*');
        }

        output += currentText;
        previousText = std::move(currentText);
    }

    if (parenthesize)
        output.push_back(')');
}

void appendCasesCall(std::string& output, const CallExpr& call, unsigned radix) {
    output += "cases[";
    for (std::size_t i = 0; i < call.arguments.size(); ++i) {
        if (i != 0)
            output += "; ";
        const Expr& branchExpression = call.arguments[i];
        if (!branchExpression.isCall()
            || branchExpression.asCall().head.view() != builtins::names::caseBranch
            || branchExpression.asCall().arguments.empty()
            || branchExpression.asCall().arguments.size() > 2) {
            appendExpr(output, branchExpression, radix, precedenceLowest);
            continue;
        }
        const CallExpr& branch = branchExpression.asCall();
        appendExpr(output, branch.arguments[0], radix, precedenceLowest);
        if (branch.arguments.size() == 2) {
            output += " if ";
            appendExpr(output, branch.arguments[1], radix, precedenceLowest);
        }
    }
    output += "]";
}

void appendGenericCall(std::string& output, const CallExpr& call, unsigned radix) {
    output += call.head.name();
    output.push_back('[');

    for (std::size_t i = 0; i < call.arguments.size(); ++i) {
        if (i != 0)
            output += ", ";

        appendExpr(output, call.arguments[i], radix, precedenceLowest);
    }

    output.push_back(']');
}

[[nodiscard]] std::string_view comparisonSeparator(std::string_view head) noexcept {
    if (head == builtins::names::less)
        return " < ";
    if (head == builtins::names::lessEqual)
        return " <= ";
    if (head == builtins::names::greater)
        return " > ";
    if (head == builtins::names::greaterEqual)
        return " >= ";
    if (head == builtins::names::equal)
        return " == ";
    if (head == builtins::names::notEqual)
        return " != ";
    return {};
}

void appendFunctionSignature(std::string& output, const CallExpr& call, unsigned radix) {
    if (call.arguments.empty() || !call.arguments.front().isSymbol()) {
        appendGenericCall(output, call, radix);
        return;
    }

    output += call.arguments.front().asSymbol().name();
    output.push_back('[');

    for (std::size_t i = 1; i < call.arguments.size(); ++i) {
        if (i != 1)
            output += ", ";
        appendExpr(output, call.arguments[i], radix, precedenceLowest);
    }

    output.push_back(']');
}

void appendCall(
    std::string& output,
    const CallExpr& call,
    unsigned radix,
    int parentPrecedence) {
    const std::string_view head = call.head.view();

    if (head == builtins::names::add && !call.arguments.empty()) {
        appendAddCall(output, call, radix, parentPrecedence);
        return;
    }

    if (head == builtins::names::subtract && call.arguments.size() == 2) {
        std::vector<SignedTerm> terms;
        terms.reserve(2);
        collectSignedTerms(call.arguments[0], false, terms);
        collectSignedTerms(call.arguments[1], true, terms);
        // a-(-b) は表示上は加法になる。ここだけAddと同じ多項式順を適用して、
        // formatter -> parser -> formatter の固定点を保つ。通常の 1-x はそのまま残す。
        if (positiveMagnitudeOfNegative(call.arguments[1]))
            orderUnivariatePolynomialTermsForDisplay(terms);
        appendSignedTerms(output, terms, radix, parentPrecedence);
        return;
    }

    if (head == builtins::names::multiply && !call.arguments.empty()) {
        appendMultiplyCall(output, call, radix, parentPrecedence);
        return;
    }

    if (head == builtins::names::divide && call.arguments.size() == 2) {
        appendBinaryCall(output, call, "/", precedenceMultiplicative, radix, parentPrecedence, true);
        return;
    }

    if (head == builtins::names::power && call.arguments.size() == 2) {
        const bool parenthesize = precedencePower < parentPrecedence;
        if (parenthesize)
            output.push_back('(');

        // ^ は右結合である。左operandがPowerなら括弧を残し，
        // ((a^b)^c) を a^b^c と誤って右結合表示しない。
        appendExpr(output, call.arguments[0], radix, precedencePower + 1);
        output.push_back('^');
        appendExpr(output, call.arguments[1], radix, precedencePower);

        if (parenthesize)
            output.push_back(')');
        return;
    }

    if (head == builtins::names::negate && call.arguments.size() == 1) {
        const bool parenthesize = precedenceUnary < parentPrecedence;
        if (parenthesize)
            output.push_back('(');

        output.push_back('-');
        const Expr& magnitude = call.arguments.front();
        const bool multiplicativeCall = magnitude.isCall()
            && (magnitude.asCall().head.view() == builtins::names::multiply
                || magnitude.asCall().head.view() == builtins::names::divide);
        const bool rationalNumber = magnitude.isNumber()
            && magnitude.asNumber().isReal()
            && !magnitude.asNumber().asReal().isInteger();
        appendExpr(output, magnitude, radix,
            multiplicativeCall || rationalNumber
                ? precedenceMultiplicative
                : precedenceUnary);

        if (parenthesize)
            output.push_back(')');
        return;
    }

    if (head == builtins::names::factorial && call.arguments.size() == 1) {
        const bool parenthesize = precedencePostfix < parentPrecedence;
        if (parenthesize)
            output.push_back('(');

        appendExpr(output, call.arguments.front(), radix, precedencePostfix);
        output.push_back('!');

        if (parenthesize)
            output.push_back(')');
        return;
    }

    if ((head == builtins::names::set || head == builtins::names::setDelayed)
        && call.arguments.size() == 2) {
        appendBinaryCall(output, call, ":=", precedenceAssignment, radix, parentPrecedence, true);
        return;
    }

    const std::string_view separator = comparisonSeparator(head);
    if (!separator.empty() && call.arguments.size() == 2) {
        appendBinaryCall(
            output,
            call,
            separator,
            precedenceComparison,
            radix,
            parentPrecedence,
            true);
        return;
    }

    if (head == builtins::names::cases) {
        appendCasesCall(output, call, radix);
        return;
    }

    if (head == builtins::names::caseBranch) {
        // 内部headが単独で表面化した場合だけgeneric表示する。
        appendGenericCall(output, call, radix);
        return;
    }

    if (head == builtins::names::functionSignature) {
        appendFunctionSignature(output, call, radix);
        return;
    }

    if (head == builtins::names::history
        && call.arguments.size() == 1
        && call.arguments.front().isNumber()
        && call.arguments.front().asNumber().isReal()
        && call.arguments.front().asNumber().asReal().isInteger()) {
        const std::string countText = call.arguments.front().asNumber().toString();
        std::size_t count = 1;
        try {
            count = static_cast<std::size_t>(std::stoull(countText));
        }
        catch (...) {
            appendGenericCall(output, call, radix);
            return;
        }
        output.append(count, '%');
        return;
    }

    if (head == builtins::names::unitApplied
        && call.arguments.size() == 2
        && (call.arguments.back().isString() || call.arguments.back().isSymbol())) {
        appendExpr(output, call.arguments.front(), radix, precedenceMultiplicative);
        output.push_back(' ');
        if (call.arguments.back().isString())
            output += call.arguments.back().asString();
        else
            output += call.arguments.back().asSymbol().name();
        return;
    }

    appendGenericCall(output, call, radix);
}

void appendExpr(
    std::string& output,
    const Expr& expression,
    unsigned radix,
    int parentPrecedence) {
    using expression::ExprKind;

    switch (expression.kind()) {
    case ExprKind::Number:
        appendNumber(output, expression.asNumber(), radix, parentPrecedence);
        break;

    case ExprKind::DecimalApproximation:
        output += expression.asDecimalApproximation().text();
        break;

    case ExprKind::ComplexDecimalApproximation:
        output += expression.asComplexDecimalApproximation().text();
        break;

    case ExprKind::Boolean:
        output += expression.asBoolean() ? "True" : "False";
        break;

    case ExprKind::String:
        appendEscapedString(output, expression.asString());
        break;

    case ExprKind::Symbol:
        output += expression.asSymbol().name();
        break;

    case ExprKind::Array: {
        const auto& array = expression.asArray();
        if (expression::braceLiteralPreservesShape(array.shape)) {
            appendArrayDimension(output, array, 0, 0, radix);
            break;
        }

        output += builtins::names::reshape;
        output += "[{}, ";
        appendShapeLiteral(output, array.shape);
        output.push_back(']');
        break;
    }

    case ExprKind::List: {
        output.push_back('{');
        const auto& list = expression.asList();
        for (std::size_t i = 0; i < list.elements.size(); ++i) {
            if (i != 0)
                output += ", ";
            appendExpr(output, list.elements[i], radix, precedenceLowest);
        }
        output.push_back('}');
        break;
    }

    case ExprKind::Call:
        appendCall(output, expression.asCall(), radix, parentPrecedence);
        break;

    case ExprKind::SolutionSet:
        appendSolutionSet(output, expression.asSolutionSet(), radix);
        break;
    }
}

} // namespace

std::string formatExpr(const expression::Expr& expression, unsigned radix) {
    std::string output;
    appendExpr(output, expression, radix, precedenceLowest);
    return output;
}

std::string formatAssumptions(
    const mathematics::AssumptionSet& assumptions,
    unsigned radix) {
    std::string output;
    appendConditions(output, assumptions, radix);
    return output;
}

std::string trimRedundantFractionalZeros(std::string_view text) {
    std::string output;
    output.reserve(text.size());

    bool inString = false;
    bool escaped = false;
    for (std::size_t i = 0; i < text.size();) {
        const char character = text[i];
        if (inString) {
            output.push_back(character);
            ++i;
            if (escaped) {
                escaped = false;
                continue;
            }
            if (character == '\\') {
                escaped = true;
                continue;
            }
            if (character == '"')
                inString = false;
            continue;
        }

        if (character == '"') {
            inString = true;
            output.push_back(character);
            ++i;
            continue;
        }

        const unsigned char byte = static_cast<unsigned char>(character);
        if (std::isdigit(byte) == 0) {
            output.push_back(character);
            ++i;
            continue;
        }

        std::size_t integerEnd = i;
        while (integerEnd < text.size()
            && std::isdigit(static_cast<unsigned char>(text[integerEnd])) != 0)
            ++integerEnd;

        if (integerEnd >= text.size() || text[integerEnd] != '.'
            || integerEnd + 1 >= text.size()
            || std::isdigit(static_cast<unsigned char>(text[integerEnd + 1])) == 0) {
            output.append(text.substr(i, integerEnd - i));
            i = integerEnd;
            continue;
        }

        std::size_t fractionEnd = integerEnd + 1;
        while (fractionEnd < text.size()
            && std::isdigit(static_cast<unsigned char>(text[fractionEnd])) != 0)
            ++fractionEnd;

        std::size_t trimmedEnd = fractionEnd;
        while (trimmedEnd > integerEnd + 1 && text[trimmedEnd - 1] == '0')
            --trimmedEnd;

        output.append(text.substr(i, integerEnd - i));
        if (trimmedEnd != integerEnd + 1) {
            output.push_back('.');
            output.append(text.substr(integerEnd + 1, trimmedEnd - integerEnd - 1));
        }
        i = fractionEnd;
    }

    return output;
}

} // namespace mmcal::formatting
