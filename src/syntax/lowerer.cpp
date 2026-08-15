// 構文木からExprへのlowering
#include "lowerer.hpp"

#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "expression/array_utils.hpp"
#include "numeric/integer_algorithms.hpp"

#include <charconv>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

namespace mmcal::syntax {
namespace {

using expression::Expr;
using expression::Symbol;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;
using numeric::RealNumber;

[[nodiscard]] Expr integerExpr(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] Expr callExpr(
    symbols::SymbolTable& symbolTable,
    std::string_view name,
    std::vector<Expr> arguments) {
    return Expr::call(symbolTable.intern(name), std::move(arguments));
}

[[nodiscard]] unsigned parseRadix(
    std::string_view text,
    source::SourceSpan span) {
    unsigned radix = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), radix);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size()
        || radix < 2 || radix > 36)
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "Radix must be an integer from 2 through 36",
            span);

    return radix;
}

[[nodiscard]] std::uint64_t parseUnsignedExponent(
    std::string_view text,
    source::SourceSpan span) {
    std::uint64_t exponent = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), exponent);
    if (result.ec == std::errc::result_out_of_range)
        error::throwCalcError(
            error::CalcErrorType::Overflow,
            "Decimal exponent is too large",
            span);
    if (result.ec != std::errc{} || result.ptr != text.data() + text.size())
        error::throwCalcError(
            error::CalcErrorType::Syntax,
            "Invalid decimal exponent",
            span);

    return exponent;
}

[[nodiscard]] Number exactDecimal(std::string_view text, source::SourceSpan span) {
    const std::size_t exponentMark = text.find_first_of("eE");
    const std::string_view mantissa = text.substr(0, exponentMark);

    bool negativeExponent = false;
    std::uint64_t exponent = 0;
    if (exponentMark != std::string_view::npos) {
        std::string_view exponentText = text.substr(exponentMark + 1);
        if (!exponentText.empty() && (exponentText.front() == '+' || exponentText.front() == '-')) {
            negativeExponent = exponentText.front() == '-';
            exponentText.remove_prefix(1);
        }
        exponent = parseUnsignedExponent(exponentText, span);
    }

    const std::size_t point = mantissa.find('.');
    const std::size_t fractionDigits = point == std::string_view::npos
        ? 0
        : mantissa.size() - point - 1;

    std::string digits;
    digits.reserve(mantissa.size());
    for (const char value : mantissa)
        if (value != '.')
            digits.push_back(value);

    BigInt numerator = BigInt::parse(digits.empty() ? "0" : digits);
    std::uint64_t denominatorPower = static_cast<std::uint64_t>(fractionDigits);

    if (negativeExponent) {
        if (denominatorPower > std::numeric_limits<std::uint64_t>::max() - exponent)
            error::throwCalcError(
                error::CalcErrorType::Overflow,
                "Decimal scale is too large",
                span);
        denominatorPower += exponent;
    }
    else if (exponent >= denominatorPower) {
        numerator *= numeric::pow(BigInt{10}, exponent - denominatorPower);
        denominatorPower = 0;
    }
    else
        denominatorPower -= exponent;

    if (denominatorPower == 0)
        return Number{std::move(numerator)};

    const BigInt denominator = numeric::pow(BigInt{10}, denominatorPower);
    return Number{Rational{std::move(numerator), denominator}};
}

[[nodiscard]] Number parseNumberLiteral(
    const NumberLiteralSyntax& number,
    source::SourceSpan span) {
    const std::string_view text = number.text;

    try {
        const std::size_t radixMark = text.find('#');
        if (radixMark != std::string_view::npos) {
            const unsigned radix = parseRadix(text.substr(0, radixMark), span);
            const std::string_view digits = text.substr(radixMark + 1);
            if (digits.find('.') == std::string_view::npos)
                return Number{BigInt::parse(digits, radix)};
            return Number{Rational::parse(digits, radix)};
        }

        if (text.size() > 2 && text[0] == '0') {
            unsigned radix = 0;
            switch (text[1]) {
            case 'b':
            case 'B': radix = 2; break;
            case 'o':
            case 'O': radix = 8; break;
            case 'x':
            case 'X': radix = 16; break;
            default: break;
            }
            if (radix != 0)
                return Number{BigInt::parse(text.substr(2), radix)};
        }
        return exactDecimal(text, span);
    }
    catch (const error::CalcError&) {
        throw;
    }
    catch (const std::overflow_error& exception) {
        error::throwCalcError(error::CalcErrorType::Overflow, exception.what(), span);
    }
    catch (const std::exception& exception) {
        error::throwCalcError(error::CalcErrorType::Syntax, exception.what(), span);
    }
}

// 構文上矩形なbraceだけを先に判定する。leafの意味評価は行わない。
[[nodiscard]] bool rectangularArrayShape(
    const ArrayLiteralSyntax& array,
    std::vector<std::size_t>& shape) {
    shape = {array.elements.size()};
    if (array.elements.empty())
        return true;

    const auto* firstChild = std::get_if<ArrayLiteralSyntax>(&array.elements.front()->data);
    if (!firstChild) {
        for (std::size_t i = 1; i < array.elements.size(); ++i)
            if (std::holds_alternative<ArrayLiteralSyntax>(array.elements[i]->data))
                return false;
        return true;
    }

    std::vector<std::size_t> childShape;
    if (!rectangularArrayShape(*firstChild, childShape))
        return false;
    for (std::size_t i = 1; i < array.elements.size(); ++i) {
        const auto* child = std::get_if<ArrayLiteralSyntax>(&array.elements[i]->data);
        if (!child)
            return false;
        std::vector<std::size_t> currentShape;
        if (!rectangularArrayShape(*child, currentShape) || currentShape != childShape)
            return false;
    }
    shape.insert(shape.end(), childShape.begin(), childShape.end());
    return true;
}

[[nodiscard]] std::string decodeString(
    std::string_view text,
    source::SourceSpan span) {
    std::string result;
    result.reserve(text.size() >= 2 ? text.size() - 2 : 0);

    for (std::size_t i = 1; i + 1 < text.size(); ++i) {
        const char value = text[i];
        if (value != '\\') {
            result.push_back(value);
            continue;
        }

        if (++i + 1 > text.size())
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Incomplete escape sequence",
                span);

        switch (text[i]) {
        case '\\':
            result.push_back('\\');
            break;
        case '"':
            result.push_back('"');
            break;
        case 'n':
            result.push_back('\n');
            break;
        case 'r':
            result.push_back('\r');
            break;
        case 't':
            result.push_back('\t');
            break;
        default:
            error::throwCalcError(
                error::CalcErrorType::Syntax,
                "Unsupported string escape sequence",
                span);
        }
    }

    return result;
}

[[nodiscard]] std::string_view binaryName(BinaryOperator operation) {
    switch (operation) {
    case BinaryOperator::Add:
        return builtins::names::add;
    case BinaryOperator::Subtract:
        return builtins::names::subtract;
    case BinaryOperator::Multiply:
    case BinaryOperator::ImplicitMultiply:
        return builtins::names::multiply;
    case BinaryOperator::Divide:
        return builtins::names::divide;
    case BinaryOperator::Power:
        return builtins::names::power;
    }

    return "UnknownBinary";
}

[[nodiscard]] std::string_view comparisonName(ComparisonOperator operation) {
    switch (operation) {
    case ComparisonOperator::Less:
        return builtins::names::less;
    case ComparisonOperator::LessEqual:
        return builtins::names::lessEqual;
    case ComparisonOperator::Greater:
        return builtins::names::greater;
    case ComparisonOperator::GreaterEqual:
        return builtins::names::greaterEqual;
    case ComparisonOperator::Equal:
        return builtins::names::equal;
    case ComparisonOperator::NotEqual:
        return builtins::names::notEqual;
    }

    return "UnknownComparison";
}


} // namespace

LoweringOptions LoweringOptions::defaults() {
    return LoweringOptions{
        symbols::defaultSymbolRegistry().sourcePredefinedNames(),
        {"sqrt"},
        {}
    };
}

bool LoweringOptions::isConstant(std::string_view name) const {
    return constants.contains(std::string{name});
}

bool LoweringOptions::isFunction(std::string_view name) const {
    return functions.contains(std::string{name});
}

bool LoweringOptions::isVariable(std::string_view name) const {
    return variables.contains(std::string{name});
}

Lowerer::Lowerer(
    LoweringOptions options,
    const symbols::SymbolRegistry& symbolRegistry,
    symbols::SymbolTable& symbolTable)
    : options_(std::move(options)),
      symbolRegistry_(symbolRegistry),
      symbolTable_(symbolTable) {}

Expr Lowerer::lower(const SyntaxTree& tree) const {
    return lowerNode(tree.root(), nullptr);
}

LoweringResult Lowerer::lowerTracked(const SyntaxTree& tree) const {
    return lowerTracked(tree, nullptr);
}

LoweringResult Lowerer::lowerTracked(
    const SyntaxTree& tree,
    std::shared_ptr<const source::SourceDocument> document) const {
    expression::OriginMap origins{std::move(document)};
    Expr expression = lowerNode(tree.root(), &origins);
    return LoweringResult{std::move(expression), std::move(origins)};
}

const LoweringOptions& Lowerer::options() const noexcept {
    return options_;
}

LoweringOptions& Lowerer::options() noexcept {
    return options_;
}

Expr Lowerer::lowerNode(
    const SyntaxNode& node,
    expression::OriginMap* origins) const {
    Expr result = [&]() -> Expr {
        if (const auto* number = std::get_if<NumberLiteralSyntax>(&node.data))
            return lowerNumber(*number, node.span);
        if (const auto* string = std::get_if<StringLiteralSyntax>(&node.data))
            return lowerString(*string, node.span);
        if (const auto* identifier = std::get_if<IdentifierSyntax>(&node.data))
            return lowerIdentifier(*identifier);
        if (const auto* history = std::get_if<HistoryReferenceSyntax>(&node.data)) {
            if (history->kind == HistoryReferenceKind::Input)
                return callExpr(
                    symbolTable_,
                    builtins::names::inputHistory,
                    {integerExpr(-static_cast<std::int64_t>(history->depth))});
            return callExpr(
                symbolTable_,
                builtins::names::history,
                {integerExpr(static_cast<std::int64_t>(history->depth))});
        }
        if (const auto* array = std::get_if<ArrayLiteralSyntax>(&node.data))
            return lowerArray(*array, node.span, origins);
        if (const auto* call = std::get_if<CallSyntax>(&node.data))
            return lowerCall(*call, origins);
        if (const auto* group = std::get_if<GroupSyntax>(&node.data))
            return lowerNode(*group->expression, origins);
        if (const auto* unary = std::get_if<UnarySyntax>(&node.data)) {
            Expr operand = lowerNode(*unary->operand, origins);
            if (unary->operation == UnaryOperator::Plus)
                return operand;
            return callExpr(symbolTable_, builtins::names::negate, {std::move(operand)});
        }
        if (const auto* binary = std::get_if<BinarySyntax>(&node.data))
            return callExpr(
                symbolTable_,
                binaryName(binary->operation),
                {lowerNode(*binary->left, origins), lowerNode(*binary->right, origins)});
        if (const auto* postfix = std::get_if<PostfixSyntax>(&node.data))
            return callExpr(
                symbolTable_,
                builtins::names::factorial,
                {lowerNode(*postfix->operand, origins)});
        if (const auto* comparison = std::get_if<ComparisonSyntax>(&node.data))
            return lowerComparison(*comparison, origins);
        if (const auto* assignment = std::get_if<AssignmentSyntax>(&node.data))
            return lowerAssignment(*assignment, origins);
        if (const auto* signature = std::get_if<FunctionSignatureSyntax>(&node.data)) {
            std::vector<Expr> arguments;
            arguments.reserve(signature->parameters.size() + 1);
            arguments.emplace_back(symbolTable_.intern(signature->name));
            for (const std::string& parameter : signature->parameters)
                arguments.emplace_back(symbolTable_.intern(parameter));
            return callExpr(symbolTable_, builtins::names::functionSignature, std::move(arguments));
        }
        if (const auto* unit = std::get_if<UnitAppliedSyntax>(&node.data))
            return callExpr(
                symbolTable_,
                builtins::names::unitApplied,
                {lowerNode(*unit->value, origins), Expr{unit->unit}});

        error::throwCalcError(
            error::CalcErrorType::Internal,
            "Unknown syntax node kind",
            node.span);
    }();

    if (origins)
        origins->record(result, node.span);

    return result;
}

Expr Lowerer::lowerNumber(
    const NumberLiteralSyntax& number,
    source::SourceSpan span) const {
    return Expr{parseNumberLiteral(number, span)};
}

Expr Lowerer::lowerString(
    const StringLiteralSyntax& string,
    source::SourceSpan span) const {
    return Expr{decodeString(string.text, span)};
}

Expr Lowerer::lowerIdentifier(const IdentifierSyntax& identifier) const {
    const Symbol symbol = symbolTable_.intern(identifier.name);
    const symbols::PredefinedSymbolDefinition* definition = symbolRegistry_.find(identifier.name);
    if (!definition)
        return Expr{symbol};

    switch (definition->kind) {
    case symbols::PredefinedSymbolKind::ImaginaryUnit:
        return Expr{Number::complex(RealNumber{}, RealNumber{BigInt{1}})};
    case symbols::PredefinedSymbolKind::BooleanTrue:
        return Expr{true};
    case symbols::PredefinedSymbolKind::BooleanFalse:
        return Expr{false};
    case symbols::PredefinedSymbolKind::SymbolicConstant:
    case symbols::PredefinedSymbolKind::MathematicalDomain:
    case symbols::PredefinedSymbolKind::EnumeratedValue:
        return Expr{symbol};
    }

    error::throwCalcError(
        error::CalcErrorType::Internal,
        "Unknown predefined symbol kind");
}

Expr Lowerer::lowerArray(
    const ArrayLiteralSyntax& array,
    source::SourceSpan span,
    expression::OriginMap* origins) const {
    try {
        std::vector<std::size_t> shape;
        if (rectangularArrayShape(array, shape)) {
            expression::ArrayBuilder builder;
            builder.reserve(expression::arrayElementCount(shape));

            const auto appendLeaves = [&](auto&& self, const ArrayLiteralSyntax& current) -> void {
                for (const SyntaxNodePtr& element : current.elements) {
                    if (const auto* child = std::get_if<ArrayLiteralSyntax>(&element->data)) {
                        self(self, *child);
                        continue;
                    }
                    // 数値literalはExpr nodeを経由せずpacked pageへ直接入れる。
                    if (const auto* number = std::get_if<NumberLiteralSyntax>(&element->data)) {
                        builder.append(parseNumberLiteral(*number, element->span));
                        continue;
                    }
                    builder.append(lowerNode(*element, origins));
                }
            };
            appendLeaves(appendLeaves, array);
            return Expr::array(builder.finish(std::move(shape)));
        }

        std::vector<Expr> elements;
        elements.reserve(array.elements.size());
        for (const SyntaxNodePtr& element : array.elements)
            elements.push_back(lowerNode(*element, origins));

        // ragged / scalar-Array混在は従来どおり一般brace Listとして保持する。
        return expression::braceValue(std::move(elements));
    }
    catch (const std::length_error& exception) {
        error::throwCalcError(error::CalcErrorType::Overflow, exception.what(), span);
    }
}

Expr Lowerer::lowerCall(
    const CallSyntax& call,
    expression::OriginMap* origins) const {
    std::vector<Expr> arguments;
    arguments.reserve(call.arguments.size());
    for (const SyntaxNodePtr& argument : call.arguments) {
        // UnDefはHoldAllなので、Iのように通常loweringで値へ変換されるpredefined名も
        // 名前そのものとして渡す。これにより全protected symbolを同じ規則で拒否できる。
        if (call.name == builtins::names::undefine) {
            if (const auto* identifier = std::get_if<IdentifierSyntax>(&argument->data)) {
                Expr lowered{symbolTable_.intern(identifier->name)};
                if (origins)
                    origins->record(lowered, argument->span);
                arguments.push_back(std::move(lowered));
                continue;
            }
        }
        arguments.push_back(lowerNode(*argument, origins));
    }

    return Expr::call(symbolTable_.intern(call.name), std::move(arguments));
}

Expr Lowerer::lowerComparison(
    const ComparisonSyntax& comparison,
    expression::OriginMap* origins) const {
    std::vector<Expr> comparisons;
    comparisons.reserve(comparison.operations.size());

    for (std::size_t i = 0; i < comparison.operations.size(); ++i)
        comparisons.push_back(callExpr(
            symbolTable_,
            comparisonName(comparison.operations[i]),
            {
                lowerNode(*comparison.operands[i], origins),
                lowerNode(*comparison.operands[i + 1], origins)
            }));

    if (comparisons.size() == 1)
        return comparisons.front();

    return callExpr(symbolTable_, builtins::names::logicalAnd, std::move(comparisons));
}

Expr Lowerer::lowerAssignment(
    const AssignmentSyntax& assignment,
    expression::OriginMap* origins) const {
    Expr target = lowerNode(*assignment.target, origins);
    Expr value = lowerNode(*assignment.value, origins);

    if (std::holds_alternative<FunctionSignatureSyntax>(assignment.target->data))
        return callExpr(symbolTable_, builtins::names::setDelayed, {std::move(target), std::move(value)});

    return callExpr(symbolTable_, builtins::names::set, {std::move(target), std::move(value)});
}

} // namespace mmcal::syntax
