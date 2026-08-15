// 評価済み値の軽量introspection
#include "explain.hpp"

#include "error/error_message.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/number.hpp"

#include <algorithm>
#include <cstddef>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace mmcal::builtins {
namespace {

using expression::ArrayExpr;
using expression::ArrayStorageKind;
using expression::Expr;
using expression::ExprKind;
using numeric::ApproximationOrigin;
using numeric::BigInt;
using numeric::ComplexDecimalApproximation;
using numeric::DecimalApproximation;
using numeric::Number;
using mathematics::ArithmeticClass;
using mathematics::ConstantDefinition;
using symbols::PredefinedSymbolDefinition;
using symbols::PredefinedSymbolId;
using symbols::PredefinedSymbolKind;

[[nodiscard]] Expr integer(std::size_t value) {
    return Expr{Number{BigInt::fromUnsigned(value)}};
}

[[nodiscard]] Expr text(std::string_view value) {
    return Expr{std::string{value}};
}

[[nodiscard]] Expr pair(std::string_view name, Expr value) {
    return Expr::list({text(name), std::move(value)});
}

void add(std::vector<Expr>& properties, std::string_view name, Expr value) {
    properties.push_back(pair(name, std::move(value)));
}

[[nodiscard]] Expr brace(std::vector<Expr> elements) {
    return Expr::list(std::move(elements));
}

[[nodiscard]] std::string_view kindName(ExprKind kind) noexcept {
    switch (kind) {
    case ExprKind::Number: return "Number";
    case ExprKind::DecimalApproximation: return "DecimalApproximation";
    case ExprKind::ComplexDecimalApproximation: return "ComplexDecimalApproximation";
    case ExprKind::Boolean: return "Boolean";
    case ExprKind::String: return "String";
    case ExprKind::Symbol: return "Symbol";
    case ExprKind::Array: return "Array";
    case ExprKind::List: return "List";
    case ExprKind::Call: return "Call";
    case ExprKind::SolutionSet: return "SolutionSet";
    }
    return "Unknown";
}

[[nodiscard]] std::string_view arithmeticClassName(ArithmeticClass value) noexcept {
    switch (value) {
    case ArithmeticClass::Unknown: return "Unknown";
    case ArithmeticClass::Algebraic: return "Algebraic";
    case ArithmeticClass::Transcendental: return "Transcendental";
    }
    return "Unknown";
}

[[nodiscard]] std::string_view semanticKind(
    const Expr& value,
    const symbols::SymbolRegistry& symbols,
    const mathematics::MathRegistry& mathematics) noexcept {
    if (!value.isSymbol())
        return kindName(value.kind());

    const auto& symbol = value.asSymbol();
    if (mathematics.findConstant(symbol))
        return "Constant";

    const PredefinedSymbolDefinition* definition = symbols.find(symbol);
    if (!definition)
        return "Symbol";

    switch (definition->kind) {
    case PredefinedSymbolKind::SymbolicConstant:
    case PredefinedSymbolKind::ImaginaryUnit:
        return "Constant";
    case PredefinedSymbolKind::MathematicalDomain:
        return "MathematicalDomain";
    case PredefinedSymbolKind::EnumeratedValue:
        switch (definition->id) {
        case PredefinedSymbolId::DegreeUnit:
        case PredefinedSymbolId::RadianUnit:
        case PredefinedSymbolId::GradianUnit:
            return "AngleUnit";
        default:
            return "EnumeratedValue";
        }
    case PredefinedSymbolKind::BooleanTrue:
    case PredefinedSymbolKind::BooleanFalse:
        return "Boolean";
    }
    return "Symbol";
}

void addPredefinedInternal(
    const PredefinedSymbolDefinition* definition,
    std::vector<Expr>& properties,
    bool internal) {
    if (!internal || !definition)
        return;
    add(properties, "Predefined", Expr{true});
    add(properties, "Protected", Expr{definition->protectedName});
}

void explainMathematicalConstant(
    const ConstantDefinition& definition,
    const PredefinedSymbolDefinition* predefined,
    std::vector<Expr>& properties,
    bool internal) {
    const auto& p = definition.properties;
    add(properties, "Domain", text(p.real ? "Real" : "Complex"));
    add(properties, "Exactness", text(p.exact ? "Exact" : "Unknown"));
    add(properties, "Name", text(definition.symbol.view()));
    add(properties, "Real", Expr{p.real});
    if (p.positive)
        add(properties, "Positive", Expr{true});
    if (p.irrational)
        add(properties, "Irrational", Expr{true});
    if (p.arithmeticClass != ArithmeticClass::Unknown)
        add(properties, "ArithmeticClass", text(arithmeticClassName(p.arithmeticClass)));
    addPredefinedInternal(predefined, properties, internal);
}

void explainPredefinedSymbol(
    const PredefinedSymbolDefinition& definition,
    std::vector<Expr>& properties,
    bool internal) {
    switch (definition.id) {
    case PredefinedSymbolId::I:
        add(properties, "Domain", text("Complex"));
        add(properties, "Exactness", text("Exact"));
        add(properties, "Name", text(definition.symbol.view()));
        add(properties, "Real", Expr{false});
        add(properties, "Zero", Expr{false});
        add(properties, "ImaginaryUnit", Expr{true});
        add(properties, "ArithmeticClass", text("Algebraic"));
        break;
    case PredefinedSymbolId::Infinity:
        add(properties, "Domain", text("ExtendedReal"));
        add(properties, "Exactness", text("Exact"));
        add(properties, "Name", text(definition.symbol.view()));
        add(properties, "Infinite", Expr{true});
        add(properties, "Finite", Expr{false});
        add(properties, "Sign", text("Positive"));
        break;
    case PredefinedSymbolId::IntegerDomain:
    case PredefinedSymbolId::RationalDomain:
    case PredefinedSymbolId::RealDomain:
    case PredefinedSymbolId::ComplexDomain:
        add(properties, "Exactness", text("Exact"));
        add(properties, "Name", text(definition.symbol.view()));
        break;
    case PredefinedSymbolId::DegreeUnit:
    case PredefinedSymbolId::RadianUnit:
    case PredefinedSymbolId::GradianUnit:
        add(properties, "Domain", text("AngleUnit"));
        add(properties, "Exactness", text("Exact"));
        add(properties, "Name", text(definition.symbol.view()));
        break;
    case PredefinedSymbolId::Pi:
    case PredefinedSymbolId::E:
    case PredefinedSymbolId::Phi:
        // MathRegistry側で説明する。
        add(properties, "Domain", text("Unknown"));
        add(properties, "Exactness", text("Exact"));
        add(properties, "Name", text(definition.symbol.view()));
        break;
    case PredefinedSymbolId::True:
    case PredefinedSymbolId::False:
        add(properties, "Domain", text("Boolean"));
        add(properties, "Exactness", text("Exact"));
        add(properties, "Name", text(definition.symbol.view()));
        break;
    }

    addPredefinedInternal(&definition, properties, internal);
}

[[nodiscard]] std::string_view exactNumberDomain(const Number& number) noexcept {
    if (number.isComplex())
        return "Complex";
    const auto& real = number.asReal();
    return real.isInteger() ? "Integer" : "Rational";
}

[[nodiscard]] std::string_view arrayDomain(ArrayStorageKind kind) noexcept {
    switch (kind) {
    case ArrayStorageKind::Integer: return "Integer";
    case ArrayStorageKind::Rational: return "Rational";
    case ArrayStorageKind::Number: return "Number";
    case ArrayStorageKind::DecimalApproximation: return "Real";
    case ArrayStorageKind::ComplexDecimalApproximation: return "Complex";
    case ArrayStorageKind::Generic: return "Expression";
    }
    return "Unknown";
}

[[nodiscard]] std::string_view arrayExactness(ArrayStorageKind kind) noexcept {
    switch (kind) {
    case ArrayStorageKind::Integer:
    case ArrayStorageKind::Rational:
    case ArrayStorageKind::Number:
        return "Exact";
    case ArrayStorageKind::DecimalApproximation:
    case ArrayStorageKind::ComplexDecimalApproximation:
        return "CertifiedApproximation";
    case ArrayStorageKind::Generic:
        return "Unknown";
    }
    return "Unknown";
}

[[nodiscard]] std::string_view storageName(ArrayStorageKind kind) noexcept {
    switch (kind) {
    case ArrayStorageKind::Integer: return "Integer";
    case ArrayStorageKind::Rational: return "Rational";
    case ArrayStorageKind::Number: return "Number";
    case ArrayStorageKind::DecimalApproximation: return "DecimalApproximation";
    case ArrayStorageKind::ComplexDecimalApproximation: return "ComplexDecimalApproximation";
    case ArrayStorageKind::Generic: return "Generic";
    }
    return "Unknown";
}

[[nodiscard]] Expr enclosure(const DecimalApproximation& value) {
    return brace({
        Expr{Number{value.certifiedLower()}},
        Expr{Number{value.certifiedUpper()}}
    });
}

[[nodiscard]] Expr enclosure(const ComplexDecimalApproximation& value) {
    return brace({
        pair("Real", enclosure(value.real())),
        pair("Imaginary", enclosure(value.imaginary()))
    });
}

void explainNumber(const Number& number, std::vector<Expr>& properties, bool internal) {
    add(properties, "Domain", text(exactNumberDomain(number)));
    add(properties, "Exactness", text("Exact"));
    add(properties, "Zero", Expr{number.isZero()});

    if (number.isReal()) {
        const auto& real = number.asReal();
        add(properties, "Sign", text(real.isZero() ? "Zero" : real.isNegative() ? "Negative" : "Positive"));
        if (real.isInteger())
            add(properties, "BitLength", integer(real.asInteger().bitLength()));
        else {
            add(properties, "NumeratorBitLength", integer(real.asRational().numerator().bitLength()));
            add(properties, "DenominatorBitLength", integer(real.asRational().denominator().bitLength()));
        }
    }

    if (internal)
        add(properties, "Representation", text(number.isReal() ? "Number<RealNumber>" : "Number<ComplexNumber>"));
}

void explainDecimal(
    const DecimalApproximation& value,
    std::vector<Expr>& properties,
    bool internal) {
    add(properties, "Domain", text("Real"));
    add(properties, "Exactness", text("CertifiedApproximation"));
    add(properties, "RequestedFractionalDigits", integer(value.requestedFractionalDigits()));
    add(properties, "DisplayedFractionalDigits", integer(value.fractionalDigits()));
    add(properties, "Rounded", Expr{value.isRounded()});
    add(properties, "Enclosure", enclosure(value));

    if (internal) {
        add(properties, "Representation", text("DecimalApproximation"));
        add(properties, "ApproximationOrigin", text(
            value.origin() == ApproximationOrigin::ExactValue
                ? "ExactValue"
                : "CertifiedInterval"));
        add(properties, "EnclosureIsPoint", Expr{value.enclosureIsPoint()});
    }
}

void explainComplexDecimal(
    const ComplexDecimalApproximation& value,
    std::vector<Expr>& properties,
    bool internal) {
    add(properties, "Domain", text("Complex"));
    add(properties, "Exactness", text("CertifiedApproximation"));
    add(properties, "RequestedFractionalDigits", integer(std::min(
        value.real().requestedFractionalDigits(),
        value.imaginary().requestedFractionalDigits())));
    add(properties, "Enclosure", enclosure(value));

    if (internal) {
        add(properties, "Representation", text("ComplexDecimalApproximation"));
        add(properties, "RealExactlyZero", Expr{value.realExactlyZero()});
        add(properties, "ImaginaryExactlyZero", Expr{value.imaginaryExactlyZero()});
    }
}

void explainArray(const ArrayExpr& array, std::vector<Expr>& properties, bool internal) {
    add(properties, "Domain", text(arrayDomain(array.storageKind())));
    add(properties, "Exactness", text(arrayExactness(array.storageKind())));
    add(properties, "ArrayRank", integer(array.rank()));

    std::vector<BigInt> dimensions;
    dimensions.reserve(array.shape.size());
    for (const std::size_t dimension : array.shape)
        dimensions.push_back(BigInt::fromUnsigned(dimension));
    const std::size_t dimensionCount = dimensions.size();
    add(properties, "Dimensions", Expr::integerArray({dimensionCount}, std::move(dimensions)));

    add(properties, "ElementCount", integer(array.size()));
    add(properties, "Rectangular", Expr{true});
    add(properties, "Empty", Expr{array.empty()});

    if (array.isVector())
        add(properties, "Vector", Expr{true});
    if (array.isMatrix()) {
        add(properties, "Matrix", Expr{true});
        const bool square = array.shape[0] == array.shape[1];
        add(properties, "Square", Expr{square});
        if (square)
            add(properties, "Order", integer(array.shape[0]));
    }

    if (internal) {
        add(properties, "Representation", text("ArrayExpr"));
        add(properties, "Storage", text(storageName(array.storageKind())));
        add(properties, "Contiguous", Expr{array.isContiguous()});
        add(properties, "StoredExpressions", Expr{array.hasStoredExpressions()});
    }
}

[[nodiscard]] bool internalMode(std::span<const Expr> arguments) {
    if (arguments.size() == 1)
        return false;
    if (!arguments[1].isString())
        error::throwCalcError(error::CalcErrorType::Type,
            "explain mode must be the string \"internal\"");
    if (arguments[1].asString() != "internal")
        error::throwCalcError(error::CalcErrorType::Domain,
            "Unknown explain mode: " + arguments[1].asString());
    return true;
}

} // namespace

expression::Expr evaluateExplain(
    std::span<const expression::Expr> arguments,
    const symbols::SymbolRegistry& symbols,
    const mathematics::MathRegistry& mathematics) {
    if (arguments.empty() || arguments.size() > 2)
        error::throwCalcError(error::CalcErrorType::Type,
            "explain expects a value and optional mode");

    const bool internal = internalMode(arguments);
    const Expr& value = arguments.front();
    std::vector<Expr> properties;
    properties.reserve(16);
    add(properties, "Kind", text(semanticKind(value, symbols, mathematics)));

    switch (value.kind()) {
    case ExprKind::Number:
        explainNumber(value.asNumber(), properties, internal);
        break;
    case ExprKind::DecimalApproximation:
        explainDecimal(value.asDecimalApproximation(), properties, internal);
        break;
    case ExprKind::ComplexDecimalApproximation:
        explainComplexDecimal(value.asComplexDecimalApproximation(), properties, internal);
        break;
    case ExprKind::Boolean:
        add(properties, "Domain", text("Boolean"));
        add(properties, "Exactness", text("Exact"));
        break;
    case ExprKind::String:
        add(properties, "Domain", text("String"));
        add(properties, "Exactness", text("Exact"));
        break;
    case ExprKind::Symbol: {
        const auto& symbol = value.asSymbol();
        const PredefinedSymbolDefinition* predefined = symbols.find(symbol);
        if (const ConstantDefinition* constant = mathematics.findConstant(symbol))
            explainMathematicalConstant(*constant, predefined, properties, internal);
        else if (predefined)
            explainPredefinedSymbol(*predefined, properties, internal);
        else {
            add(properties, "Domain", text("Unknown"));
            add(properties, "Exactness", text("Exact"));
            add(properties, "Name", text(symbol.view()));
        }
        break;
    }
    case ExprKind::Array:
        explainArray(value.asArray(), properties, internal);
        break;
    case ExprKind::List:
        add(properties, "Domain", text("Expression"));
        add(properties, "Exactness", text("Unknown"));
        add(properties, "ElementCount", integer(value.asList().size()));
        add(properties, "Rectangular", Expr{false});
        if (internal)
            add(properties, "Representation", text("ListExpr"));
        break;
    case ExprKind::Call:
        add(properties, "Domain", text("Unknown"));
        add(properties, "Exactness", text("Unknown"));
        add(properties, "Head", text(value.asCall().head.view()));
        add(properties, "ArgumentCount", integer(value.asCall().arguments.size()));
        if (internal)
            add(properties, "Representation", text("CallExpr"));
        break;
    case ExprKind::SolutionSet:
        add(properties, "Domain", text("SolutionSet"));
        add(properties, "Exactness", text("Unknown"));
        if (internal)
            add(properties, "Representation", text("SolutionSet"));
        break;
    }

    return Expr::list(std::move(properties));
}

} // namespace mmcal::builtins
