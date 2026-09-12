// Plot sampling専用IRとBigFloat executor
#include "plot_program.hpp"

#include "approximation/certification_error.hpp"
#include "approximation/certified_atan.hpp"
#include "approximation/certified_constants.hpp"
#include "approximation/certified_elementary_functions.hpp"
#include "approximation/certified_exponential.hpp"
#include "approximation/certified_logarithm.hpp"
#include "approximation/certified_sqrt.hpp"
#include "approximation/certified_special_functions.hpp"
#include "approximation/certified_trigonometry.hpp"
#include "approximation/real_interval.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rounding_mode.hpp"

#include <algorithm>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace mmcal::plot {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using mathematics::ConstantId;
using numeric::BigFloat;
using numeric::BigInt;
using numeric::Rational;
using numeric::RoundingMode;

inline constexpr std::size_t maxPlotInstructions = 4096;
inline constexpr std::size_t maxPlotConstants = 1024;
inline constexpr std::size_t elementaryGuardBits = 16;

struct CompileContext final {
    const expression::Symbol& variable;
    const evaluation::BuiltinRegistry& builtins;
    const mathematics::MathRegistry& mathematics;
    PlotProgram program;
    std::unordered_map<const void*, PlotRegister> registersByIdentity;
    PlotCompileStatus failure = PlotCompileStatus::Success;
};

[[nodiscard]] PlotRegister nextRegister(CompileContext& context, bool variableDependent) {
    const std::size_t index = context.program.variableDependent.size();
    if (index >= maxPlotInstructions || index > std::numeric_limits<PlotRegister>::max()) {
        context.failure = PlotCompileStatus::ResourceLimit;
        return 0;
    }
    context.program.variableDependent.push_back(variableDependent);
    return static_cast<PlotRegister>(index);
}

[[nodiscard]] std::optional<std::int64_t> exactIntegerExponent(const Expr& expression) {
    if (!expression.isNumber() || !expression.asNumber().isReal())
        return std::nullopt;
    const auto& value = expression.asNumber().asReal();
    if (!value.isInteger() || value.asInteger().bitLength() > 31)
        return std::nullopt;
    try {
        return std::stoll(value.asInteger().toString());
    }
    catch (...) {
        return std::nullopt;
    }
}

[[nodiscard]] std::optional<std::uint32_t> addConstant(
    CompileContext& context,
    PlotConstant constant) {
    if (context.program.constants.size() >= maxPlotConstants) {
        context.failure = PlotCompileStatus::ResourceLimit;
        return std::nullopt;
    }
    const std::size_t index = context.program.constants.size();
    context.program.constants.push_back(std::move(constant));
    return static_cast<std::uint32_t>(index);
}

[[nodiscard]] std::optional<PlotRegister> emitConstant(
    CompileContext& context,
    PlotConstant constant,
    const void* identity) {
    const auto index = addConstant(context, std::move(constant));
    if (!index)
        return std::nullopt;
    const PlotRegister destination = nextRegister(context, false);
    if (context.failure != PlotCompileStatus::Success)
        return std::nullopt;
    context.program.instructions.push_back(PlotInstruction{
        PlotOpcode::Constant, destination, 0, 0, *index, 0});
    context.registersByIdentity.emplace(identity, destination);
    return destination;
}

[[nodiscard]] bool isUnary(PlotOpcode opcode) noexcept {
    switch (opcode) {
    case PlotOpcode::Negate:
    case PlotOpcode::AngleConvert:
    case PlotOpcode::Abs:
    case PlotOpcode::Floor:
    case PlotOpcode::Ceil:
    case PlotOpcode::Trunc:
    case PlotOpcode::Round:
    case PlotOpcode::Frac:
    case PlotOpcode::Sign:
    case PlotOpcode::Sqrt:
    case PlotOpcode::Cbrt:
    case PlotOpcode::Exp:
    case PlotOpcode::Expm1:
    case PlotOpcode::Sinc:
    case PlotOpcode::Cosc:
    case PlotOpcode::Tanc:
    case PlotOpcode::Sinhc:
    case PlotOpcode::Tanhc:
    case PlotOpcode::Expc:
    case PlotOpcode::Log:
    case PlotOpcode::Log1p:
    case PlotOpcode::Log2:
    case PlotOpcode::Log10:
    case PlotOpcode::Sin:
    case PlotOpcode::Cos:
    case PlotOpcode::Tan:
    case PlotOpcode::Cot:
    case PlotOpcode::Sec:
    case PlotOpcode::Csc:
    case PlotOpcode::Asin:
    case PlotOpcode::Acos:
    case PlotOpcode::Atan:
    case PlotOpcode::Sinh:
    case PlotOpcode::Cosh:
    case PlotOpcode::Tanh:
    case PlotOpcode::Asinh:
    case PlotOpcode::Acosh:
    case PlotOpcode::Atanh:
    case PlotOpcode::Sech:
    case PlotOpcode::Erf:
    case PlotOpcode::Erfc:
    case PlotOpcode::FresnelC:
    case PlotOpcode::FresnelS:
    case PlotOpcode::ExponentialIntegralEi:
    case PlotOpcode::SineIntegralSi:
    case PlotOpcode::CosineIntegralCi:
        return true;
    default:
        return false;
    }
}

[[nodiscard]] std::optional<PlotOpcode> opcodeForBuiltin(BuiltinId id) {
    switch (id) {
    case BuiltinId::Subtract: return PlotOpcode::Subtract;
    case BuiltinId::Divide: return PlotOpcode::Divide;
    case BuiltinId::Negate: return PlotOpcode::Negate;
    case BuiltinId::Abs: return PlotOpcode::Abs;
    case BuiltinId::Floor: return PlotOpcode::Floor;
    case BuiltinId::Ceil: return PlotOpcode::Ceil;
    case BuiltinId::Trunc: return PlotOpcode::Trunc;
    case BuiltinId::Round: return PlotOpcode::Round;
    case BuiltinId::Frac: return PlotOpcode::Frac;
    case BuiltinId::Sign: return PlotOpcode::Sign;
    case BuiltinId::Clamp: return PlotOpcode::Clamp;
    case BuiltinId::Sqrt: return PlotOpcode::Sqrt;
    case BuiltinId::Cbrt: return PlotOpcode::Cbrt;
    case BuiltinId::Exp: return PlotOpcode::Exp;
    case BuiltinId::Expm1: return PlotOpcode::Expm1;
    case BuiltinId::Sinc: return PlotOpcode::Sinc;
    case BuiltinId::Cosc: return PlotOpcode::Cosc;
    case BuiltinId::Tanc: return PlotOpcode::Tanc;
    case BuiltinId::Sinhc: return PlotOpcode::Sinhc;
    case BuiltinId::Tanhc: return PlotOpcode::Tanhc;
    case BuiltinId::Expc: return PlotOpcode::Expc;
    case BuiltinId::Log: return PlotOpcode::Log;
    case BuiltinId::Log1p: return PlotOpcode::Log1p;
    case BuiltinId::Log2: return PlotOpcode::Log2;
    case BuiltinId::Log10: return PlotOpcode::Log10;
    case BuiltinId::Sin: return PlotOpcode::Sin;
    case BuiltinId::Cos: return PlotOpcode::Cos;
    case BuiltinId::Tan: return PlotOpcode::Tan;
    case BuiltinId::Cot: return PlotOpcode::Cot;
    case BuiltinId::Sec: return PlotOpcode::Sec;
    case BuiltinId::Csc: return PlotOpcode::Csc;
    case BuiltinId::Asin: return PlotOpcode::Asin;
    case BuiltinId::Acos: return PlotOpcode::Acos;
    case BuiltinId::Atan: return PlotOpcode::Atan;
    case BuiltinId::Sinh: return PlotOpcode::Sinh;
    case BuiltinId::Cosh: return PlotOpcode::Cosh;
    case BuiltinId::Tanh: return PlotOpcode::Tanh;
    case BuiltinId::Asinh: return PlotOpcode::Asinh;
    case BuiltinId::Acosh: return PlotOpcode::Acosh;
    case BuiltinId::Atanh: return PlotOpcode::Atanh;
    case BuiltinId::Sech: return PlotOpcode::Sech;
    case BuiltinId::Erf: return PlotOpcode::Erf;
    case BuiltinId::Erfc: return PlotOpcode::Erfc;
    case BuiltinId::FresnelC: return PlotOpcode::FresnelC;
    case BuiltinId::FresnelS: return PlotOpcode::FresnelS;
    case BuiltinId::ExponentialIntegralEi: return PlotOpcode::ExponentialIntegralEi;
    case BuiltinId::SineIntegralSi: return PlotOpcode::SineIntegralSi;
    case BuiltinId::CosineIntegralCi: return PlotOpcode::CosineIntegralCi;
    default: return std::nullopt;
    }
}

[[nodiscard]] std::optional<PlotRegister> compileExpression(
    CompileContext& context,
    const Expr& expression);

[[nodiscard]] std::optional<PlotRegister> emitBinary(
    CompileContext& context,
    PlotOpcode opcode,
    PlotRegister lhs,
    PlotRegister rhs,
    const void* identity) {
    const bool depends = context.program.variableDependent[lhs]
        || context.program.variableDependent[rhs];
    const PlotRegister destination = nextRegister(context, depends);
    if (context.failure != PlotCompileStatus::Success)
        return std::nullopt;
    context.program.instructions.push_back(
        PlotInstruction{opcode, destination, lhs, rhs, 0, 0});
    if (identity)
        context.registersByIdentity.emplace(identity, destination);
    return destination;
}

[[nodiscard]] std::optional<PlotRegister> compileVariadic(
    CompileContext& context,
    const Expr& expression,
    PlotOpcode opcode,
    const Rational& identityValue) {
    const auto& arguments = expression.asCall().arguments;
    if (arguments.empty())
        return emitConstant(context, PlotConstant{identityValue}, expression.identity());

    auto result = compileExpression(context, arguments.front());
    if (!result)
        return std::nullopt;
    for (std::size_t i = 1; i < arguments.size(); ++i) {
        const auto rhs = compileExpression(context, arguments[i]);
        if (!rhs)
            return std::nullopt;
        const auto combined = emitBinary(
            context, opcode, *result, *rhs,
            i + 1 == arguments.size() ? expression.identity() : nullptr);
        if (!combined)
            return std::nullopt;
        result = *combined;
    }
    if (arguments.size() == 1)
        context.registersByIdentity.emplace(expression.identity(), *result);
    return result;
}

[[nodiscard]] std::optional<PlotRegister> compileExpression(
    CompileContext& context,
    const Expr& expression) {
    if (const auto found = context.registersByIdentity.find(expression.identity());
        found != context.registersByIdentity.end())
        return found->second;

    if (expression.isNumber()) {
        const auto& number = expression.asNumber();
        if (!number.isReal()) {
            context.failure = PlotCompileStatus::NonRealConstant;
            return std::nullopt;
        }
        return emitConstant(
            context, PlotConstant{number.asReal().toRational()}, expression.identity());
    }

    if (expression.isDecimalApproximation())
        return emitConstant(
            context,
            PlotConstant{expression.asDecimalApproximation().displayedValue()},
            expression.identity());

    if (expression.isSymbol()) {
        if (expression.asSymbol().sameIdentity(context.variable)) {
            const PlotRegister destination = nextRegister(context, true);
            if (context.failure != PlotCompileStatus::Success)
                return std::nullopt;
            context.program.instructions.push_back(
                PlotInstruction{PlotOpcode::Variable, destination});
            context.registersByIdentity.emplace(expression.identity(), destination);
            return destination;
        }
        if (const auto* constant = context.mathematics.findConstant(expression.asSymbol()))
            return emitConstant(
                context, PlotConstant{constant->id}, expression.identity());
        context.failure = PlotCompileStatus::UnsupportedExpression;
        return std::nullopt;
    }

    if (!expression.isCall()) {
        context.failure = PlotCompileStatus::UnsupportedExpression;
        return std::nullopt;
    }

    const auto* definition = context.builtins.find(expression.asCall().head);
    if (!definition) {
        context.failure = PlotCompileStatus::UnsupportedExpression;
        return std::nullopt;
    }

    if (definition->id == BuiltinId::Add)
        return compileVariadic(context, expression, PlotOpcode::Add, Rational{BigInt{0}});
    if (definition->id == BuiltinId::Multiply)
        return compileVariadic(context, expression, PlotOpcode::Multiply, Rational{BigInt{1}});

    if (definition->id == BuiltinId::Power) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 2) {
            context.failure = PlotCompileStatus::UnsupportedExpression;
            return std::nullopt;
        }

        const auto base = compileExpression(context, arguments[0]);
        if (!base)
            return std::nullopt;

        // exact integer指数は積和だけで評価でき，負の底も実軸上で正確に扱える。
        if (const auto exponent = exactIntegerExponent(arguments[1])) {
            const PlotRegister destination = nextRegister(
                context, context.program.variableDependent[*base]);
            if (context.failure != PlotCompileStatus::Success)
                return std::nullopt;
            context.program.instructions.push_back(PlotInstruction{
                PlotOpcode::IntegerPower, destination, *base, 0, 0, *exponent});
            context.registersByIdentity.emplace(expression.identity(), destination);
            return destination;
        }

        // 非整数・変数指数もprincipal Powerとして専用opcodeへ落とす。
        // 実Plotとしてどのxで値が実数になるかはPlotAnalysis側がexactに証明し，
        // executorはそのdomain外をNonReal/Undefinedとして拒否する。
        const auto exponent = compileExpression(context, arguments[1]);
        if (!exponent)
            return std::nullopt;
        const PlotRegister destination = nextRegister(
            context,
            context.program.variableDependent[*base]
                || context.program.variableDependent[*exponent]);
        if (context.failure != PlotCompileStatus::Success)
            return std::nullopt;
        context.program.instructions.push_back(PlotInstruction{
            PlotOpcode::RealPower, destination, *base, *exponent});
        context.registersByIdentity.emplace(expression.identity(), destination);
        return destination;
    }

    if (definition->id == BuiltinId::UnitApplied) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 2 || !arguments[1].isString()) {
            context.failure = PlotCompileStatus::UnsupportedExpression;
            return std::nullopt;
        }
        const auto unit = mathematics::AngleSemantics::parseUnit(arguments[1].asString());
        if (!unit) {
            context.failure = PlotCompileStatus::UnsupportedExpression;
            return std::nullopt;
        }
        const auto magnitude = compileExpression(context, arguments[0]);
        if (!magnitude)
            return std::nullopt;
        const PlotRegister destination = nextRegister(
            context, context.program.variableDependent[*magnitude]);
        if (context.failure != PlotCompileStatus::Success)
            return std::nullopt;
        context.program.instructions.push_back(PlotInstruction{
            PlotOpcode::AngleConvert, destination, *magnitude, 0,
            static_cast<std::uint32_t>(*unit), 0});
        context.registersByIdentity.emplace(expression.identity(), destination);
        return destination;
    }

    if (definition->id == BuiltinId::Log && expression.asCall().arguments.size() == 2) {
        const auto& arguments = expression.asCall().arguments;
        const auto base = compileExpression(context, arguments[0]);
        const auto value = compileExpression(context, arguments[1]);
        if (!base || !value)
            return std::nullopt;

        const PlotRegister logBase = nextRegister(
            context, context.program.variableDependent[*base]);
        const PlotRegister logValue = nextRegister(
            context, context.program.variableDependent[*value]);
        if (context.failure != PlotCompileStatus::Success)
            return std::nullopt;
        context.program.instructions.push_back(
            PlotInstruction{PlotOpcode::Log, logBase, *base});
        context.program.instructions.push_back(
            PlotInstruction{PlotOpcode::Log, logValue, *value});
        return emitBinary(
            context, PlotOpcode::Divide, logValue, logBase, expression.identity());
    }

    if (definition->id == BuiltinId::Clamp) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() != 3) {
            context.failure = PlotCompileStatus::UnsupportedExpression;
            return std::nullopt;
        }
        const auto value = compileExpression(context, arguments[0]);
        const auto lower = compileExpression(context, arguments[1]);
        const auto upper = compileExpression(context, arguments[2]);
        if (!value || !lower || !upper)
            return std::nullopt;
        const bool depends = context.program.variableDependent[*value]
            || context.program.variableDependent[*lower]
            || context.program.variableDependent[*upper];
        const PlotRegister destination = nextRegister(context, depends);
        if (context.failure != PlotCompileStatus::Success)
            return std::nullopt;
        context.program.instructions.push_back(PlotInstruction{
            PlotOpcode::Clamp, destination, *value, *lower, *upper, 0});
        context.registersByIdentity.emplace(expression.identity(), destination);
        return destination;
    }

    const auto opcode = opcodeForBuiltin(definition->id);
    if (!opcode) {
        context.failure = PlotCompileStatus::UnsupportedExpression;
        return std::nullopt;
    }
    const auto& arguments = expression.asCall().arguments;
    if (isUnary(*opcode)) {
        if (arguments.size() != 1) {
            context.failure = PlotCompileStatus::UnsupportedExpression;
            return std::nullopt;
        }
        const auto argument = compileExpression(context, arguments[0]);
        if (!argument)
            return std::nullopt;
        const PlotRegister destination = nextRegister(
            context, context.program.variableDependent[*argument]);
        if (context.failure != PlotCompileStatus::Success)
            return std::nullopt;
        context.program.instructions.push_back(
            PlotInstruction{*opcode, destination, *argument});
        context.registersByIdentity.emplace(expression.identity(), destination);
        return destination;
    }

    if (arguments.size() != 2) {
        context.failure = PlotCompileStatus::UnsupportedExpression;
        return std::nullopt;
    }
    const auto lhs = compileExpression(context, arguments[0]);
    const auto rhs = compileExpression(context, arguments[1]);
    if (!lhs || !rhs)
        return std::nullopt;
    return emitBinary(context, *opcode, *lhs, *rhs, expression.identity());
}

[[nodiscard]] PlotCurveGeometryKind classifyGeometry(const PlotProgram& program) {
    enum class Degree : unsigned char { Constant = 0, Affine = 1, Quadratic = 2, Cubic = 3, Nonlinear = 4 };
    std::vector<Degree> degree(program.registerCount(), Degree::Nonlinear);

    const auto combineAdd = [](Degree lhs, Degree rhs) {
        if (lhs == Degree::Nonlinear || rhs == Degree::Nonlinear)
            return Degree::Nonlinear;
        return static_cast<Degree>(std::max(
            static_cast<unsigned char>(lhs),
            static_cast<unsigned char>(rhs)));
    };

    const auto scaleDegree = [](Degree value, std::int64_t factor) {
        if (value == Degree::Nonlinear)
            return Degree::Nonlinear;
        if (factor < 0)
            return value == Degree::Constant ? Degree::Constant : Degree::Nonlinear;
        const auto scaled = static_cast<std::int64_t>(value) * factor;
        if (scaled > static_cast<std::int64_t>(Degree::Cubic))
            return Degree::Nonlinear;
        return static_cast<Degree>(scaled);
    };

    for (const PlotInstruction& instruction : program.instructions) {
        switch (instruction.opcode) {
        case PlotOpcode::Constant:
            degree[instruction.destination] = Degree::Constant;
            break;
        case PlotOpcode::Variable:
            degree[instruction.destination] = Degree::Affine;
            break;
        case PlotOpcode::Add:
        case PlotOpcode::Subtract:
            degree[instruction.destination] = combineAdd(
                degree[instruction.lhs], degree[instruction.rhs]);
            break;
        case PlotOpcode::Multiply: {
            const Degree lhs = degree[instruction.lhs];
            const Degree rhs = degree[instruction.rhs];
            if (lhs == Degree::Constant)
                degree[instruction.destination] = rhs;
            else if (rhs == Degree::Constant)
                degree[instruction.destination] = lhs;
            else if (lhs == Degree::Nonlinear || rhs == Degree::Nonlinear)
                degree[instruction.destination] = Degree::Nonlinear;
            else {
                const auto combined = static_cast<unsigned char>(lhs)
                    + static_cast<unsigned char>(rhs);
                degree[instruction.destination] = combined <= static_cast<unsigned char>(Degree::Cubic)
                    ? static_cast<Degree>(combined)
                    : Degree::Nonlinear;
            }
            break;
        }
        case PlotOpcode::Divide:
            degree[instruction.destination] = degree[instruction.rhs] == Degree::Constant
                ? degree[instruction.lhs]
                : Degree::Nonlinear;
            break;
        case PlotOpcode::Negate:
        case PlotOpcode::AngleConvert:
            degree[instruction.destination] = degree[instruction.lhs];
            break;
        case PlotOpcode::IntegerPower:
            if (instruction.integer == 0 || degree[instruction.lhs] == Degree::Constant)
                degree[instruction.destination] = Degree::Constant;
            else
                degree[instruction.destination] = scaleDegree(
                    degree[instruction.lhs], instruction.integer);
            break;
        case PlotOpcode::RealPower:
            degree[instruction.destination] =
                degree[instruction.lhs] == Degree::Constant
                    && degree[instruction.rhs] == Degree::Constant
                ? Degree::Constant
                : Degree::Nonlinear;
            break;
        default:
            // 初等函数は引数が不変なら定数，変数依存なら一般曲線とする。
            degree[instruction.destination] = degree[instruction.lhs] == Degree::Constant
                ? Degree::Constant
                : Degree::Nonlinear;
            break;
        }
    }

    if (program.resultRegister < program.variableDependent.size()
        && program.variableDependent[program.resultRegister]) {
        for (auto it = program.instructions.rbegin(); it != program.instructions.rend(); ++it) {
            if (it->destination != program.resultRegister)
                continue;
            if (it->opcode == PlotOpcode::Floor
                || it->opcode == PlotOpcode::Ceil
                || it->opcode == PlotOpcode::Trunc
                || it->opcode == PlotOpcode::Round
                || it->opcode == PlotOpcode::Sign)
                return PlotCurveGeometryKind::PiecewiseConstant;
            break;
        }
    }

    if (program.resultRegister >= degree.size())
        return PlotCurveGeometryKind::Generic;
    switch (degree[program.resultRegister]) {
    case Degree::Constant: return PlotCurveGeometryKind::Constant;
    case Degree::Affine: return PlotCurveGeometryKind::Affine;
    case Degree::Quadratic: return PlotCurveGeometryKind::QuadraticPolynomial;
    case Degree::Cubic: return PlotCurveGeometryKind::CubicPolynomial;
    case Degree::Nonlinear: return PlotCurveGeometryKind::Generic;
    }
    return PlotCurveGeometryKind::Generic;
}

[[nodiscard]] BigFloat midpoint(
    const approximation::RealInterval& interval,
    std::size_t precisionBits) {
    if (interval.isPoint())
        return interval.lower().rounded(precisionBits, RoundingMode::NearestEven);
    const BigFloat sum = numeric::add(
        interval.lower(), interval.upper(), precisionBits + 2, RoundingMode::NearestEven);
    const BigFloat two = BigFloat::fromBigInt(
        BigInt{2}, precisionBits + 2, RoundingMode::NearestEven);
    return numeric::divide(sum, two, precisionBits, RoundingMode::NearestEven);
}

[[nodiscard]] BigFloat one(std::size_t bits) {
    return BigFloat::fromBigInt(BigInt{1}, bits, RoundingMode::NearestEven);
}

[[nodiscard]] BigInt floorRational(const Rational& value) {
    auto result = numeric::divmod(value.numerator(), value.denominator());
    if (!result.remainder.isZero() && value.numerator().isNegative())
        result.quotient -= BigInt{1};
    return result.quotient;
}

[[nodiscard]] BigInt ceilRational(const Rational& value) {
    auto result = numeric::divmod(value.numerator(), value.denominator());
    if (!result.remainder.isZero() && value.numerator().isPositive())
        result.quotient += BigInt{1};
    return result.quotient;
}

[[nodiscard]] BigInt truncRational(const Rational& value) {
    return value.numerator() / value.denominator();
}

[[nodiscard]] BigInt roundNearestEven(const Rational& value) {
    auto result = numeric::divmod(value.numerator(), value.denominator());
    if (result.remainder.isZero())
        return result.quotient;

    const BigInt twiceRemainder = result.remainder.abs() * BigInt{2};
    const auto comparison = twiceRemainder <=> value.denominator();
    bool awayFromZero = comparison == std::strong_ordering::greater;
    if (comparison == std::strong_ordering::equal)
        awayFromZero = !(result.quotient.abs() % BigInt{2}).isZero();
    if (awayFromZero)
        result.quotient += value.numerator().isNegative() ? BigInt{-1} : BigInt{1};
    return result.quotient;
}

[[nodiscard]] BigFloat integerPower(
    BigFloat base,
    std::int64_t exponent,
    std::size_t bits) {
    if (exponent == 0)
        return one(bits);

    const bool negative = exponent < 0;
    std::uint64_t remaining = static_cast<std::uint64_t>(negative ? -exponent : exponent);
    BigFloat result = one(bits);
    while (remaining != 0) {
        if ((remaining & 1U) != 0)
            result = numeric::multiply(result, base, bits, RoundingMode::NearestEven);
        remaining >>= 1U;
        if (remaining != 0)
            base = numeric::multiply(base, base, bits, RoundingMode::NearestEven);
    }
    if (!negative)
        return result;
    if (result.isZero())
        throw std::domain_error("Division by zero");
    return numeric::divide(one(bits), result, bits, RoundingMode::NearestEven);
}


[[nodiscard]] std::optional<std::int64_t> runtimeIntegerExponent(
    const BigFloat& value) {
    const Rational rational = value.toRational();
    if (!rational.isInteger() || rational.numerator().bitLength() > 31)
        return std::nullopt;
    try {
        return std::stoll(rational.numerator().toString());
    }
    catch (...) {
        return std::nullopt;
    }
}

[[nodiscard]] approximation::RealInterval radianArgument(
    const BigFloat& value,
    std::size_t bits,
    mathematics::AngleUnit unit) {
    const approximation::RealInterval input = approximation::RealInterval::point(value);
    if (unit == mathematics::AngleUnit::Radian)
        return input;

    const auto pi = approximation::enclosePi(bits).interval;
    const Rational scale = unit == mathematics::AngleUnit::Degree
        ? Rational{BigInt{1}, BigInt{180}}
        : Rational{BigInt{1}, BigInt{200}};
    const auto scaled = approximation::multiply(
        input, approximation::RealInterval::fromRational(scale, bits), bits);
    return approximation::multiply(scaled, pi, bits);
}


[[nodiscard]] approximation::RealInterval radiansToAngleInterval(
    const approximation::RealInterval& radians,
    std::size_t bits,
    mathematics::AngleUnit unit) {
    if (unit == mathematics::AngleUnit::Radian)
        return radians;
    const auto pi = approximation::enclosePi(bits).interval;
    const Rational multiplier = unit == mathematics::AngleUnit::Degree
        ? Rational{BigInt{180}}
        : Rational{BigInt{200}};
    const auto factor = approximation::divide(
        approximation::RealInterval::fromRational(multiplier, bits), pi, bits);
    return approximation::multiply(radians, factor, bits);
}

[[nodiscard]] approximation::RealInterval cardinalTrigInterval(
    PlotOpcode opcode,
    const BigFloat& value,
    std::size_t bits,
    mathematics::AngleUnit unit) {
    if (value.isZero()) {
        const Rational exact = opcode == PlotOpcode::Cosc
            ? Rational{BigInt{0}} : Rational{BigInt{1}};
        return approximation::RealInterval::fromRational(exact, bits);
    }

    const approximation::RealInterval radians = radianArgument(value, bits, unit);
    const auto sine = approximation::encloseSinRadianInterval(radians, bits).interval;
    const auto sinc = approximation::divide(sine, radians, bits);
    if (opcode == PlotOpcode::Sinc)
        return sinc;

    const auto cosine = approximation::encloseCosRadianInterval(radians, bits).interval;
    if (opcode == PlotOpcode::Tanc)
        return approximation::divide(sinc, cosine, bits);

    // (1-cos x)/x = sin(x/2) * sinc(x/2)。1-cosの桁落ちを避ける。
    const auto half = approximation::multiply(
        radians,
        approximation::RealInterval::fromRational(
            Rational{BigInt{1}, BigInt{2}}, bits),
        bits);
    const auto sineHalf = approximation::encloseSinRadianInterval(half, bits).interval;
    const auto sincHalf = approximation::divide(sineHalf, half, bits);
    return approximation::multiply(sineHalf, sincHalf, bits);
}

[[nodiscard]] approximation::RealInterval cardinalHyperbolicInterval(
    PlotOpcode opcode,
    const BigFloat& value,
    std::size_t bits) {
    if (value.isZero())
        return approximation::RealInterval::fromRational(Rational{BigInt{1}}, bits);

    const auto input = approximation::RealInterval::point(value);
    const auto sinh = approximation::encloseSinhReal(input, bits);
    const auto sinhc = approximation::divide(sinh, input, bits);
    if (opcode == PlotOpcode::Sinhc)
        return sinhc;
    if (opcode == PlotOpcode::Tanhc) {
        const auto cosh = approximation::encloseCoshReal(input, bits);
        return approximation::divide(sinhc, cosh, bits);
    }

    // expm1(x)/x = exp(x/2) * sinh(x/2)/(x/2)。subtractionを避ける。
    const auto half = approximation::multiply(
        input,
        approximation::RealInterval::fromRational(
            Rational{BigInt{1}, BigInt{2}}, bits),
        bits);
    const auto sinhHalf = approximation::encloseSinhReal(half, bits);
    const auto sinhcHalf = approximation::divide(sinhHalf, half, bits);
    const auto expHalf = approximation::encloseExp(half, bits).interval;
    return approximation::multiply(expHalf, sinhcHalf, bits);
}

[[nodiscard]] approximation::RealInterval reciprocalInterval(
    const approximation::RealInterval& value,
    std::size_t bits) {
    if (value.containsZero())
        throw std::domain_error("Division by zero");
    return approximation::divide(
        approximation::RealInterval::fromRational(Rational{BigInt{1}}, bits),
        value, bits);
}

[[nodiscard]] BigFloat constantValue(
    const PlotConstant& constant,
    std::size_t bits) {
    if (const auto* rational = std::get_if<Rational>(&constant.value))
        return BigFloat::fromRational(*rational, bits, RoundingMode::NearestEven);
    const auto enclosure = approximation::encloseConstant(
        std::get<ConstantId>(constant.value), bits + elementaryGuardBits);
    if (!enclosure)
        throw approximation::CertifiedBackendUnsupported{"Plot constant backend is unavailable"};
    return midpoint(enclosure->interval, bits);
}

} // namespace

PlotCompileResult compilePlotProgram(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    CompileContext context{variable, builtins, mathematics, {}, {}, PlotCompileStatus::Success};
    const auto result = compileExpression(context, expression);
    if (!result || context.failure != PlotCompileStatus::Success)
        return PlotCompileResult{context.failure, std::nullopt};
    context.program.resultRegister = *result;
    context.program.geometryKind = classifyGeometry(context.program);
    return PlotCompileResult{PlotCompileStatus::Success, std::move(context.program)};
}

BigFloatPlotExecutor::BigFloatPlotExecutor(
    const PlotProgram& program,
    std::size_t precisionBits,
    mathematics::AngleSemantics angles)
    : program_(program),
      precisionBits_(precisionBits),
      angles_(angles),
      registers_(program.registerCount()),
      invariantReady_(program.registerCount(), false) {
    if (precisionBits_ < 8)
        throw std::invalid_argument("Plot BigFloat precision must be at least 8 bits");
    initializeInvariants();
}

void BigFloatPlotExecutor::initializeInvariants() {
    for (const PlotInstruction& instruction : program_.instructions) {
        if (program_.variableDependent[instruction.destination])
            continue;
        const PlotNumericResult result = executeInstruction(instruction);
        if (!result.finite())
            throw std::domain_error("Plot invariant expression is not a finite real value");
        registers_[instruction.destination] = result.value;
        invariantReady_[instruction.destination] = true;
    }
}

BigFloat BigFloatPlotExecutor::evaluateFresnelCached(
    bool cosineIntegral,
    const BigFloat& input,
    std::size_t workBits) {
    if (input.isZero())
        return BigFloat::fromBigInt(
            BigInt{0}, precisionBits_, RoundingMode::NearestEven);

    const bool negative = input.isNegative();
    const BigFloat magnitude = negative ? -input : input;
    auto& cache = cosineIntegral ? fresnelCCache_ : fresnelSCache_;
    if (const auto found = cache.find(magnitude); found != cache.end())
        return negative ? -found->second : found->second;

    const auto point = approximation::RealInterval::point(magnitude);
    const auto enclosure = cosineIntegral
        ? approximation::encloseFresnelCReal(point, workBits)
        : approximation::encloseFresnelSReal(point, workBits);
    BigFloat value = midpoint(enclosure, precisionBits_);
    cache.emplace(magnitude, value);
    return negative ? -value : value;
}

PlotNumericResult BigFloatPlotExecutor::executeInstruction(
    const PlotInstruction& instruction) {
    const auto finite = [&](BigFloat value) {
        registers_[instruction.destination] = value;
        return PlotNumericResult{PlotNumericStatus::Finite, std::move(value)};
    };
    const std::size_t workBits = precisionBits_ + elementaryGuardBits;

    try {
        switch (instruction.opcode) {
        case PlotOpcode::Constant:
            return finite(constantValue(program_.constants[instruction.auxiliary], precisionBits_));
        case PlotOpcode::Variable:
            return PlotNumericResult{PlotNumericStatus::Unsupported, {}};
        case PlotOpcode::Add:
            return finite(numeric::add(
                registers_[instruction.lhs], registers_[instruction.rhs],
                precisionBits_, RoundingMode::NearestEven));
        case PlotOpcode::Subtract:
            return finite(numeric::subtract(
                registers_[instruction.lhs], registers_[instruction.rhs],
                precisionBits_, RoundingMode::NearestEven));
        case PlotOpcode::Multiply:
            return finite(numeric::multiply(
                registers_[instruction.lhs], registers_[instruction.rhs],
                precisionBits_, RoundingMode::NearestEven));
        case PlotOpcode::Divide:
            if (registers_[instruction.rhs].isZero())
                return PlotNumericResult{PlotNumericStatus::DivisionByZero, {}};
            return finite(numeric::divide(
                registers_[instruction.lhs], registers_[instruction.rhs],
                precisionBits_, RoundingMode::NearestEven));
        case PlotOpcode::Negate:
            return finite(-registers_[instruction.lhs]);
        case PlotOpcode::IntegerPower:
            if (registers_[instruction.lhs].isZero() && instruction.integer < 0)
                return PlotNumericResult{PlotNumericStatus::DivisionByZero, {}};
            if (registers_[instruction.lhs].isZero() && instruction.integer == 0)
                return PlotNumericResult{PlotNumericStatus::Undefined, {}};
            return finite(integerPower(
                registers_[instruction.lhs], instruction.integer, precisionBits_));
        case PlotOpcode::RealPower: {
            const BigFloat& base = registers_[instruction.lhs];
            const BigFloat& exponent = registers_[instruction.rhs];

            // mmCalのprincipal Power: 0^qはq>0だけ定義する。
            if (base.isZero()) {
                if (exponent.isZero())
                    return PlotNumericResult{PlotNumericStatus::Undefined, {}};
                if (exponent.isNegative())
                    return PlotNumericResult{PlotNumericStatus::DivisionByZero, {}};
                return finite(BigFloat::fromBigInt(
                    BigInt{0}, precisionBits_, RoundingMode::NearestEven));
            }

            // 負実数底のprincipal Powerが実数へ戻るのは実整数指数だけ。
            // floor/ceil等のinteger-valued指数もここでexactなBigFloat整数として扱える。
            if (base.isNegative()) {
                const auto integerExponent = runtimeIntegerExponent(exponent);
                if (!integerExponent)
                    return PlotNumericResult{PlotNumericStatus::NonReal, {}};
                return finite(integerPower(base, *integerExponent, precisionBits_));
            }

            if (exponent.isZero() || base == one(precisionBits_))
                return finite(one(precisionBits_));

            const auto logarithm = approximation::encloseLogPositive(
                approximation::RealInterval::point(base), workBits).interval;
            const auto scaled = approximation::multiply(
                logarithm,
                approximation::RealInterval::point(exponent),
                workBits);
            const auto powered = approximation::encloseExp(scaled, workBits).interval;
            return finite(midpoint(powered, precisionBits_));
        }
        case PlotOpcode::AngleConvert: {
            const auto sourceUnit = static_cast<mathematics::AngleUnit>(instruction.auxiliary);
            if (sourceUnit == angles_.defaultUnit())
                return finite(registers_[instruction.lhs]);
            const auto radians = radianArgument(
                registers_[instruction.lhs], workBits, sourceUnit);
            const auto converted = radiansToAngleInterval(
                radians, workBits, angles_.defaultUnit());
            return finite(midpoint(converted, precisionBits_));
        }
        case PlotOpcode::Abs:
            return finite(registers_[instruction.lhs].isNegative()
                ? -registers_[instruction.lhs]
                : registers_[instruction.lhs]);
        case PlotOpcode::Floor: {
            const BigInt value = floorRational(registers_[instruction.lhs].toRational());
            return finite(BigFloat::fromBigInt(value, precisionBits_, RoundingMode::NearestEven));
        }
        case PlotOpcode::Ceil: {
            const BigInt value = ceilRational(registers_[instruction.lhs].toRational());
            return finite(BigFloat::fromBigInt(value, precisionBits_, RoundingMode::NearestEven));
        }
        case PlotOpcode::Trunc: {
            const BigInt value = truncRational(registers_[instruction.lhs].toRational());
            return finite(BigFloat::fromBigInt(value, precisionBits_, RoundingMode::NearestEven));
        }
        case PlotOpcode::Round: {
            const BigInt value = roundNearestEven(registers_[instruction.lhs].toRational());
            return finite(BigFloat::fromBigInt(value, precisionBits_, RoundingMode::NearestEven));
        }
        case PlotOpcode::Frac: {
            const Rational value = registers_[instruction.lhs].toRational();
            const Rational fraction = value - Rational{floorRational(value)};
            return finite(BigFloat::fromRational(
                fraction, precisionBits_, RoundingMode::NearestEven));
        }
        case PlotOpcode::Sign: {
            const BigFloat& value = registers_[instruction.lhs];
            const BigInt sign = value.isZero() ? BigInt{0}
                : (value.isNegative() ? BigInt{-1} : BigInt{1});
            return finite(BigFloat::fromBigInt(sign, precisionBits_, RoundingMode::NearestEven));
        }
        case PlotOpcode::Clamp: {
            const BigFloat& value = registers_[instruction.lhs];
            const BigFloat& lower = registers_[instruction.rhs];
            const auto upperRegister = static_cast<PlotRegister>(instruction.auxiliary);
            if (upperRegister >= registers_.size())
                return PlotNumericResult{PlotNumericStatus::Unsupported, {}};
            const BigFloat& upper = registers_[upperRegister];
            if (upper < lower)
                return PlotNumericResult{PlotNumericStatus::Undefined, {}};
            return finite(value < lower ? lower : (value > upper ? upper : value));
        }
        case PlotOpcode::Sqrt: {
            if (registers_[instruction.lhs].isNegative())
                return PlotNumericResult{PlotNumericStatus::NonReal, {}};
            const auto enclosure = approximation::encloseSqrt(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits).interval;
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Cbrt: {
            const auto enclosure = approximation::encloseRealCubeRoot(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Exp: {
            const auto enclosure = approximation::encloseExp(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits).interval;
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Expm1: {
            const auto exponential = approximation::encloseExp(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits).interval;
            const auto enclosure = approximation::subtract(
                exponential,
                approximation::RealInterval::fromRational(Rational{BigInt{1}}, workBits),
                workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Sinc:
        case PlotOpcode::Cosc:
        case PlotOpcode::Tanc: {
            const auto enclosure = cardinalTrigInterval(
                instruction.opcode, registers_[instruction.lhs], workBits,
                angles_.defaultUnit());
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Sinhc:
        case PlotOpcode::Tanhc:
        case PlotOpcode::Expc: {
            const auto enclosure = cardinalHyperbolicInterval(
                instruction.opcode, registers_[instruction.lhs], workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Log:
        case PlotOpcode::Log2:
        case PlotOpcode::Log10: {
            const BigFloat& argument = registers_[instruction.lhs];
            if (argument.isZero())
                return PlotNumericResult{PlotNumericStatus::Undefined, {}};
            if (argument.isNegative())
                return PlotNumericResult{PlotNumericStatus::NonReal, {}};
            auto enclosure = approximation::encloseLogPositive(
                approximation::RealInterval::point(argument), workBits).interval;
            if (instruction.opcode != PlotOpcode::Log) {
                const Rational base = instruction.opcode == PlotOpcode::Log2
                    ? Rational{BigInt{2}} : Rational{BigInt{10}};
                const auto denominator = approximation::encloseLogPositive(
                    approximation::RealInterval::fromRational(base, workBits), workBits).interval;
                enclosure = approximation::divide(enclosure, denominator, workBits);
            }
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Log1p: {
            const auto shifted = numeric::add(
                registers_[instruction.lhs], one(workBits), workBits, RoundingMode::NearestEven);
            if (shifted.isZero())
                return PlotNumericResult{PlotNumericStatus::Undefined, {}};
            if (shifted.isNegative())
                return PlotNumericResult{PlotNumericStatus::NonReal, {}};
            const auto enclosure = approximation::encloseLogPositive(
                approximation::RealInterval::point(shifted), workBits).interval;
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Sin:
        case PlotOpcode::Cos:
        case PlotOpcode::Tan:
        case PlotOpcode::Cot:
        case PlotOpcode::Sec:
        case PlotOpcode::Csc: {
            const auto argument = radianArgument(
                registers_[instruction.lhs], workBits, angles_.defaultUnit());
            if (instruction.opcode == PlotOpcode::Sin) {
                const auto sine = approximation::encloseSinRadianInterval(argument, workBits).interval;
                return finite(midpoint(sine, precisionBits_));
            }
            if (instruction.opcode == PlotOpcode::Cos) {
                const auto cosine = approximation::encloseCosRadianInterval(argument, workBits).interval;
                return finite(midpoint(cosine, precisionBits_));
            }
            const auto sine = approximation::encloseSinRadianInterval(argument, workBits).interval;
            const auto cosine = approximation::encloseCosRadianInterval(argument, workBits).interval;
            approximation::RealInterval enclosure = [&] {
                switch (instruction.opcode) {
                case PlotOpcode::Tan: return approximation::divide(sine, cosine, workBits);
                case PlotOpcode::Cot: return approximation::divide(cosine, sine, workBits);
                case PlotOpcode::Sec: return reciprocalInterval(cosine, workBits);
                case PlotOpcode::Csc: return reciprocalInterval(sine, workBits);
                default: return sine;
                }
            }();
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Asin:
        case PlotOpcode::Acos:
        case PlotOpcode::Atan: {
            const auto input = approximation::RealInterval::point(registers_[instruction.lhs]);
            const Rational value = registers_[instruction.lhs].toRational();
            approximation::RealInterval radians = [&] {
                if (instruction.opcode == PlotOpcode::Atan)
                    return approximation::encloseAtan(input, workBits).interval;
                if (value < Rational{BigInt{-1}} || value > Rational{BigInt{1}})
                    throw std::domain_error("Inverse trigonometric argument is non-real");
                return instruction.opcode == PlotOpcode::Asin
                    ? approximation::encloseAsinRealRadian(input, workBits)
                    : approximation::encloseAcosRealRadian(input, workBits);
            }();
            const auto enclosure = radiansToAngleInterval(
                radians, workBits, angles_.defaultUnit());
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Sinh: {
            const auto enclosure = approximation::encloseSinhReal(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Cosh: {
            const auto enclosure = approximation::encloseCoshReal(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Tanh: {
            const auto enclosure = approximation::encloseTanhReal(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Asinh: {
            const auto enclosure = approximation::encloseAsinhReal(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Acosh: {
            if (registers_[instruction.lhs].toRational() < Rational{BigInt{1}})
                return PlotNumericResult{PlotNumericStatus::NonReal, {}};
            const auto enclosure = approximation::encloseAcoshReal(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Atanh: {
            const Rational value = registers_[instruction.lhs].toRational();
            if (value == Rational{BigInt{-1}} || value == Rational{BigInt{1}})
                return PlotNumericResult{PlotNumericStatus::Undefined, {}};
            if (value < Rational{BigInt{-1}} || value > Rational{BigInt{1}})
                return PlotNumericResult{PlotNumericStatus::NonReal, {}};
            const auto enclosure = approximation::encloseAtanhReal(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits);
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::Sech: {
            const auto cosine = approximation::encloseCoshReal(
                approximation::RealInterval::point(registers_[instruction.lhs]), workBits);
            return finite(midpoint(reciprocalInterval(cosine, workBits), precisionBits_));
        }
        case PlotOpcode::Erf:
        case PlotOpcode::Erfc:
        case PlotOpcode::ExponentialIntegralEi:
        case PlotOpcode::SineIntegralSi:
        case PlotOpcode::CosineIntegralCi: {
            const auto input = approximation::RealInterval::point(registers_[instruction.lhs]);
            if (instruction.opcode == PlotOpcode::ExponentialIntegralEi
                && registers_[instruction.lhs].isZero())
                return PlotNumericResult{PlotNumericStatus::Undefined, {}};
            if (instruction.opcode == PlotOpcode::CosineIntegralCi) {
                if (registers_[instruction.lhs].isZero())
                    return PlotNumericResult{PlotNumericStatus::Undefined, {}};
                if (registers_[instruction.lhs].isNegative())
                    return PlotNumericResult{PlotNumericStatus::NonReal, {}};
            }
            approximation::RealInterval enclosure = [&] {
                switch (instruction.opcode) {
                case PlotOpcode::Erf: return approximation::encloseErfReal(input, workBits);
                case PlotOpcode::Erfc: return approximation::encloseErfcReal(input, workBits);
                case PlotOpcode::ExponentialIntegralEi:
                    return approximation::encloseExponentialIntegralEiReal(input, workBits);
                case PlotOpcode::SineIntegralSi:
                    return approximation::encloseSineIntegralSiReal(input, workBits);
                case PlotOpcode::CosineIntegralCi:
                    return approximation::encloseCosineIntegralCiPositive(input, workBits);
                default: return input;
                }
            }();
            return finite(midpoint(enclosure, precisionBits_));
        }
        case PlotOpcode::FresnelC:
            return finite(evaluateFresnelCached(
                true, registers_[instruction.lhs], workBits));
        case PlotOpcode::FresnelS:
            return finite(evaluateFresnelCached(
                false, registers_[instruction.lhs], workBits));
        }
    }
    catch (const approximation::PrecisionInsufficient&) {
        return PlotNumericResult{PlotNumericStatus::PrecisionInsufficient, {}};
    }
    catch (const approximation::CertifiedBackendUnsupported&) {
        return PlotNumericResult{PlotNumericStatus::Unsupported, {}};
    }
    catch (const std::overflow_error&) {
        return PlotNumericResult{PlotNumericStatus::Overflow, {}};
    }
    catch (const std::domain_error&) {
        return PlotNumericResult{PlotNumericStatus::Undefined, {}};
    }
    return PlotNumericResult{PlotNumericStatus::Unsupported, {}};
}

PlotNumericResult BigFloatPlotExecutor::evaluate(const BigFloat& variableValue) {
    for (const PlotInstruction& instruction : program_.instructions) {
        if (!program_.variableDependent[instruction.destination])
            continue;
        if (instruction.opcode == PlotOpcode::Variable) {
            registers_[instruction.destination] = variableValue.rounded(
                precisionBits_, RoundingMode::NearestEven);
            continue;
        }
        const PlotNumericResult result = executeInstruction(instruction);
        if (!result.finite())
            return result;
    }
    return PlotNumericResult{
        PlotNumericStatus::Finite,
        registers_[program_.resultRegister]};
}

} // namespace mmcal::plot
