#pragma once

#include "evaluation/builtin_registry.hpp"
#include "expression/expr.hpp"
#include "expression/symbol.hpp"
#include "mathematics/angle.hpp"
#include "mathematics/math_ids.hpp"
#include "mathematics/math_registry.hpp"
#include "numeric/big_float.hpp"
#include "numeric/rational.hpp"

#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <variant>
#include <vector>

namespace mmcal::plot {

using PlotRegister = std::uint32_t;

enum class PlotCurveGeometryKind {
    Generic,
    Constant,
    Affine,
    QuadraticPolynomial,
    CubicPolynomial,
    PiecewiseConstant
};

// PlotProgramは数値型に依存しないsampling用IR。
// BigFloat/double等のexecutorはこのopcode集合を共有し，Expr/evaluatorへ戻らない。
enum class PlotOpcode {
    Constant,
    Variable,
    Add,
    Subtract,
    Multiply,
    Divide,
    Negate,
    IntegerPower,
    RealPower,
    AngleConvert,
    Abs,
    Floor,
    Ceil,
    Trunc,
    Round,
    Frac,
    Sign,
    Clamp,
    Sqrt,
    Cbrt,
    Exp,
    Expm1,
    Sinc,
    Cosc,
    Tanc,
    Sinhc,
    Tanhc,
    Expc,
    Log,
    Log1p,
    Log2,
    Log10,
    Sin,
    Cos,
    Tan,
    Cot,
    Sec,
    Csc,
    Asin,
    Acos,
    Atan,
    Sinh,
    Cosh,
    Tanh,
    Asinh,
    Acosh,
    Atanh,
    Sech,
    Erf,
    Erfc,
    FresnelC,
    FresnelS,
    ExponentialIntegralEi,
    SineIntegralSi,
    CosineIntegralCi
};

struct PlotConstant final {
    std::variant<numeric::Rational, mathematics::ConstantId> value;

    [[nodiscard]] bool operator==(const PlotConstant&) const = default;
};

struct PlotInstruction final {
    PlotOpcode opcode = PlotOpcode::Constant;
    PlotRegister destination = 0;
    PlotRegister lhs = 0;
    PlotRegister rhs = 0;
    std::uint32_t auxiliary = 0;
    std::int64_t integer = 0;

    [[nodiscard]] bool operator==(const PlotInstruction&) const = default;
};

struct PlotProgram final {
    std::vector<PlotInstruction> instructions;
    std::vector<PlotConstant> constants;
    // registerごとの依存性。falseならexecutor生成時に一度だけ評価できる。
    std::vector<bool> variableDependent;
    PlotRegister resultRegister = 0;
    PlotCurveGeometryKind geometryKind = PlotCurveGeometryKind::Generic;

    [[nodiscard]] std::size_t registerCount() const noexcept {
        return variableDependent.size();
    }
};

enum class PlotCompileStatus {
    Success,
    UnsupportedExpression,
    NonRealConstant,
    UnsupportedExponent,
    ResourceLimit
};

struct PlotCompileResult final {
    PlotCompileStatus status = PlotCompileStatus::UnsupportedExpression;
    std::optional<PlotProgram> program;

    [[nodiscard]] explicit operator bool() const noexcept {
        return status == PlotCompileStatus::Success && program.has_value();
    }
};

// Exprをsampling専用IRへlowerする。未対応函数は式を壊さずcompile failureへ返す。
[[nodiscard]] PlotCompileResult compilePlotProgram(
    const expression::Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics);

enum class PlotNumericStatus {
    Finite,
    Undefined,
    NonReal,
    DivisionByZero,
    PrecisionInsufficient,
    Overflow,
    Unsupported
};

struct PlotNumericResult final {
    PlotNumericStatus status = PlotNumericStatus::Undefined;
    numeric::BigFloat value;

    [[nodiscard]] bool finite() const noexcept {
        return status == PlotNumericStatus::Finite;
    }
};

// 初版はBigFloat executor。IR自体にはBigFloatを埋め込まないため，後からdouble backendを追加できる。
// instanceごとにregister scratchを保持するため，evaluateはallocation-freeだがthread-safeではない。
class BigFloatPlotExecutor final {
public:
    BigFloatPlotExecutor(
        const PlotProgram& program,
        std::size_t precisionBits,
        mathematics::AngleSemantics angles = mathematics::defaultAngleSemantics());

    [[nodiscard]] PlotNumericResult evaluate(const numeric::BigFloat& variableValue);
    [[nodiscard]] std::size_t precisionBits() const noexcept { return precisionBits_; }

private:
    const PlotProgram& program_;
    std::size_t precisionBits_ = 0;
    mathematics::AngleSemantics angles_;
    std::vector<numeric::BigFloat> registers_;
    std::vector<bool> invariantReady_;
    // Fresnel C/S は奇函数で、対称rangeのadaptive samplingでは ±x をほぼ同じ順序で
    // 多数評価する。正側の64-bit plot値をexecutor lifetimeだけ再利用し、certified
    // backendを同じ |x| へ二度走らせない。
    std::map<numeric::BigFloat, numeric::BigFloat> fresnelCCache_;
    std::map<numeric::BigFloat, numeric::BigFloat> fresnelSCache_;

    [[nodiscard]] PlotNumericResult executeInstruction(const PlotInstruction& instruction);
    [[nodiscard]] numeric::BigFloat evaluateFresnelCached(
        bool cosineIntegral,
        const numeric::BigFloat& input,
        std::size_t workBits);
    void initializeInvariants();
};

} // namespace mmcal::plot
