// 記号微分D
#include "differentiation.hpp"

#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "evaluation/iterator_spec.hpp"
#include "expression/array_utils.hpp"
#include "symbolic/substitution.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id) {
    return builtins.isCallTo(expression, id);
}

[[nodiscard]] Expr call(
    const evaluation::BuiltinRegistry& builtins,
    BuiltinId id,
    std::vector<Expr> arguments) {
    return Expr::call(builtins.symbol(id), std::move(arguments));
}

[[nodiscard]] Expr add(
    const evaluation::BuiltinRegistry& builtins,
    std::vector<Expr> arguments) {
    return call(builtins, BuiltinId::Add, std::move(arguments));
}

[[nodiscard]] Expr multiply(
    const evaluation::BuiltinRegistry& builtins,
    std::vector<Expr> arguments) {
    return call(builtins, BuiltinId::Multiply, std::move(arguments));
}

[[nodiscard]] Expr negate(
    const evaluation::BuiltinRegistry& builtins,
    Expr value) {
    return call(builtins, BuiltinId::Negate, {std::move(value)});
}

[[nodiscard]] Expr subtract(
    const evaluation::BuiltinRegistry& builtins,
    Expr lhs,
    Expr rhs) {
    return call(builtins, BuiltinId::Subtract, {std::move(lhs), std::move(rhs)});
}

[[nodiscard]] Expr divide(
    const evaluation::BuiltinRegistry& builtins,
    Expr numerator,
    Expr denominator) {
    return call(builtins, BuiltinId::Divide, {std::move(numerator), std::move(denominator)});
}

[[nodiscard]] Expr power(
    const evaluation::BuiltinRegistry& builtins,
    Expr base,
    Expr exponent) {
    return call(builtins, BuiltinId::Power, {std::move(base), std::move(exponent)});
}

[[nodiscard]] Expr pi(const mathematics::MathRegistry& mathematics) {
    const auto* definition = mathematics.findConstant(mathematics::ConstantId::Pi);
    if (!definition)
        error::throwCalcError(error::CalcErrorType::Internal, "Pi is not registered");
    return Expr{definition->symbol};
}

[[nodiscard]] Expr imaginaryUnit() {
    return Expr{Number::complex(numeric::RealNumber{}, numeric::RealNumber{BigInt{1}})};
}

[[nodiscard]] Expr angleConversionScale(
    BuiltinId id,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics) {
    switch (id) {
    case BuiltinId::DegreeToRadian:
        return divide(builtins, pi(mathematics), integer(180));
    case BuiltinId::DegreeToGradian:
        return divide(builtins, integer(10), integer(9));
    case BuiltinId::RadianToDegree:
        return divide(builtins, integer(180), pi(mathematics));
    case BuiltinId::RadianToGradian:
        return divide(builtins, integer(200), pi(mathematics));
    case BuiltinId::GradianToDegree:
        return divide(builtins, integer(9), integer(10));
    case BuiltinId::GradianToRadian:
        return divide(builtins, pi(mathematics), integer(200));
    default:
        error::throwCalcError(error::CalcErrorType::Internal, "Not an angle-conversion builtin");
    }
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

[[nodiscard]] bool containsVariable(
    const Expr& root,
    const expression::Symbol& variable) {
    std::vector<Expr> pending{root};
    while (!pending.empty()) {
        Expr current = std::move(pending.back());
        pending.pop_back();
        if (current.isSymbol() && current.asSymbol().sameIdentity(variable))
            return true;
        if (current.isCall()) {
            for (const Expr& argument : current.asCall().arguments)
                pending.push_back(argument);
        }
        else if (current.isArray()) {
            for (const Expr& element : current.asArray().elements)
                pending.push_back(element);
        }
        else if (current.isList()) {
            for (const Expr& element : current.asList().elements)
                pending.push_back(element);
        }
    }
    return false;
}

[[nodiscard]] std::optional<mathematics::AngleUnit> explicitAngleUnit(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (!isHead(expression, builtins, BuiltinId::UnitApplied))
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (arguments.size() != 2 || !arguments[1].isString())
        return std::nullopt;
    return mathematics::AngleSemantics::parseUnit(arguments[1].asString());
}

[[nodiscard]] Expr directTrigScale(
    const Expr& argument,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    mathematics::AngleUnit unit = angles.defaultUnit();
    if (const auto explicitUnit = explicitAngleUnit(argument, builtins))
        unit = *explicitUnit;

    switch (unit) {
    case mathematics::AngleUnit::Radian:
        return integer(1);
    case mathematics::AngleUnit::Degree:
        return divide(builtins, pi(mathematics), integer(180));
    case mathematics::AngleUnit::Gradian:
        return divide(builtins, pi(mathematics), integer(200));
    }
    return integer(1);
}

[[nodiscard]] Expr inverseScaledQuotient(
    Expr numerator,
    Expr denominator,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    switch (angles.defaultUnit()) {
    case mathematics::AngleUnit::Radian:
        return divide(builtins, std::move(numerator), std::move(denominator));
    case mathematics::AngleUnit::Degree:
        return divide(builtins,
            multiply(builtins, {integer(180), std::move(numerator)}),
            multiply(builtins, {pi(mathematics), std::move(denominator)}));
    case mathematics::AngleUnit::Gradian:
        return divide(builtins,
            multiply(builtins, {integer(200), std::move(numerator)}),
            multiply(builtins, {pi(mathematics), std::move(denominator)}));
    }
    return divide(builtins, std::move(numerator), std::move(denominator));
}

[[nodiscard]] Expr derivativeCore(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles);

[[nodiscard]] Expr unresolvedDerivative(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins) {
    return call(builtins, BuiltinId::Derivative, {expression, Expr{variable}});
}

[[nodiscard]] Expr chain(
    Expr outerDerivative,
    const Expr& argument,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    Expr innerDerivative = derivativeCore(
        argument, variable, builtins, mathematics, angles);
    return multiply(builtins, {std::move(outerDerivative), std::move(innerDerivative)});
}

[[nodiscard]] Expr derivativeCore(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!containsVariable(expression, variable))
        return integer(0);
    if (expression.isSymbol())
        return expression.asSymbol().sameIdentity(variable) ? integer(1) : integer(0);
    if (expression.isNumber() || expression.isBoolean() || expression.isString())
        return integer(0);

    if (expression.isArray()) {
        std::vector<Expr> elements;
        elements.reserve(expression.asArray().elements.size());
        for (const Expr& element : expression.asArray().elements)
            elements.push_back(derivativeCore(
                element, variable, builtins, mathematics, angles));
        return Expr::array(expression.asArray().shape, std::move(elements));
    }
    if (expression.isList()) {
        std::vector<Expr> elements;
        elements.reserve(expression.asList().elements.size());
        for (const Expr& element : expression.asList().elements)
            elements.push_back(derivativeCore(
                element, variable, builtins, mathematics, angles));
        return expression::braceValue(std::move(elements));
    }

    const auto* definition = builtins.find(expression.asCall().head);
    if (!definition)
        return unresolvedDerivative(expression, variable, builtins);
    const auto& a = expression.asCall().arguments;

    switch (definition->id) {
    case BuiltinId::Add: {
        std::vector<Expr> terms;
        terms.reserve(a.size());
        for (const Expr& term : a)
            terms.push_back(derivativeCore(term, variable, builtins, mathematics, angles));
        return add(builtins, std::move(terms));
    }
    case BuiltinId::Subtract:
        if (a.size() == 2)
            return subtract(builtins,
                derivativeCore(a[0], variable, builtins, mathematics, angles),
                derivativeCore(a[1], variable, builtins, mathematics, angles));
        break;
    case BuiltinId::Negate:
        if (a.size() == 1)
            return negate(builtins,
                derivativeCore(a[0], variable, builtins, mathematics, angles));
        break;
    case BuiltinId::Multiply: {
        std::vector<Expr> sumTerms;
        sumTerms.reserve(a.size());
        for (std::size_t i = 0; i < a.size(); ++i) {
            Expr di = derivativeCore(a[i], variable, builtins, mathematics, angles);
            if (di.isNumber() && di.asNumber().isZero())
                continue;
            std::vector<Expr> factors;
            factors.reserve(a.size());
            for (std::size_t j = 0; j < a.size(); ++j)
                factors.push_back(i == j ? di : a[j]);
            sumTerms.push_back(multiply(builtins, std::move(factors)));
        }
        return sumTerms.empty() ? integer(0) : add(builtins, std::move(sumTerms));
    }
    case BuiltinId::Divide:
        if (a.size() == 2) {
            // 分母が微分変数に依存しない場合は quotient rule を展開しない。
            // これにより d(f/c)=f'/c をそのまま保ち、不要な c/c^2 を生成しない。
            if (!containsVariable(a[1], variable))
                return divide(builtins,
                    derivativeCore(a[0], variable, builtins, mathematics, angles), a[1]);

            Expr du = derivativeCore(a[0], variable, builtins, mathematics, angles);
            Expr dv = derivativeCore(a[1], variable, builtins, mathematics, angles);
            Expr numerator = subtract(builtins,
                multiply(builtins, {std::move(du), a[1]}),
                multiply(builtins, {a[0], std::move(dv)}));
            return divide(builtins, std::move(numerator), power(builtins, a[1], integer(2)));
        }
        break;
    case BuiltinId::Power:
        if (a.size() == 2) {
            const bool baseDepends = containsVariable(a[0], variable);
            const bool exponentDepends = containsVariable(a[1], variable);
            if (!exponentDepends) {
                Expr du = derivativeCore(a[0], variable, builtins, mathematics, angles);
                return multiply(builtins, {
                    a[1],
                    power(builtins, a[0], subtract(builtins, a[1], integer(1))),
                    std::move(du)});
            }
            if (!baseDepends) {
                Expr dv = derivativeCore(a[1], variable, builtins, mathematics, angles);
                return multiply(builtins, {
                    expression,
                    call(builtins, BuiltinId::Log, {a[0]}),
                    std::move(dv)});
            }

            Expr du = derivativeCore(a[0], variable, builtins, mathematics, angles);
            Expr dv = derivativeCore(a[1], variable, builtins, mathematics, angles);
            return multiply(builtins, {
                expression,
                add(builtins, {
                    multiply(builtins, {std::move(dv), call(builtins, BuiltinId::Log, {a[0]})}),
                    divide(builtins, multiply(builtins, {a[1], std::move(du)}), a[0])
                })
            });
        }
        break;
    case BuiltinId::Cbrt:
        if (a.size() == 1)
            return divide(builtins,
                derivativeCore(a[0], variable, builtins, mathematics, angles),
                multiply(builtins, {
                    integer(3),
                    power(builtins, call(builtins, BuiltinId::Cbrt, {a[0]}), integer(2))}));
        break;

    case BuiltinId::Hypot:
        if (a.size() == 2) {
            Expr du = derivativeCore(a[0], variable, builtins, mathematics, angles);
            Expr dv = derivativeCore(a[1], variable, builtins, mathematics, angles);
            return divide(builtins,
                add(builtins, {
                    multiply(builtins, {a[0], std::move(du)}),
                    multiply(builtins, {a[1], std::move(dv)})}),
                expression);
        }
        break;

    case BuiltinId::Cis:
        if (a.size() == 1) {
            const Expr scale = directTrigScale(a[0], builtins, mathematics, angles);
            const Expr& differentialArgument = explicitAngleUnit(a[0], builtins)
                ? a[0].asCall().arguments[0]
                : a[0];
            return multiply(builtins, {
                scale,
                imaginaryUnit(),
                expression,
                derivativeCore(
                    differentialArgument, variable, builtins, mathematics, angles)});
        }
        break;

    case BuiltinId::Polar:
        if (a.size() == 2)
            return derivativeCore(
                multiply(builtins, {a[0], call(builtins, BuiltinId::Cis, {a[1]})}),
                variable, builtins, mathematics, angles);
        break;

    case BuiltinId::DegreeToRadian:
    case BuiltinId::DegreeToGradian:
    case BuiltinId::RadianToDegree:
    case BuiltinId::RadianToGradian:
    case BuiltinId::GradianToDegree:
    case BuiltinId::GradianToRadian:
        if (a.size() == 1)
            return multiply(builtins, {
                angleConversionScale(definition->id, builtins, mathematics),
                derivativeCore(a[0], variable, builtins, mathematics, angles)});
        break;

    case BuiltinId::Sqrt:
        if (a.size() == 1)
            return divide(builtins,
                derivativeCore(a[0], variable, builtins, mathematics, angles),
                multiply(builtins, {integer(2), call(builtins, BuiltinId::Sqrt, {a[0]})}));
        break;

    case BuiltinId::Sin:
    case BuiltinId::Cos:
    case BuiltinId::Tan:
    case BuiltinId::Cot:
    case BuiltinId::Sec:
    case BuiltinId::Csc:
        if (a.size() == 1) {
            const Expr scale = directTrigScale(a[0], builtins, mathematics, angles);
            Expr outer = integer(1);
            switch (definition->id) {
            case BuiltinId::Sin:
                outer = call(builtins, BuiltinId::Cos, {a[0]});
                break;
            case BuiltinId::Cos:
                outer = negate(builtins, call(builtins, BuiltinId::Sin, {a[0]}));
                break;
            case BuiltinId::Tan:
                outer = power(builtins, call(builtins, BuiltinId::Sec, {a[0]}), integer(2));
                break;
            case BuiltinId::Cot:
                outer = negate(builtins,
                    power(builtins, call(builtins, BuiltinId::Csc, {a[0]}), integer(2)));
                break;
            case BuiltinId::Sec:
                outer = multiply(builtins, {
                    call(builtins, BuiltinId::Sec, {a[0]}),
                    call(builtins, BuiltinId::Tan, {a[0]})});
                break;
            case BuiltinId::Csc:
                outer = negate(builtins, multiply(builtins, {
                    call(builtins, BuiltinId::Csc, {a[0]}),
                    call(builtins, BuiltinId::Cot, {a[0]})}));
                break;
            default:
                break;
            }
            const Expr& differentialArgument = explicitAngleUnit(a[0], builtins)
                ? a[0].asCall().arguments[0]
                : a[0];
            return multiply(builtins, {
                scale,
                std::move(outer),
                derivativeCore(
                    differentialArgument, variable, builtins, mathematics, angles)});
        }
        break;

    case BuiltinId::Asin:
    case BuiltinId::Acos:
    case BuiltinId::Atan:
        if (a.size() == 1) {
            Expr denominator = integer(1);
            if (definition->id == BuiltinId::Atan)
                denominator = add(builtins, {integer(1), power(builtins, a[0], integer(2))});
            else
                denominator = call(builtins, BuiltinId::Sqrt, {
                    subtract(builtins, integer(1), power(builtins, a[0], integer(2)))});
            Expr outer = inverseScaledQuotient(
                integer(1), std::move(denominator), builtins, mathematics, angles);
            if (definition->id == BuiltinId::Acos)
                outer = negate(builtins, std::move(outer));
            return chain(std::move(outer), a[0], variable, builtins, mathematics, angles);
        }
        break;

    case BuiltinId::Atan2:
        if (a.size() == 2) {
            // atan2[y,x]' = unitScale * (x y' - y x') / (x^2+y^2)
            Expr dy = derivativeCore(a[0], variable, builtins, mathematics, angles);
            Expr dx = derivativeCore(a[1], variable, builtins, mathematics, angles);
            Expr numerator = subtract(builtins,
                multiply(builtins, {a[1], std::move(dy)}),
                multiply(builtins, {a[0], std::move(dx)}));
            Expr denominator = add(builtins, {
                multiply(builtins, {a[1], a[1]}),
                multiply(builtins, {a[0], a[0]})});
            return inverseScaledQuotient(
                std::move(numerator), std::move(denominator),
                builtins, mathematics, angles);
        }
        break;

    case BuiltinId::Sinh:
        if (a.size() == 1)
            return chain(call(builtins, BuiltinId::Cosh, {a[0]}),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Cosh:
        if (a.size() == 1)
            return chain(call(builtins, BuiltinId::Sinh, {a[0]}),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Tanh:
        if (a.size() == 1)
            return chain(power(builtins, call(builtins, BuiltinId::Sech, {a[0]}), integer(2)),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Csch:
        if (a.size() == 1)
            return chain(negate(builtins, multiply(builtins, {
                    call(builtins, BuiltinId::Csch, {a[0]}),
                    call(builtins, BuiltinId::Coth, {a[0]})})),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Sech:
        if (a.size() == 1)
            return chain(negate(builtins, multiply(builtins, {
                    call(builtins, BuiltinId::Sech, {a[0]}),
                    call(builtins, BuiltinId::Tanh, {a[0]})})),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Coth:
        if (a.size() == 1)
            return chain(negate(builtins,
                    power(builtins, call(builtins, BuiltinId::Csch, {a[0]}), integer(2))),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Asinh:
        if (a.size() == 1)
            return chain(divide(builtins, integer(1), call(builtins, BuiltinId::Sqrt, {
                    add(builtins, {power(builtins, a[0], integer(2)), integer(1)})})),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Acosh:
        if (a.size() == 1)
            return chain(divide(builtins, integer(1), multiply(builtins, {
                    call(builtins, BuiltinId::Sqrt, {subtract(builtins, a[0], integer(1))}),
                    call(builtins, BuiltinId::Sqrt, {add(builtins, {a[0], integer(1)})})})),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Atanh:
        if (a.size() == 1)
            return chain(divide(builtins, integer(1),
                    subtract(builtins, integer(1), power(builtins, a[0], integer(2)))),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Log:
        if (a.size() == 1)
            return divide(builtins,
                derivativeCore(a[0], variable, builtins, mathematics, angles), a[0]);
        if (a.size() == 2) {
            // 底が微分変数に依存しない通常ケースは直接 u'/(u Log[b]) とし、
            // 後段Simplifierへ不要な Log[b]/Log[b]^2 の約分を押し付けない。
            if (!containsVariable(a[0], variable))
                return divide(
                    builtins,
                    derivativeCore(a[1], variable, builtins, mathematics, angles),
                    multiply(builtins, {a[1], call(builtins, BuiltinId::Log, {a[0]})}));

            // 底も変数なら change-of-base へ一度だけ展開して既存規則へ委譲する。
            const Expr changeOfBase = divide(
                builtins,
                call(builtins, BuiltinId::Log, {a[1]}),
                call(builtins, BuiltinId::Log, {a[0]}));
            return derivativeCore(changeOfBase, variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::Log2:
    case BuiltinId::Log10:
        if (a.size() == 1) {
            const Expr base = integer(definition->id == BuiltinId::Log2 ? 2 : 10);
            return divide(
                builtins,
                derivativeCore(a[0], variable, builtins, mathematics, angles),
                multiply(builtins, {a[0], call(builtins, BuiltinId::Log, {base})}));
        }
        break;
    case BuiltinId::Erf:
    case BuiltinId::Erfc:
        if (a.size() == 1) {
            Expr scale = divide(
                builtins, integer(2), call(builtins, BuiltinId::Sqrt, {pi(mathematics)}));
            if (definition->id == BuiltinId::Erfc)
                scale = negate(builtins, std::move(scale));
            Expr exponential = call(builtins, BuiltinId::Exp, {
                negate(builtins, power(builtins, a[0], integer(2)))});
            return chain(
                multiply(builtins, {std::move(scale), std::move(exponential)}),
                a[0], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::FresnelC:
    case BuiltinId::FresnelS:
        if (a.size() == 1) {
            // Fresnel C/S の定義に現れる角度は常にRadian。
            // sessionの既定角度単位へ依存させないため、内部sin/cosには明示Radを付ける。
            Expr phase = divide(
                builtins,
                multiply(builtins, {pi(mathematics), power(builtins, a[0], integer(2))}),
                integer(2));
            Expr radianPhase = call(builtins, BuiltinId::UnitApplied, {
                std::move(phase), Expr{std::string{"Rad"}}});
            Expr kernel = call(
                builtins,
                definition->id == BuiltinId::FresnelC ? BuiltinId::Cos : BuiltinId::Sin,
                {std::move(radianPhase)});
            return chain(std::move(kernel), a[0], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::Hypergeometric1F1:
        if (a.size() == 3
            && !containsVariable(a[0], variable)
            && !containsVariable(a[1], variable)) {
            // d/dz M(a,b,z) = (a/b) M(a+1,b+1,z)。
            // parameter微分は別の特殊函数知識を要するため、a/bが変数に依存する場合は保持する。
            Expr shifted = call(builtins, BuiltinId::Hypergeometric1F1, {
                add(builtins, {a[0], integer(1)}),
                add(builtins, {a[1], integer(1)}),
                a[2]});
            Expr kernel = multiply(builtins, {
                divide(builtins, a[0], a[1]),
                std::move(shifted)});
            return chain(std::move(kernel), a[2], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::Hypergeometric2F1:
        if (a.size() == 4
            && !containsVariable(a[0], variable)
            && !containsVariable(a[1], variable)
            && !containsVariable(a[2], variable)) {
            // d/dz 2F1(a,b;c;z)=(ab/c)2F1(a+1,b+1;c+1;z)。
            Expr shifted = call(builtins, BuiltinId::Hypergeometric2F1, {
                add(builtins, {a[0], integer(1)}),
                add(builtins, {a[1], integer(1)}),
                add(builtins, {a[2], integer(1)}),
                a[3]});
            Expr kernel = multiply(builtins, {
                divide(builtins, multiply(builtins, {a[0], a[1]}), a[2]),
                std::move(shifted)});
            return chain(std::move(kernel), a[3], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::EllipticF:
    case BuiltinId::EllipticE:
        if (a.size() == 2 && !containsVariable(a[1], variable)) {
            // Legendre incomplete elliptic integrals use a Radian amplitude independent of session angle mode.
            Expr radianAmplitude = call(builtins, BuiltinId::UnitApplied, {
                a[0], Expr{std::string{"Rad"}}});
            Expr sine = call(builtins, BuiltinId::Sin, {std::move(radianAmplitude)});
            Expr radicand = subtract(builtins, integer(1),
                multiply(builtins, {a[1], power(builtins, sine, integer(2))}));
            Expr root = call(builtins, BuiltinId::Sqrt, {std::move(radicand)});
            Expr kernel = definition->id == BuiltinId::EllipticF
                ? divide(builtins, integer(1), std::move(root))
                : std::move(root);
            return chain(std::move(kernel), a[0], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::EllipticPi:
        if (a.size() == 3
            && !containsVariable(a[0], variable)
            && !containsVariable(a[2], variable)) {
            Expr radianAmplitude = call(builtins, BuiltinId::UnitApplied, {
                a[1], Expr{std::string{"Rad"}}});
            Expr sine = call(builtins, BuiltinId::Sin, {std::move(radianAmplitude)});
            Expr sineSquared = power(builtins, sine, integer(2));
            Expr characteristic = subtract(builtins, integer(1),
                multiply(builtins, {a[0], sineSquared}));
            Expr radicand = subtract(builtins, integer(1),
                multiply(builtins, {a[2], sineSquared}));
            Expr denominator = multiply(builtins, {
                std::move(characteristic), call(builtins, BuiltinId::Sqrt, {std::move(radicand)})});
            return chain(divide(builtins, integer(1), std::move(denominator)),
                a[1], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::ExponentialIntegralEi:
        if (a.size() == 1) {
            Expr kernel = divide(builtins,
                call(builtins, BuiltinId::Exp, {a[0]}), a[0]);
            return chain(std::move(kernel), a[0], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::SineIntegralSi:
    case BuiltinId::CosineIntegralCi:
        if (a.size() == 1) {
            // Si/Ciの定義核はsession angle modeではなく常にRadian。
            Expr radian = angles.defaultUnit() == mathematics::AngleUnit::Radian
                ? a[0]
                : call(builtins, BuiltinId::UnitApplied, {a[0], Expr{std::string{"Rad"}}});
            const BuiltinId trig = definition->id == BuiltinId::SineIntegralSi
                ? BuiltinId::Sin : BuiltinId::Cos;
            Expr kernel = divide(builtins, call(builtins, trig, {std::move(radian)}), a[0]);
            return chain(std::move(kernel), a[0], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::LogarithmicIntegralLi:
        if (a.size() == 1) {
            Expr kernel = divide(builtins, integer(1), call(builtins, BuiltinId::Log, {a[0]}));
            return chain(std::move(kernel), a[0], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::Polylog:
        if (a.size() == 2 && !containsVariable(a[0], variable)) {
            Expr numerator = integer(0);
            if (a[0].isNumber() && a[0].asNumber().isReal()
                && a[0].asNumber().asReal().toRational() == Rational{BigInt{2}}) {
                // Li_1(z)=-Log(1-z) をここで直接使い、D結果を未評価polylog[1,z]へ戻さない。
                numerator = negate(builtins, call(builtins, BuiltinId::Log, {
                    subtract(builtins, integer(1), a[1])}));
            }
            else {
                Expr lowerOrder = subtract(builtins, a[0], integer(1));
                numerator = call(builtins, BuiltinId::Polylog, {std::move(lowerOrder), a[1]});
            }
            Expr kernel = divide(builtins, std::move(numerator), a[1]);
            return chain(std::move(kernel), a[1], variable, builtins, mathematics, angles);
        }
        break;

    case BuiltinId::Exp:
        if (a.size() == 1)
            return chain(expression, a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Expm1:
        if (a.size() == 1)
            return chain(call(builtins, BuiltinId::Exp, {a[0]}),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Log1p:
        if (a.size() == 1)
            return divide(builtins,
                derivativeCore(a[0], variable, builtins, mathematics, angles),
                add(builtins, {integer(1), a[0]}));
        break;

    case BuiltinId::SymbolicIntegral:
    case BuiltinId::Limit:
        if (a.size() == 2) {
            // 不定積分は原始函数代表元を意味するため、同じ積分変数での微分は被積分函数へ戻す。
            if (a[1].isSymbol() && a[1].asSymbol().sameIdentity(variable))
                return a[0];

            // Leibniz rule:
            // d/dx ∫_{lo(x)}^{hi(x)} f(t,x)dt
            // = f(hi,x)hi' - f(lo,x)lo' + ∫ ∂f/∂x dt
            if (const auto iterator = evaluation::parseRangeIteratorSpec(a[1])) {
                Expr upperDerivative = derivativeCore(
                    iterator->upper, variable, builtins, mathematics, angles);
                Expr lowerDerivative = derivativeCore(
                    iterator->lower, variable, builtins, mathematics, angles);
                Expr upperValue = substituteSymbol(a[0], iterator->variable, iterator->upper);
                Expr lowerValue = substituteSymbol(a[0], iterator->variable, iterator->lower);
                Expr boundary = subtract(builtins,
                    multiply(builtins, {std::move(upperValue), std::move(upperDerivative)}),
                    multiply(builtins, {std::move(lowerValue), std::move(lowerDerivative)}));

                Expr parameterDerivative = derivativeCore(
                    a[0], variable, builtins, mathematics, angles);
                if (parameterDerivative.isNumber() && parameterDerivative.asNumber().isZero())
                    return boundary;

                Expr remaining = call(builtins, BuiltinId::SymbolicIntegral, {
                    std::move(parameterDerivative), a[1]});
                return add(builtins, {std::move(boundary), std::move(remaining)});
            }
        }
        break;

    // これらは一般複素変数に関して正則でない、または現段階で
    // 厳密な微分規則を定義していない。誤ったscalar derivativeを返さない。
    case BuiltinId::Abs:
    case BuiltinId::Sign:
    case BuiltinId::Re:
    case BuiltinId::Im:
    case BuiltinId::Conj:
    case BuiltinId::Arg:
    case BuiltinId::Gamma:
    case BuiltinId::LogGamma:
    case BuiltinId::Beta:
    case BuiltinId::BetaLog:
    case BuiltinId::GeneralizedBinomial:
    case BuiltinId::FallingFactorial:
    case BuiltinId::RisingFactorial:
    case BuiltinId::RandSeed:
    case BuiltinId::Rand:
    case BuiltinId::RandInt:
    case BuiltinId::Choice:
    case BuiltinId::RandN:
    case BuiltinId::Factorial:
    case BuiltinId::Derivative:
    case BuiltinId::Floor:
    case BuiltinId::Ceil:
    case BuiltinId::Trunc:
    case BuiltinId::Round:
    case BuiltinId::Frac:
    case BuiltinId::Gcd:
    case BuiltinId::Lcm:
    case BuiltinId::Mod:
    case BuiltinId::Rem:
    case BuiltinId::Quotient:
    case BuiltinId::Permutation:
    case BuiltinId::Combination:
    case BuiltinId::Fibonacci:
    case BuiltinId::DiscreteFourierTransform:
    case BuiltinId::FastFourierTransform:
    case BuiltinId::InverseFourierTransform:
    case BuiltinId::Convolution:
    case BuiltinId::NextPow2:
    // Cardinal函数の通常の商微分はx=0に偽のholeを導入するため、piecewise/limit表現を持つまでは未評価Dとして保持する。
    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc:
    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc:
    case BuiltinId::Sum:
    case BuiltinId::Product:
    case BuiltinId::Min:
    case BuiltinId::Max:
    case BuiltinId::Mean:
    case BuiltinId::Median:
    case BuiltinId::Mode:
    case BuiltinId::Quantile:
    case BuiltinId::Percentile:
    case BuiltinId::VariancePopulation:
    case BuiltinId::VarianceSample:
    case BuiltinId::StddevPopulation:
    case BuiltinId::StddevSample:
    case BuiltinId::GeometricMean:
    case BuiltinId::HarmonicMean:
    case BuiltinId::Rms:
    case BuiltinId::MedianAbsoluteDeviation:
    case BuiltinId::MeanAbsoluteDeviation:
    case BuiltinId::Skewness:
    case BuiltinId::KurtosisPopulation:
    case BuiltinId::KurtosisSample:
    case BuiltinId::CoefficientVariation:
    case BuiltinId::StandardError:
    case BuiltinId::ZScore:
    case BuiltinId::Iqr:
    case BuiltinId::TrimMean:
    case BuiltinId::WinsorMean:
    case BuiltinId::Winsorized:
    case BuiltinId::Covariance:
    case BuiltinId::Correlation:
    case BuiltinId::SpearmanCorrelation:
    case BuiltinId::PercentRank:
    case BuiltinId::Dimensions:
    case BuiltinId::ArrayRank:
    case BuiltinId::ArrayGet:
    case BuiltinId::Reshape:
    case BuiltinId::Identity:
    case BuiltinId::Zeros:
    case BuiltinId::MatrixGet:
    case BuiltinId::Trace:
    case BuiltinId::Rows:
    case BuiltinId::Cols:
    case BuiltinId::Diag:
    case BuiltinId::VectorAdd:
    case BuiltinId::VectorSubtract:
    case BuiltinId::VectorScale:
    case BuiltinId::VectorDot:
    case BuiltinId::VectorCross:
    case BuiltinId::VectorNorm:
    case BuiltinId::VectorManhattan:
    case BuiltinId::VectorEuclidean:
    case BuiltinId::VectorNormalize:
    case BuiltinId::VectorProject:
    case BuiltinId::VectorAngle:
    case BuiltinId::VectorReflect:
    case BuiltinId::VectorReflectAxis:
    case BuiltinId::VectorSum:
    case BuiltinId::Transpose:
    case BuiltinId::MatrixAdd:
    case BuiltinId::MatrixMultiply:
    case BuiltinId::Determinant:
    case BuiltinId::Inverse:
    case BuiltinId::Rref:
    case BuiltinId::Rank:
    case BuiltinId::SolveLinear:
    case BuiltinId::NullSpace:
    case BuiltinId::LuDecomposition:
    case BuiltinId::QrDecomposition:
    case BuiltinId::SingularValueDecomposition:
    case BuiltinId::Eigenvalues:
    case BuiltinId::Eigenvectors:
    case BuiltinId::Eigensystem:
    case BuiltinId::ConjugateTranspose:
    case BuiltinId::Length:
    case BuiltinId::NumericDerivative:
    case BuiltinId::NumericIntegral:
    case BuiltinId::NumericalApproximation:
    case BuiltinId::Precision:
    case BuiltinId::Accuracy:
    case BuiltinId::Rationalize:
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify:
    case BuiltinId::Expand:
    case BuiltinId::Factor:
    case BuiltinId::Collect:
    case BuiltinId::Solve:
    case BuiltinId::Set:
    case BuiltinId::SetDelayed:
    case BuiltinId::Less:
    case BuiltinId::LessEqual:
    case BuiltinId::Greater:
    case BuiltinId::GreaterEqual:
    case BuiltinId::Equal:
    case BuiltinId::NotEqual:
    case BuiltinId::LogicalAnd:
    case BuiltinId::Element:
    case BuiltinId::If:
    case BuiltinId::History:
    case BuiltinId::InputHistory:
    case BuiltinId::OutputHistory:
    case BuiltinId::Exit:
    case BuiltinId::Clear:
    case BuiltinId::Definitions:
    case BuiltinId::Undefine:
    case BuiltinId::AngleMode:
    case BuiltinId::UnitApplied:
        break;
    }

    return unresolvedDerivative(expression, variable, builtins);
}

} // namespace

Expr differentiateExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    return simplify(
        derivativeCore(expression, variable, builtins, mathematics, angles),
        builtins, mathematics, angles);
}

} // namespace mmcal::symbolic
