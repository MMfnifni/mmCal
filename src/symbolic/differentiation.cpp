// 記号微分D
#include "symbolic/cases.hpp"
#include "differentiation.hpp"

#include "builtins/names.hpp"
#include "error/error_message.hpp"
#include "numeric/big_int.hpp"
#include "numeric/integer_algorithms.hpp"
#include "numeric/number.hpp"
#include "numeric/rational.hpp"
#include "simplification/simplification_context.hpp"
#include "simplification/simplifier.hpp"
#include "evaluation/iterator_spec.hpp"
#include "expression/array_utils.hpp"
#include "expression/exact_value.hpp"
#include "symbolic/algebra_transforms.hpp"
#include "symbolic/polynomial.hpp"
#include "symbolic/substitution.hpp"
#include "symbolic/series.hpp"

#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::symbolic {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using expression::exact::integer;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

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
        else if (current.isArray()
            && current.asArray().storageKind() == expression::ArrayStorageKind::Generic) {
            for (const Expr& element : current.asArray().storedExpressions())
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


[[nodiscard]] std::optional<std::uint64_t> boundedNonnegativeIntegerOrder(
    const Expr& expression,
    std::uint64_t maximum = 128) {
    if (!expression.isNumber() || !expression.asNumber().isReal()
        || !expression.asNumber().asReal().isInteger())
        return std::nullopt;
    const BigInt& value = expression.asNumber().asReal().asInteger();
    if (value.isNegative())
        return std::nullopt;
    const auto order = numeric::tryToUint64(value);
    if (!order || *order > maximum)
        return std::nullopt;
    return order;
}

[[nodiscard]] Expr finiteFactorialProductForDerivative(
    const Expr& x,
    std::uint64_t order,
    bool rising,
    const evaluation::BuiltinRegistry& builtins) {
    if (order == 0)
        return integer(1);
    std::vector<Expr> factors;
    factors.reserve(static_cast<std::size_t>(order));
    for (std::uint64_t k = 0; k < order; ++k) {
        Expr offset = integer(static_cast<std::int64_t>(k));
        factors.push_back(rising
            ? add(builtins, {x, std::move(offset)})
            : subtract(builtins, x, std::move(offset)));
    }
    return multiply(builtins, std::move(factors));
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
        const auto& array = expression.asArray();
        expression::ArrayBuilder builder;
        builder.reserve(array.size());
        for (std::size_t i = 0; i < array.size(); ++i) {
            if (array.storedKindAt(i) == expression::ArrayStorageKind::Generic)
                builder.append(derivativeCore(
                    array.expressionAt(i), variable, builtins, mathematics, angles));
            else
                builder.append(numeric::BigInt{0});
        }
        return Expr::array(builder.finish(array.shape));
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

            // u^u is defined only where the principal Power itself is defined.  On that
            // domain d(u^u)=u^u u' (Log[u]+1), so constructing u/u is unnecessary and
            // would leave a removable-looking but semantically delicate quotient in output.
            if (a[0] == a[1]) {
                Expr du = derivativeCore(a[0], variable, builtins, mathematics, angles);
                return multiply(builtins, {
                    expression, std::move(du),
                    add(builtins, {call(builtins, BuiltinId::Log, {a[0]}), integer(1)})});
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
    case BuiltinId::Sinc:
    case BuiltinId::Cosc:
    case BuiltinId::Tanc:
        if (a.size() == 1) {
            // Cardinal三角函数は0のremovable singularityを埋めた函数である。
            // 商微分だけでは0に偽のholeを作るため，各点極限をscalar ifで明示する。
            const Expr scale = directTrigScale(a[0], builtins, mathematics, angles);
            const Expr& differentialArgument = explicitAngleUnit(a[0], builtins)
                ? a[0].asCall().arguments[0]
                : a[0];
            Expr radians = multiply(builtins, {scale, differentialArgument});

            Expr trigArgument = angles.defaultUnit() == mathematics::AngleUnit::Radian
                    && !explicitAngleUnit(a[0], builtins)
                ? radians
                : call(builtins, BuiltinId::UnitApplied, {
                    radians, Expr{std::string{"Rad"}}});
            Expr sine = call(builtins, BuiltinId::Sin, {trigArgument});
            Expr cosine = call(builtins, BuiltinId::Cos, {trigArgument});

            Expr numerator = integer(0);
            Expr zeroDerivative = integer(0);
            if (definition->id == BuiltinId::Sinc) {
                numerator = subtract(builtins,
                    multiply(builtins, {radians, std::move(cosine)}),
                    std::move(sine));
                zeroDerivative = integer(0);
            } else if (definition->id == BuiltinId::Cosc) {
                numerator = add(builtins, {
                    multiply(builtins, {radians, std::move(sine)}),
                    std::move(cosine), integer(-1)});
                zeroDerivative = divide(builtins, scale, integer(2));
            } else {
                Expr secSquared = power(builtins,
                    call(builtins, BuiltinId::Sec, {std::move(trigArgument)}), integer(2));
                numerator = subtract(builtins,
                    multiply(builtins, {radians, std::move(secSquared)}),
                    call(builtins, BuiltinId::Tan, {angles.defaultUnit() == mathematics::AngleUnit::Radian
                            && !explicitAngleUnit(a[0], builtins)
                        ? radians
                        : call(builtins, BuiltinId::UnitApplied, {radians, Expr{std::string{"Rad"}}})}));
                zeroDerivative = integer(0);
            }

            Expr ordinary = divide(builtins,
                std::move(numerator), power(builtins, radians, integer(2)));
            const Expr innerDerivative = derivativeCore(
                differentialArgument, variable, builtins, mathematics, angles);
            ordinary = multiply(builtins, {scale, innerDerivative, std::move(ordinary)});
            zeroDerivative = multiply(builtins, {innerDerivative, std::move(zeroDerivative)});
            Expr nonzero = call(builtins, BuiltinId::NotEqual, {radians, integer(0)});
            Expr zeroCondition = call(builtins, BuiltinId::Equal, {radians, integer(0)});
            return detail::makeCases(builtins, {
                detail::makeCaseBranch(builtins, std::move(ordinary), std::move(nonzero)),
                detail::makeCaseBranch(builtins, std::move(zeroDerivative), std::move(zeroCondition))});
        }
        break;
    case BuiltinId::Sinhc:
    case BuiltinId::Tanhc:
    case BuiltinId::Expc:
        if (a.size() == 1) {
            // 非三角Cardinal函数も0での導函数を極限値として明示する。
            const Expr& u = a[0];
            Expr numerator = integer(0);
            Expr atZero = integer(0);
            if (definition->id == BuiltinId::Sinhc) {
                numerator = subtract(builtins,
                    multiply(builtins, {u, call(builtins, BuiltinId::Cosh, {u})}),
                    call(builtins, BuiltinId::Sinh, {u}));
                atZero = integer(0);
            } else if (definition->id == BuiltinId::Tanhc) {
                numerator = subtract(builtins,
                    multiply(builtins, {u, power(builtins, call(builtins, BuiltinId::Sech, {u}), integer(2))}),
                    call(builtins, BuiltinId::Tanh, {u}));
                atZero = integer(0);
            } else {
                Expr exponential = call(builtins, BuiltinId::Exp, {u});
                numerator = add(builtins, {
                    multiply(builtins, {u, exponential}),
                    negate(builtins, std::move(exponential)), integer(1)});
                atZero = divide(builtins, integer(1), integer(2));
            }
            const Expr innerDerivative = derivativeCore(
                u, variable, builtins, mathematics, angles);
            Expr ordinary = multiply(builtins, {
                innerDerivative,
                divide(builtins, std::move(numerator), power(builtins, u, integer(2)))});
            atZero = multiply(builtins, {innerDerivative, std::move(atZero)});
            Expr nonzero = call(builtins, BuiltinId::NotEqual, {u, integer(0)});
            Expr zeroCondition = call(builtins, BuiltinId::Equal, {u, integer(0)});
            return detail::makeCases(builtins, {
                detail::makeCaseBranch(builtins, std::move(ordinary), std::move(nonzero)),
                detail::makeCaseBranch(builtins, std::move(atZero), std::move(zeroCondition))});
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
            // Si/Ciの定義核は常にRadian。sincはsession angle modeを持つためRadを明示し，
            // Ciのcosだけは既定がRadianなら不要なUnitAppliedを省く。
            if (definition->id == BuiltinId::SineIntegralSi) {
                Expr radian = call(builtins, BuiltinId::UnitApplied, {
                    a[0], Expr{std::string{"Rad"}}});
                // sincがremovable singularityを埋めるので Si'(0)=1 も保持できる。
                return chain(call(builtins, BuiltinId::Sinc, {std::move(radian)}),
                    a[0], variable, builtins, mathematics, angles);
            }
            Expr radian = angles.defaultUnit() == mathematics::AngleUnit::Radian
                    && !explicitAngleUnit(a[0], builtins)
                ? a[0]
                : call(builtins, BuiltinId::UnitApplied, {
                    a[0], Expr{std::string{"Rad"}}});
            Expr kernel = divide(builtins,
                call(builtins, BuiltinId::Cos, {std::move(radian)}), a[0]);
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
            const Expr innerDerivative = derivativeCore(
                a[1], variable, builtins, mathematics, angles);
            Expr ordinary = multiply(builtins, {innerDerivative, std::move(kernel)});
            Expr nonzero = call(builtins, BuiltinId::NotEqual, {a[1], integer(0)});
            Expr zeroCondition = call(builtins, BuiltinId::Equal, {a[1], integer(0)});
            return detail::makeCases(builtins, {
                detail::makeCaseBranch(builtins, std::move(ordinary), std::move(nonzero)),
                detail::makeCaseBranch(builtins, innerDerivative, std::move(zeroCondition))});
        }
        break;

    case BuiltinId::LambertW:
        if ((a.size() == 1 || a.size() == 2)
            && (a.size() == 1 || !containsVariable(a[0], variable))) {
            const Expr& z = a.back();
            const Expr w = expression;
            // DLMF 4.13.4 gives two equivalent forms.  Exp[-W]/(1+W) is preferable
            // symbolically because it remains regular at W_0(0)=0 and therefore lets
            // higher derivatives proceed without an artificial x!=0 / x==0 Cases split.
            Expr kernel = divide(builtins,
                call(builtins, BuiltinId::Exp, {negate(builtins, w)}),
                add(builtins, {integer(1), w}));
            return multiply(builtins, {
                derivativeCore(z, variable, builtins, mathematics, angles),
                std::move(kernel)});
        }
        break;

    case BuiltinId::Gamma:
        if (a.size() == 1) {
            Expr kernel = multiply(builtins, {
                call(builtins, BuiltinId::Gamma, {a[0]}),
                call(builtins, BuiltinId::Digamma, {a[0]})});
            return chain(std::move(kernel), a[0], variable, builtins, mathematics, angles);
        }
        break;
    case BuiltinId::LogGamma:
        if (a.size() == 1)
            return chain(call(builtins, BuiltinId::Digamma, {a[0]}),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::Digamma:
        if (a.size() == 1)
            return chain(call(builtins, BuiltinId::Trigamma, {a[0]}),
                a[0], variable, builtins, mathematics, angles);
        break;
    case BuiltinId::IncompleteBeta:
        if (a.size() == 3
            && !containsVariable(a[0], variable)
            && !containsVariable(a[1], variable)) {
            Expr xPower = power(builtins, a[2], subtract(builtins, a[0], integer(1)));
            Expr oneMinusXPower = power(builtins,
                subtract(builtins, integer(1), a[2]),
                subtract(builtins, a[1], integer(1)));
            Expr kernel = divide(builtins,
                multiply(builtins, {std::move(xPower), std::move(oneMinusXPower)}),
                call(builtins, BuiltinId::Beta, {a[0], a[1]}));
            return chain(std::move(kernel), a[2], variable, builtins, mathematics, angles);
        }
        break;

    case BuiltinId::Beta:
    case BuiltinId::BetaLog:
        if (a.size() == 2) {
            // d log B(a,b) = da (psi(a)-psi(a+b)) + db (psi(b)-psi(a+b)).
            // Beta自体はこれへB(a,b)を掛ける。現在のBetaの正実数domainでも，
            // symbolic analytic continuation上でも同じ局所微分式を使える。
            Expr da = derivativeCore(a[0], variable, builtins, mathematics, angles);
            Expr db = derivativeCore(a[1], variable, builtins, mathematics, angles);
            Expr sum = add(builtins, {a[0], a[1]});
            Expr psiSum = call(builtins, BuiltinId::Digamma, {sum});
            Expr termA = multiply(builtins, {
                std::move(da),
                subtract(builtins,
                    call(builtins, BuiltinId::Digamma, {a[0]}), psiSum)});
            Expr termB = multiply(builtins, {
                std::move(db),
                subtract(builtins,
                    call(builtins, BuiltinId::Digamma, {a[1]}), std::move(psiSum))});
            Expr logarithmicDerivative = add(
                builtins, {std::move(termA), std::move(termB)});
            if (definition->id == BuiltinId::BetaLog)
                return logarithmicDerivative;
            return multiply(builtins, {
                expression, std::move(logarithmicDerivative)});
        }
        break;

    case BuiltinId::Fma:
        if (a.size() == 3) {
            // symbolic fma[a,b,c]は数学的にはa*b+c。数値評価時の単一丸め契約を
            // 微分ASTへ持ち込まず，exact algebraとして微分する。
            Expr expanded = add(builtins, {multiply(builtins, {a[0], a[1]}), a[2]});
            return derivativeCore(expanded, variable, builtins, mathematics, angles);
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
    case BuiltinId::GeneralizedBinomial:
    case BuiltinId::FallingFactorial:
    case BuiltinId::RisingFactorial:
        if (a.size() == 2) {
            const auto order = boundedNonnegativeIntegerOrder(a[1]);
            if (order) {
                Expr expanded = finiteFactorialProductForDerivative(
                    a[0], *order, definition->id == BuiltinId::RisingFactorial, builtins);
                if (definition->id == BuiltinId::GeneralizedBinomial)
                    expanded = divide(builtins, std::move(expanded),
                        Expr{Number{numeric::factorial(*order)}});
                return derivativeCore(
                    simplify(std::move(expanded), builtins, mathematics, angles),
                    variable, builtins, mathematics, angles);
            }
        }
        break;

    case BuiltinId::Zeta:
    case BuiltinId::Trigamma:
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
    case BuiltinId::BitAnd:
    case BuiltinId::BitOr:
    case BuiltinId::BitXor:
    case BuiltinId::BitNot:
    case BuiltinId::BitShiftLeft:
    case BuiltinId::BitShiftRight:
    case BuiltinId::BitLength:
    case BuiltinId::BitCount:
    case BuiltinId::BitGet:
    case BuiltinId::Clamp:
    case BuiltinId::Proj:
    case BuiltinId::Gcd:
    case BuiltinId::Lcm:
    case BuiltinId::Mod:
    case BuiltinId::Rem:
    case BuiltinId::Quotient:
    case BuiltinId::IsPrime:
    case BuiltinId::NextPrime:
    case BuiltinId::PreviousPrime:
    case BuiltinId::FactorInteger:
    case BuiltinId::Totient:
    case BuiltinId::Permutation:
    case BuiltinId::Combination:
    case BuiltinId::Fibonacci:
    case BuiltinId::DiscreteFourierTransform:
    case BuiltinId::FastFourierTransform:
    case BuiltinId::InverseFourierTransform:
    case BuiltinId::Convolution:
    case BuiltinId::NextPow2:
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
    case BuiltinId::Trace:
    case BuiltinId::Rows:
    case BuiltinId::Cols:
    case BuiltinId::Diag:
    case BuiltinId::VectorAdd:
    case BuiltinId::VectorSubtract:
    case BuiltinId::VectorScale:
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
    case BuiltinId::VectorInner:
    case BuiltinId::VectorOuter:
    case BuiltinId::VectorRejection:
    case BuiltinId::OrthogonalQ:
    case BuiltinId::OrthonormalQ:
    case BuiltinId::LinearIndependentQ:
    case BuiltinId::GramSchmidt:
    case BuiltinId::Gradient:
    case BuiltinId::Divergence:
    case BuiltinId::Curl:
    case BuiltinId::Laplacian:
    case BuiltinId::Jacobian:
    case BuiltinId::Hessian:
    case BuiltinId::DirectionalDerivative:
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
    case BuiltinId::ConditionNumber:
    case BuiltinId::LeastSquares:
    case BuiltinId::PseudoInverse:
    case BuiltinId::Eigenvalues:
    case BuiltinId::Eigenvectors:
    case BuiltinId::Eigensystem:
    case BuiltinId::ConjugateTranspose:
    case BuiltinId::Length:
    case BuiltinId::NumericDerivative:
    case BuiltinId::NumericIntegral:
    case BuiltinId::NumericalApproximation:
    case BuiltinId::Precision:
    case BuiltinId::Map:
    case BuiltinId::Range:
    case BuiltinId::Table:
    case BuiltinId::Accuracy:
    case BuiltinId::Explain:
    case BuiltinId::Rationalize:
    case BuiltinId::Root:
    case BuiltinId::Simplify:
    case BuiltinId::FullSimplify:
    case BuiltinId::Expand:
    case BuiltinId::Factor:
    case BuiltinId::Collect:
    case BuiltinId::Cases:
        {
            std::vector<Expr> branches;
            std::vector<Expr> boundaryConditions;
            branches.reserve(a.size());
            bool safe = true;
            bool strictVariableConditionSeen = false;
            for (const Expr& branchExpression : a) {
                if (!isHead(branchExpression, builtins, BuiltinId::CaseBranch)
                    || branchExpression.asCall().arguments.empty()
                    || branchExpression.asCall().arguments.size() > 2) {
                    safe = false;
                    break;
                }
                const auto& branch = branchExpression.asCall().arguments;
                std::optional<Expr> derivativeCondition;
                if (branch.size() == 2) {
                    derivativeCondition = branch[1];
                    if (containsVariable(branch[1], variable)) {
                        if (!branch[1].isCall()) {
                            safe = false;
                            break;
                        }
                        const auto* conditionDefinition = builtins.find(branch[1].asCall().head);
                        if (!conditionDefinition
                            || branch[1].asCall().arguments.size() != 2) {
                            safe = false;
                            break;
                        }

                        // casesの各枝をそのまま微分できるのは領域内部だけである。
                        // <=/>=の閉境界では隣接枝との接続を別途証明する必要があるため，
                        // 証明できない境界には値を捏造せずD[...]を明示的に残す。
                        if (conditionDefinition->id == BuiltinId::GreaterEqual
                            || conditionDefinition->id == BuiltinId::LessEqual) {
                            const auto& conditionArguments = branch[1].asCall().arguments;
                            const BuiltinId strictId = conditionDefinition->id == BuiltinId::GreaterEqual
                                ? BuiltinId::Greater : BuiltinId::Less;
                            derivativeCondition = call(builtins, strictId, {
                                conditionArguments[0], conditionArguments[1]});
                            Expr boundary = call(builtins, BuiltinId::Equal, {
                                conditionArguments[0], conditionArguments[1]});
                            if (std::find(boundaryConditions.begin(), boundaryConditions.end(), boundary)
                                == boundaryConditions.end())
                                boundaryConditions.push_back(std::move(boundary));
                        }
                        else if (conditionDefinition->id == BuiltinId::Greater
                            || conditionDefinition->id == BuiltinId::Less) {
                            strictVariableConditionSeen = true;
                        }
                        else {
                            // Equal/NotEqualはpolylogやcardinal函数の可除特異点を
                            // 表すことがあり，孤立枝の値を微分しても周囲の函数の導函数には
                            // ならない。複合条件も境界解析を要するため未評価のまま保持する。
                            safe = false;
                            break;
                        }
                    }
                }

                if (branch.size() == 1 && !boundaryConditions.empty()) {
                    for (Expr& boundary : boundaryConditions)
                        branches.push_back(detail::makeCaseBranch(
                            builtins, unresolvedDerivative(expression, variable, builtins),
                            std::move(boundary)));
                    boundaryConditions.clear();
                }

                Expr branchDerivative = branch.size() == 1 && strictVariableConditionSeen
                    ? unresolvedDerivative(expression, variable, builtins)
                    : derivativeCore(branch[0], variable, builtins, mathematics, angles);
                // x依存の狭義不等式の後にあるdefault枝にはその境界点が含まれ得る。
                // そこでdefault値を直接微分すると，abs[x]型の非微分可能点へ偽の値を
                // 与えるため，補集合を内部と境界へ分割できない場合はD[...]を保持する。
                branches.push_back(detail::makeCaseBranch(
                    builtins, std::move(branchDerivative), std::move(derivativeCondition)));
            }
            if (safe) {
                for (Expr& boundary : boundaryConditions)
                    branches.push_back(detail::makeCaseBranch(
                        builtins, unresolvedDerivative(expression, variable, builtins),
                        std::move(boundary)));
                return detail::makeCases(builtins, std::move(branches));
            }
        }
        break;
    case BuiltinId::CaseBranch:
    case BuiltinId::Solve:
    case BuiltinId::GroebnerBasis:
    case BuiltinId::PolynomialReduce:
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
    case BuiltinId::Series:
        break;
    case BuiltinId::SeriesData:
        if (const auto series = parseSeriesData(expression, builtins))
            if (auto derivative = differentiateSeriesExpression(
                    *series, variable, builtins, mathematics, angles))
                return *derivative;
        break;
    case BuiltinId::Normal:
    case BuiltinId::ToNormal:
    case BuiltinId::UnitApplied:
        break;
    }

    return unresolvedDerivative(expression, variable, builtins);
}

[[nodiscard]] Expr polylogOrderShift(
    const Expr& order,
    std::uint64_t shift,
    const Expr& z,
    const evaluation::BuiltinRegistry& builtins) {
    Expr lowered = order;
    if (shift != 0) {
        if (order.isNumber() && order.asNumber().isReal()) {
            const Rational shifted = order.asNumber().asReal().toRational()
                - Rational{BigInt::fromUnsigned(shift)};
            lowered = Expr{Number{numeric::RealNumber{shifted}}};
        }
        else
            lowered = subtract(builtins, order, integer(static_cast<std::int64_t>(shift)));
    }
    if (lowered.isNumber() && lowered.asNumber().isReal()) {
        const auto rationalOrder = lowered.asNumber().asReal().toRational();
        if (rationalOrder == Rational{BigInt{1}})
            return negate(builtins, call(builtins, BuiltinId::Log, {
                subtract(builtins, integer(1), z)}));
        if (rationalOrder == Rational{BigInt{0}})
            return divide(builtins, z, subtract(builtins, integer(1), z));
    }
    return call(builtins, BuiltinId::Polylog, {std::move(lowered), z});
}

[[nodiscard]] std::optional<Expr> repeatedQuadraticExponentialDerivative(
    const Expr& expression,
    const expression::Symbol& variable,
    std::uint64_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (order == 0 || order > 64
        || !isHead(expression, builtins, BuiltinId::Exp)
        || expression.asCall().arguments.size() != 1)
        return std::nullopt;
    const Expr& exponent = expression.asCall().arguments[0];
    if (!containsVariable(exponent, variable))
        return std::nullopt;

    // exact quadratic exponentではP_nを式木として微分・展開せず，
    // 係数vector上で P_(n+1)=P'_n+q'P_n を更新する。
    if (const auto q = toRationalPolynomial(
            exponent, variable, builtins, PolynomialConversionOptions{3, 16});
        q && q->degree() <= 2) {
        const Rational q0 = q->coefficient(0);
        static_cast<void>(q0);
        const Rational linear = q->coefficient(1);
        const Rational quadratic = q->coefficient(2);
        const Rational slope = Rational{BigInt{2}} * quadratic;
        std::vector<Rational> coefficients{Rational{BigInt{1}}};
        for (std::uint64_t n = 0; n < order; ++n) {
            std::vector<Rational> next(coefficients.size() + 1, Rational{BigInt{0}});
            for (std::size_t i = 0; i < coefficients.size(); ++i) {
                if (i != 0)
                    next[i - 1] += Rational{BigInt::fromUnsigned(i)} * coefficients[i];
                next[i] += linear * coefficients[i];
                next[i + 1] += slope * coefficients[i];
            }
            while (next.size() > 1 && next.back().isZero())
                next.pop_back();
            coefficients = std::move(next);
        }
        RationalPolynomial polynomial{std::move(coefficients)};
        Expr factor = polynomialToExpandedExpr(polynomial, variable, builtins);
        return multiply(builtins, {std::move(factor), expression});
    }

    Expr first = simplify(
        derivativeCore(exponent, variable, builtins, mathematics, angles),
        builtins, mathematics, angles);
    Expr second = simplify(
        derivativeCore(first, variable, builtins, mathematics, angles),
        builtins, mathematics, angles);
    Expr third = simplify(
        derivativeCore(second, variable, builtins, mathematics, angles),
        builtins, mathematics, angles);
    if (!(third.isNumber() && third.asNumber().isZero()))
        return std::nullopt;

    // exp[q(x)] with third derivative zero: P_0=1, P_(n+1)=P'_n+q'P_n.
    // The recurrence keeps the common exponential factor outside instead of
    // repeatedly expanding it through the generic product rule.
    Expr polynomial = integer(1);
    for (std::uint64_t n = 0; n < order; ++n) {
        Expr derivative = simplify(
            derivativeCore(polynomial, variable, builtins, mathematics, angles),
            builtins, mathematics, angles);
        polynomial = simplify(
            add(builtins, {
                std::move(derivative),
                multiply(builtins, {first, polynomial})}),
            builtins, mathematics, angles);
    }
    polynomial = expandExpression(polynomial, builtins, mathematics, angles, {256});
    return multiply(builtins, {std::move(polynomial), expression});
}

[[nodiscard]] std::optional<Expr> repeatedDirectLambertDerivative(
    const Expr& expression,
    const expression::Symbol& variable,
    std::uint64_t order,
    const evaluation::BuiltinRegistry& builtins) {
    if (order == 0 || order > 64
        || !isHead(expression, builtins, BuiltinId::LambertW))
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if ((arguments.size() != 1 && arguments.size() != 2)
        || !arguments.back().isSymbol()
        || !arguments.back().asSymbol().sameIdentity(variable)
        || (arguments.size() == 2 && containsVariable(arguments[0], variable)))
        return std::nullopt;

    // DLMF 4.13.4_1--4.13.4_2:
    // D^n W = exp[-n W] p_(n-1)(W)/(1+W)^(2n-1),
    // p_0=1, p_n=(1+W)p'_(n-1)+(1-n(W+3))p_(n-1).
    // 係数をexact BigIntで更新し，x=0のremovable singularityを人工的に作らない。
    std::vector<BigInt> coefficients{BigInt{1}};
    for (std::uint64_t n = 1; n < order; ++n) {
        std::vector<BigInt> next(coefficients.size() + 1, BigInt{0});
        const BigInt nBig = BigInt::fromUnsigned(n);
        const BigInt constantFactor = BigInt{1} - BigInt::fromUnsigned(3 * n);
        for (std::size_t i = 0; i < coefficients.size(); ++i) {
            const BigInt& coefficient = coefficients[i];
            next[i] += constantFactor * coefficient;
            next[i + 1] -= nBig * coefficient;
            if (i == 0)
                continue;
            const BigInt iBig = BigInt::fromUnsigned(i);
            next[i - 1] += iBig * coefficient;
            next[i] += iBig * coefficient;
        }
        coefficients = std::move(next);
    }

    const Expr w = expression;
    std::vector<Expr> polynomialTerms;
    polynomialTerms.reserve(coefficients.size());
    for (std::size_t i = 0; i < coefficients.size(); ++i) {
        const BigInt& coefficient = coefficients[i];
        if (coefficient.isZero())
            continue;
        Expr monomial = i == 0
            ? integer(1)
            : (i == 1 ? w : power(builtins, w, integer(static_cast<std::int64_t>(i))));
        if (coefficient == BigInt{-1})
            monomial = negate(builtins, std::move(monomial));
        else if (!(coefficient == BigInt{1}))
            monomial = multiply(builtins, {
                Expr{Number{coefficient}}, std::move(monomial)});
        polynomialTerms.push_back(std::move(monomial));
    }
    Expr polynomial = polynomialTerms.size() == 1
        ? std::move(polynomialTerms.front())
        : add(builtins, std::move(polynomialTerms));
    Expr exponential = call(builtins, BuiltinId::Exp, {
        multiply(builtins, {
            integer(-static_cast<std::int64_t>(order)), w})});
    Expr denominator = power(builtins,
        add(builtins, {integer(1), w}),
        integer(static_cast<std::int64_t>(2 * order - 1)));
    return divide(builtins,
        multiply(builtins, {std::move(polynomial), std::move(exponential)}),
        std::move(denominator));
}

[[nodiscard]] bool isAdditiveExpression(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    return isHead(expression, builtins, BuiltinId::Add)
        || isHead(expression, builtins, BuiltinId::Subtract)
        || isHead(expression, builtins, BuiltinId::Negate);
}

void appendDerivativeLinearTerms(
    const Expr& expression,
    Rational coefficient,
    std::vector<Expr>& terms,
    const evaluation::BuiltinRegistry& builtins,
    std::size_t depth = 0) {
    if (coefficient.isZero())
        return;

    if (depth <= 64 && isHead(expression, builtins, BuiltinId::Add)) {
        for (const Expr& argument : expression.asCall().arguments)
            appendDerivativeLinearTerms(
                argument, coefficient, terms, builtins, depth + 1);
        return;
    }
    if (depth <= 64 && isHead(expression, builtins, BuiltinId::Subtract)
        && expression.asCall().arguments.size() == 2) {
        appendDerivativeLinearTerms(
            expression.asCall().arguments[0], coefficient, terms, builtins, depth + 1);
        appendDerivativeLinearTerms(
            expression.asCall().arguments[1], -coefficient, terms, builtins, depth + 1);
        return;
    }
    if (depth <= 64 && isHead(expression, builtins, BuiltinId::Negate)
        && expression.asCall().arguments.size() == 1) {
        appendDerivativeLinearTerms(
            expression.asCall().arguments[0], -coefficient, terms, builtins, depth + 1);
        return;
    }

    // exactな有理scalarだけを加法子へ分配する。一般expandを呼ばず，
    // 高階product ruleで生じる2(A-B)のような形だけを平坦化する。
    if (depth <= 64 && isHead(expression, builtins, BuiltinId::Multiply)) {
        Rational scalar = coefficient;
        const Expr* additive = nullptr;
        bool compatible = true;
        for (const Expr& factor : expression.asCall().arguments) {
            if (factor.isNumber() && factor.asNumber().isReal()) {
                scalar *= factor.asNumber().asReal().toRational();
                continue;
            }
            if (!additive && isAdditiveExpression(factor, builtins)) {
                additive = &factor;
                continue;
            }
            compatible = false;
            break;
        }
        if (compatible && additive) {
            appendDerivativeLinearTerms(
                *additive, scalar, terms, builtins, depth + 1);
            return;
        }
    }

    if (coefficient == Rational{BigInt{1}})
        terms.push_back(expression);
    else if (coefficient == Rational{BigInt{-1}})
        terms.push_back(negate(builtins, expression));
    else
        terms.push_back(multiply(builtins, {
            Expr{Number{coefficient}}, expression}));
}

[[nodiscard]] std::optional<Expr> repeatedDirectPolylogDerivative(
    const Expr& expression,
    const expression::Symbol& variable,
    std::uint64_t order,
    const evaluation::BuiltinRegistry& builtins) {
    if (order == 0 || order > 64
        || !isHead(expression, builtins, BuiltinId::Polylog)
        || expression.asCall().arguments.size() != 2)
        return std::nullopt;
    const auto& arguments = expression.asCall().arguments;
    if (containsVariable(arguments[0], variable)
        || !arguments[1].isSymbol()
        || !arguments[1].asSymbol().sameIdentity(variable))
        return std::nullopt;

    // D^n = z^-n * theta(theta-1)...(theta-n+1), theta=z D.
    // theta^k Li_s = Li_{s-k}; coefficients are signed Stirling numbers s(n,k).
    std::vector<BigInt> stirling(static_cast<std::size_t>(order + 1), BigInt{0});
    stirling[0] = BigInt{1};
    for (std::uint64_t n = 1; n <= order; ++n) {
        std::vector<BigInt> next(static_cast<std::size_t>(order + 1), BigInt{0});
        for (std::uint64_t k = 1; k <= n; ++k)
            next[static_cast<std::size_t>(k)] =
                stirling[static_cast<std::size_t>(k - 1)]
                - BigInt::fromUnsigned(n - 1) * stirling[static_cast<std::size_t>(k)];
        stirling = std::move(next);
    }

    std::vector<Expr> terms;
    terms.reserve(static_cast<std::size_t>(order));
    const Expr z = arguments[1];
    for (std::uint64_t k = 1; k <= order; ++k) {
        const BigInt& coefficient = stirling[static_cast<std::size_t>(k)];
        if (coefficient.isZero())
            continue;
        Expr term = polylogOrderShift(arguments[0], k, z, builtins);
        if (!(coefficient == BigInt{1}))
            term = multiply(builtins, {Expr{Number{coefficient}}, std::move(term)});
        terms.push_back(std::move(term));
    }
    Expr ordinary = divide(
        builtins, add(builtins, std::move(terms)),
        power(builtins, z, integer(static_cast<std::int64_t>(order))));

    BigInt factorial{1};
    for (std::uint64_t k = 2; k <= order; ++k)
        factorial *= BigInt::fromUnsigned(k);
    Expr atZero = divide(
        builtins, Expr{Number{factorial}},
        power(builtins, Expr{Number{BigInt::fromUnsigned(order)}}, arguments[0]));
    return detail::makeCases(builtins, {
        detail::makeCaseBranch(builtins, std::move(ordinary),
            call(builtins, BuiltinId::NotEqual, {z, integer(0)})),
        detail::makeCaseBranch(builtins, std::move(atZero),
            call(builtins, BuiltinId::Equal, {z, integer(0)}))});
}

} // namespace

std::optional<Expr> differentiateKnownRepeatedExpression(
    const Expr& expression,
    const expression::Symbol& variable,
    std::uint64_t order,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (auto exponential = repeatedQuadraticExponentialDerivative(
            expression, variable, order, builtins, mathematics, angles))
        return simplify(std::move(*exponential), builtins, mathematics, angles);
    if (auto lambert = repeatedDirectLambertDerivative(
            expression, variable, order, builtins))
        return simplify(std::move(*lambert), builtins, mathematics, angles);
    if (auto polylog = repeatedDirectPolylogDerivative(
            expression, variable, order, builtins))
        return simplify(std::move(*polylog), builtins, mathematics, angles);
    return std::nullopt;
}

Expr canonicalizeDerivativeOutput(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const mathematics::MathRegistry& mathematics,
    const mathematics::AngleSemantics& angles) {
    if (!isAdditiveExpression(expression, builtins))
        return expression;

    std::vector<Expr> terms;
    appendDerivativeLinearTerms(
        expression, Rational{BigInt{1}}, terms, builtins);
    if (terms.empty())
        return integer(0);
    Expr flattened = terms.size() == 1
        ? std::move(terms.front())
        : add(builtins, std::move(terms));
    return simplify(std::move(flattened), builtins, mathematics, angles);
}

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
