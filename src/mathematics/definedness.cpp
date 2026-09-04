// 式の定義条件生成
#include "definedness.hpp"

#include "evaluation/builtin_registry.hpp"
#include "knowledge_context.hpp"
#include "math_registry.hpp"
#include "numeric/big_int.hpp"
#include "numeric/number.hpp"
#include "predicate.hpp"

#include <cstdint>
#include <span>
#include <utility>

namespace mmcal::mathematics {
namespace {

using evaluation::BuiltinId;
using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

[[nodiscard]] Expr integer(std::int64_t value) {
    return Expr{Number{BigInt{value}}};
}

class DefinednessCollector final {
public:
    DefinednessCollector(
        const evaluation::BuiltinRegistry& builtins,
        const MathRegistry& mathematics)
        : builtins_(builtins), mathematics_(mathematics), knowledge_(builtins, mathematics, none_) {}

    [[nodiscard]] std::optional<AssumptionSet> collect(const Expr& expression) {
        if (!collectExpression(expression))
            return std::nullopt;
        return conditions_;
    }

private:
    const evaluation::BuiltinRegistry& builtins_;
    const MathRegistry& mathematics_;
    AssumptionSet none_;
    KnowledgeContext knowledge_;
    AssumptionSet conditions_;

    [[nodiscard]] bool collectExpression(const Expr& expression) {
        if (expression.isNumber() || expression.isSymbol())
            return true;
        if (expression.isArray()) {
            for (std::size_t i = 0; i < expression.asArray().size(); ++i)
                if (!collectExpression(expression.asArray().element(i)))
                    return false;
            return true;
        }
        if (expression.isList()) {
            for (const Expr& element : expression.asList().elements)
                if (!collectExpression(element))
                    return false;
            return true;
        }
        if (!expression.isCall())
            return false;

        const auto* builtin = builtins_.find(expression.asCall().head);
        if (!builtin)
            return false;
        const auto& arguments = expression.asCall().arguments;

        switch (builtin->id) {
        case BuiltinId::Add:
        case BuiltinId::Subtract:
        case BuiltinId::Multiply:
        case BuiltinId::Negate:
            return collectAll(arguments);

        case BuiltinId::Divide:
            if (arguments.size() != 2 || !collectAll(arguments))
                return false;
            requireNonZero(arguments[1]);
            return true;

        default:
            break;
        }

        const FunctionDefinition* function = mathematics_.findFunction(expression.asCall().head);
        if (!function || !function->acceptsArity(arguments.size()))
            return false;
        return collectFunction(*function, arguments);
    }

    [[nodiscard]] bool collectAll(std::span<const Expr> arguments) {
        for (const Expr& argument : arguments)
            if (!collectExpression(argument))
                return false;
        return true;
    }

    void requireNonZero(Expr expression) {
        // product/powerの非零条件は，可能なら因子・底へexactに分解する。
        // x(x+1)!=0 や u^2!=0 をそのまま保持するより，x!=0 && x+1!=0，u!=0
        // の方がKnowledgeContextで再利用しやすく，数学的にも同値である。
        if (expression.isCall()) {
            const auto* builtin = builtins_.find(expression.asCall().head);
            const auto& arguments = expression.asCall().arguments;
            if (builtin && builtin->id == BuiltinId::Multiply) {
                for (const Expr& factor : arguments)
                    requireNonZero(factor);
                return;
            }
            if (builtin && builtin->id == BuiltinId::Negate && arguments.size() == 1) {
                requireNonZero(arguments[0]);
                return;
            }
            if (builtin && builtin->id == BuiltinId::Power && arguments.size() == 2
                && arguments[1].isNumber() && arguments[1].asNumber().isReal()
                && arguments[1].asNumber().asReal().isInteger()
                && !arguments[1].asNumber().asReal().asInteger().isZero()) {
                requireNonZero(arguments[0]);
                return;
            }
        }

        Predicate predicate = relation(RelationKind::NotEqual, std::move(expression), integer(0));
        if (knowledge_.prove(predicate) != TruthValue::True)
            conditions_.add(std::move(predicate));
    }

    void requireDomain(Expr expression, NumericDomain domain) {
        Predicate predicate = elementOf(std::move(expression), domain);
        if (knowledge_.prove(predicate) != TruthValue::True)
            conditions_.add(std::move(predicate));
    }

    [[nodiscard]] Expr square(const Expr& expression) const {
        return Expr::call(
            builtins_.symbol(BuiltinId::Power),
            {expression, integer(2)});
    }

    [[nodiscard]] bool collectFunction(
        const FunctionDefinition& function,
        std::span<const Expr> arguments) {
        if (function.definednessRule == FunctionDefinednessRule::PrincipalPower)
            return collectPrincipalPower(arguments);

        // tan[atan[z]] / tanh[atanh[z]] はinverseが有限に定義される点では外側函数も必ず定義される。
        // genericなCosNonZero/CoshNonZeroを重ねると，恒等式Solveに冗長な条件だけが残る。
        if (arguments.size() == 1 && arguments[0].isCall()) {
            const FunctionDefinition* inner = mathematics_.findFunction(arguments[0].asCall().head);
            if (inner && ((function.id == FunctionId::Tan && inner->id == FunctionId::Atan)
                || (function.id == FunctionId::Tanh && inner->id == FunctionId::Atanh)))
                return collectExpression(arguments[0]);
        }

        if (!collectAll(arguments))
            return false;

        switch (function.definednessRule) {
        case FunctionDefinednessRule::Everywhere:
            return true;

        case FunctionDefinednessRule::ArgumentReal:
            requireDomain(arguments[0], NumericDomain::Real);
            return true;

        case FunctionDefinednessRule::ArgumentsReal:
            for (const Expr& argument : arguments)
                requireDomain(argument, NumericDomain::Real);
            return true;

        case FunctionDefinednessRule::ArgumentNonZero:
            requireNonZero(arguments[0]);
            return true;

        case FunctionDefinednessRule::Logarithm:
            if (arguments.size() == 1) {
                requireNonZero(arguments[0]);
                return true;
            }
            if (arguments.size() == 2) {
                requireNonZero(arguments[0]);
                requireNonZero(Expr::call(
                    builtins_.symbol(BuiltinId::Subtract),
                    {arguments[0], integer(1)}));
                requireNonZero(arguments[1]);
                return true;
            }
            return false;

        case FunctionDefinednessRule::OnePlusArgumentNonZero:
            requireNonZero(Expr::call(
                builtins_.symbol(BuiltinId::Add),
                {integer(1), arguments[0]}));
            return true;

        case FunctionDefinednessRule::SinNonZero:
            requireNonZero(Expr::call(builtins_.symbol(BuiltinId::Sin), {arguments[0]}));
            return true;

        case FunctionDefinednessRule::CosNonZero:
            requireNonZero(Expr::call(builtins_.symbol(BuiltinId::Cos), {arguments[0]}));
            return true;

        case FunctionDefinednessRule::SinhNonZero:
            requireNonZero(Expr::call(builtins_.symbol(BuiltinId::Sinh), {arguments[0]}));
            return true;

        case FunctionDefinednessRule::CoshNonZero:
            requireNonZero(Expr::call(builtins_.symbol(BuiltinId::Cosh), {arguments[0]}));
            return true;

        case FunctionDefinednessRule::OnePlusSquareNonZero:
            requireNonZero(Expr::call(
                builtins_.symbol(BuiltinId::Add),
                {integer(1), square(arguments[0])}));
            return true;

        case FunctionDefinednessRule::OneMinusSquareNonZero:
            requireNonZero(Expr::call(
                builtins_.symbol(BuiltinId::Subtract),
                {integer(1), square(arguments[0])}));
            return true;

        case FunctionDefinednessRule::ArgumentsPositiveReal:
            for (const Expr& argument : arguments) {
                requireDomain(argument, NumericDomain::Real);
                Predicate positive = relation(RelationKind::Greater, argument, integer(0));
                if (knowledge_.prove(positive) != TruthValue::True)
                    conditions_.add(std::move(positive));
            }
            return true;

        case FunctionDefinednessRule::GammaPoles:
            // {0,-1,-2,...} の補集合は一般symbolic Predicateでは有限個に表せない。
            // ただしexact numeric引数ならpoleか否かを完全に判定できるので，
            // 定数項の相殺まで不必要に阻害しない。
            if (arguments.size() != 1 || !arguments[0].isNumber())
                return false;
            if (!arguments[0].asNumber().isReal())
                return true;
            if (!arguments[0].asNumber().asReal().isInteger())
                return true;
            return arguments[0].asNumber().asReal().asInteger().isPositive();

        case FunctionDefinednessRule::ZetaPole: {
            // Riemann zetaの有限平面上の唯一のpoleはs=1。branch cutを持たないため，
            // symbolicでもs!=1をexactなdefinedness条件として共有できる。
            Predicate predicate = relation(RelationKind::NotEqual, arguments[0], integer(1));
            if (knowledge_.prove(predicate) != TruthValue::True)
                conditions_.add(std::move(predicate));
            return true;
        }

        case FunctionDefinednessRule::IncompleteBetaPrincipal:
            // 第一版ibetaはa,b>0かつ0<=x<=1のregularized real branchだけを実装する。
            // この三条件をgeneric definedness collectorへ半端に複製せず専用builtin側で検証する。
            return false;

        case FunctionDefinednessRule::Hypergeometric1F1Poles:
            // 1F1(a;b;z) の b=0,-1,-2,... pole集合も有限Predicateでは正確に表せない。
            // terminating seriesによる可除ケースもあるため、弱い条件を捏造せず未解決とする。
            return false;

        case FunctionDefinednessRule::Hypergeometric2F1Poles:
            // 2F1(a,b;c;z) も c=0,-1,-2,... にparameter poleを持つ。
            // principal branch cutとterminating caseを同時に有限Predicateへ落とさない。
            return false;

        case FunctionDefinednessRule::EllipticPrincipal:
            // Legendre楕円積分のbranch/singularity条件はparameterとamplitudeの双方に依存する。
            // 一般形をEverywhereと誤認するより、現段階ではdefinedness証明を保守的に保留する。
            return false;

        case FunctionDefinednessRule::SpecialPrincipal:
            // li/polylog等のprincipal branchはbranch cutとparameter依存特異点を持つ。
            // 有限Predicateへ弱く近似せず、証明不能として保守的に扱う。
            return false;

        case FunctionDefinednessRule::RealPairNotBothZero:
        case FunctionDefinednessRule::PrincipalPower:
            return false;
        }
        return false;
    }

    [[nodiscard]] bool collectPrincipalPower(std::span<const Expr> arguments) {
        if (arguments.size() != 2 || !collectExpression(arguments[0]))
            return false;
        if (!arguments[1].isNumber() || !arguments[1].asNumber().isReal())
            return false;

        // exact real exponentはRationalとして扱える。principal powerでは0^qはq>0なら
        // 定義され，q<=0ならbase!=0を要する。一般複素指数はここでは証明しない。
        const Rational exponent = arguments[1].asNumber().asReal().toRational();
        if (exponent.isZero() || exponent.numerator().isNegative())
            requireNonZero(arguments[0]);
        return true;
    }
};

} // namespace

std::optional<AssumptionSet> expressionDomainConditions(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    const MathRegistry& mathematics) {
    return DefinednessCollector{builtins, mathematics}.collect(expression);
}

} // namespace mmcal::mathematics
