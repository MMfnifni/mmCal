// 代数計算
#include "exact_algebra.hpp"
#include "expression/exact_value.hpp"

#include "numeric/big_int.hpp"
#include "numeric/number.hpp"

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::mathematics {
namespace {

using expression::Expr;
using numeric::BigInt;
using numeric::Number;
using numeric::Rational;

[[nodiscard]] Rational rational(std::int64_t value) {
    return Rational{BigInt{value}};
}

[[nodiscard]] Expr numberExpr(const BigInt& value) {
    return Expr{Number{value}};
}


[[nodiscard]] bool isHead(
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins,
    evaluation::BuiltinId id) {
    return expression.isCall()
        && expression.asCall().head.sameIdentity(builtins.symbol(id));
}

[[nodiscard]] Expr buildScaledAtom(
    Rational coefficient,
    Expr atom,
    const evaluation::BuiltinRegistry& builtins) {
    if (coefficient.isZero())
        return numberExpr(BigInt{0});

    const bool negative = coefficient.numerator().isNegative();
    const BigInt numeratorMagnitude = coefficient.numerator().abs();
    const BigInt denominator = coefficient.denominator();

    Expr result = std::move(atom);

    // |分子|が1を超える係数は、符号も数値係数そのものへ持たせる。
    // -2*sqrt[2] を -(2*sqrt[2]) と保持すると表示に不要な括弧が生じるため、exact coefficientのcanonical formを -2 sqrt[2] とする。
    // |分子|==1 のときだけ従来どおり単項Negateを使う。
    bool needsOuterNegate = negative;
    if (!(numeratorMagnitude == BigInt{1})) {
        const BigInt signedNumerator = negative ? -numeratorMagnitude : numeratorMagnitude;
        result = Expr::call(
            builtins.symbol(evaluation::BuiltinId::Multiply),
            {numberExpr(signedNumerator), std::move(result)});
        needsOuterNegate = false;
    }

    // |分子|==1の負係数は、分母の外側ではなく分子側へ符号を置く。-(x/2) より -x/2 をcanonical表示にでき、Divide側の有理係数整理とも循環しない。
    if (needsOuterNegate) {
        result = Expr::call(
            builtins.symbol(evaluation::BuiltinId::Negate),
            {std::move(result)});
    }

    // 分母1なら除算ノード自体を作らない。
    if (!(denominator == BigInt{1})) {
        result = Expr::call(
            builtins.symbol(evaluation::BuiltinId::Divide),
            {std::move(result), numberExpr(denominator)});
    }

    return result;
}

} // namespace

Expr scaleExactExpression(
    const Rational& coefficient,
    const Expr& expression,
    const evaluation::BuiltinRegistry& builtins) {
    if (coefficient.isZero())
        return numberExpr(BigInt{0});
    if (coefficient == rational(1))
        return expression;

    // 最も頻出する「有理係数 × 既に分母を持つ記号式」はここで一度に約分する。
    // 例:
    //   2 * ((sqrt[6] - sqrt[2]) / 4)
    //     coefficient = 2
    //     denominator = 4
    //     combined    = 1/2
    //   -> (sqrt[6] - sqrt[2]) / 2
    //
    // この処理は式の数学的内容を推測していない。分母が厳密実有理数である場合だけRationalの四則演算を使うため、近似誤差や不正な恒等変形は発生しない。
    if (isHead(expression, builtins, evaluation::BuiltinId::Divide)) {
        const auto& arguments = expression.asCall().arguments;
        if (arguments.size() == 2) {
            if (const auto denominator = expression::exact::realRational(arguments[1]); denominator && !denominator->isZero())
                return buildScaledAtom(coefficient / *denominator, arguments[0], builtins);
        }
    }

    return buildScaledAtom(coefficient, expression, builtins);
}

std::vector<Expr> combineStructurallyIdenticalTerms(
    std::span<const Expr> terms,
    const evaluation::BuiltinRegistry& builtins) {
    struct Group final {
        Expr expression;
        std::size_t count = 0;
    };

    std::vector<Group> groups;
    groups.reserve(terms.size());

    // 現段階ではhash canonical formをまだ持っていないため線形探索にしている。
    // 加算の項数は通常小さく、何よりoperator==による「構造完全一致」だけを同類項として扱うことを優先する。
    // 式hash導入後はこの部分だけ置換できる。
    for (const Expr& term : terms) {
        auto iterator = groups.begin();
        for (; iterator != groups.end(); ++iterator) {
            if (iterator->expression == term)
                break;
        }

        if (iterator == groups.end()) {
            groups.push_back(Group{term, 1});
            continue;
        }

        ++iterator->count;
    }

    std::vector<Expr> result;
    result.reserve(groups.size());
    for (Group& group : groups) {
        if (group.count == 1) {
            result.push_back(std::move(group.expression));
            continue;
        }

        // size_tからBigIntへ直接のconstructorを持たないため、文字列を経由する。
        // これは「同一項が何個あるか」という構文上の小さな整数にだけ使う。
        const BigInt count = BigInt::parse(std::to_string(group.count));
        result.push_back(scaleExactExpression(
            Rational{count},
            group.expression,
            builtins));
    }

    return result;
}

} // namespace mmcal::mathematics
