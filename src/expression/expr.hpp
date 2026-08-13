#pragma once

#include "symbol.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/number.hpp"

#include <cstddef>
#include <memory>
#include <span>
#include <string>
#include <vector>

namespace mmcal::solver { class SolutionSet; }

namespace mmcal::expression {

class Expr;
struct ArrayExpr;
struct ListExpr;
struct CallExpr;

enum class ExprKind {
    Number,
    DecimalApproximation,
    ComplexDecimalApproximation,
    Boolean,
    String,
    Symbol,
    Array,
    List,
    Call,
    SolutionSet
};

// 数値、真偽値、文字列、記号、配列、函数呼び出しを同じ形で保持する不変な式ノード。
class Expr final {
public:
    Expr(numeric::Number number);
    Expr(numeric::DecimalApproximation decimalApproximation);
    Expr(numeric::ComplexDecimalApproximation complexDecimalApproximation);
    Expr(bool boolean);
    Expr(std::string string);
    Expr(Symbol symbol);

    [[nodiscard]] static Expr solutionSet(solver::SolutionSet value);

    [[nodiscard]] static Expr array(
        std::vector<std::size_t> shape,
        std::vector<Expr> elements);
    [[nodiscard]] static Expr list(std::vector<Expr> elements);
    [[nodiscard]] static Expr call(Symbol head, std::vector<Expr> arguments);

    [[nodiscard]] ExprKind kind() const noexcept;
    [[nodiscard]] bool isNumber() const noexcept;
    [[nodiscard]] bool isDecimalApproximation() const noexcept;
    [[nodiscard]] bool isComplexDecimalApproximation() const noexcept;
    [[nodiscard]] bool isBoolean() const noexcept;
    [[nodiscard]] bool isString() const noexcept;
    [[nodiscard]] bool isSymbol() const noexcept;
    [[nodiscard]] bool isArray() const noexcept;
    [[nodiscard]] bool isList() const noexcept;
    [[nodiscard]] bool isCall() const noexcept;
    [[nodiscard]] bool isSolutionSet() const noexcept;

    [[nodiscard]] const numeric::Number& asNumber() const;
    [[nodiscard]] const numeric::DecimalApproximation& asDecimalApproximation() const;
    [[nodiscard]] const numeric::ComplexDecimalApproximation& asComplexDecimalApproximation() const;
    [[nodiscard]] bool asBoolean() const;
    [[nodiscard]] const std::string& asString() const;
    [[nodiscard]] const Symbol& asSymbol() const;
    [[nodiscard]] const ArrayExpr& asArray() const;
    [[nodiscard]] const ListExpr& asList() const;
    [[nodiscard]] const CallExpr& asCall() const;
    [[nodiscard]] const solver::SolutionSet& asSolutionSet() const;

    // 出自情報などを外部表へ関連付けるための、式ノード固有の識別子。
    [[nodiscard]] const void* identity() const noexcept;
    [[nodiscard]] bool operator==(const Expr& rhs) const;

private:
    struct Node;
    std::shared_ptr<const Node> node_;

    explicit Expr(std::shared_ptr<const Node> node);
};

struct ArrayExpr final {
    std::vector<std::size_t> shape;
    std::vector<Expr> elements;

    [[nodiscard]] std::size_t rank() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::size_t extent(std::size_t dimension) const;
    [[nodiscard]] std::size_t flatIndex(std::span<const std::size_t> indices) const;
    [[nodiscard]] bool isVector() const noexcept;
    [[nodiscard]] bool isMatrix() const noexcept;
    [[nodiscard]] bool operator==(const ArrayExpr&) const = default;
};

// {}構文のうち矩形dense Arrayへ正規化できない一般brace値。
// Q/RやU/S/Vのように互いにshapeが異なる複数Arrayを一つの値として保持する。
// 行列算法はこの型をMatrixとして受理せず，境界で矩形性を監査する。
struct ListExpr final {
    std::vector<Expr> elements;

    [[nodiscard]] std::size_t size() const noexcept { return elements.size(); }
    [[nodiscard]] bool operator==(const ListExpr&) const = default;
};

struct CallExpr final {
    Symbol head;
    std::vector<Expr> arguments;

    [[nodiscard]] bool operator==(const CallExpr&) const = default;
};

} // namespace mmcal::expression
