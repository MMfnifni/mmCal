#pragma once

#include "symbol.hpp"
#include "numeric/decimal_approximation.hpp"
#include "numeric/complex_decimal_approximation.hpp"
#include "numeric/number.hpp"

#include <cstddef>
#include <memory>
#include <span>
#include <string>
#include <utility>
#include <vector>

namespace mmcal::solver { class SolutionSet; }
namespace mmcal::symbolic { class AlgebraicNumber; }

namespace mmcal::expression {

class Expr;
struct ArrayExpr;
struct ListExpr;
struct CallExpr;
class ArrayBuilder;

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

enum class ArrayStorageKind {
    Integer,
    Rational,
    Number,
    DecimalApproximation,
    ComplexDecimalApproximation,
    Generic
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
    [[nodiscard]] static Expr array(ArrayExpr array);
    [[nodiscard]] static Expr integerArray(
        std::vector<std::size_t> shape,
        std::vector<numeric::BigInt> elements);
    [[nodiscard]] static Expr rationalArray(
        std::vector<std::size_t> shape,
        std::vector<numeric::Rational> elements);
    [[nodiscard]] static Expr numberArray(
        std::vector<std::size_t> shape,
        std::vector<numeric::Number> elements);
    [[nodiscard]] static Expr decimalArray(
        std::vector<std::size_t> shape,
        std::vector<numeric::DecimalApproximation> elements);
    [[nodiscard]] static Expr complexDecimalArray(
        std::vector<std::size_t> shape,
        std::vector<numeric::ComplexDecimalApproximation> elements);
    [[nodiscard]] static Expr list(std::vector<Expr> elements);
    [[nodiscard]] static Expr call(
        Symbol head,
        std::vector<Expr> arguments,
        std::shared_ptr<const symbolic::AlgebraicNumber> algebraicValue = {});
    // 同じCallを子だけ再構築する。引数が構造的に不変ならRootの内部算術cacheも保持する。
    [[nodiscard]] static Expr rebuildCall(
        const CallExpr& source,
        std::vector<Expr> arguments);

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
    template <ExprKind Kind, class Value>
    struct TypedNode;

    std::shared_ptr<const Node> node_;

    explicit Expr(std::shared_ptr<const Node> node);
};

struct ArrayExpressionEntry final {
    std::size_t index = 0;
    Expr expression;
};

// dense Arrayの永続表現。値は固定ページ単位でpacked保持し，shape/layoutは不変backingを共有する。
struct ArrayExpr final {
    std::vector<std::size_t> shape;

    ArrayExpr(std::vector<std::size_t> shape, std::vector<Expr> elements);
    ArrayExpr(std::vector<std::size_t> shape, std::vector<numeric::BigInt> elements);
    ArrayExpr(std::vector<std::size_t> shape, std::vector<numeric::Rational> elements);
    ArrayExpr(std::vector<std::size_t> shape, std::vector<numeric::Number> elements);
    ArrayExpr(std::vector<std::size_t> shape, std::vector<numeric::DecimalApproximation> elements);
    ArrayExpr(
        std::vector<std::size_t> shape,
        std::vector<numeric::ComplexDecimalApproximation> elements);

    [[nodiscard]] std::size_t rank() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t extent(std::size_t dimension) const;
    [[nodiscard]] std::size_t flatIndex(std::span<const std::size_t> indices) const;
    [[nodiscard]] bool isVector() const noexcept;
    [[nodiscard]] bool isMatrix() const noexcept;
    [[nodiscard]] bool isContiguous() const noexcept;

    // storageKind()はArray全体の最小共通domain。storedKindAt()は実際のpage表現を返す。
    [[nodiscard]] ArrayStorageKind storageKind() const noexcept;
    [[nodiscard]] ArrayStorageKind storedKindAt(std::size_t index) const;
    [[nodiscard]] bool hasExactNumberStorage() const noexcept;
    [[nodiscard]] bool hasExactRealStorage() const noexcept;
    [[nodiscard]] bool hasStoredExpressions() const noexcept;

    [[nodiscard]] numeric::Number exactNumber(std::size_t index) const;
    [[nodiscard]] Expr element(std::size_t index) const;
    [[nodiscard]] const numeric::BigInt& integerAt(std::size_t index) const;
    [[nodiscard]] const numeric::Rational& rationalAt(std::size_t index) const;
    [[nodiscard]] const numeric::Number& numberAt(std::size_t index) const;
    [[nodiscard]] const numeric::DecimalApproximation& decimalAt(std::size_t index) const;
    [[nodiscard]] const numeric::ComplexDecimalApproximation& complexDecimalAt(
        std::size_t index) const;
    [[nodiscard]] const Expr& expressionAt(std::size_t index) const;

    [[nodiscard]] std::vector<ArrayExpressionEntry> expressionEntries() const;
    [[nodiscard]] std::vector<Expr> storedExpressions() const;
    [[nodiscard]] std::vector<Expr> materialize() const;
    void appendMaterialized(std::vector<Expr>& output) const;

    [[nodiscard]] ArrayExpr reshaped(std::vector<std::size_t> newShape) const;
    [[nodiscard]] ArrayExpr sliced(
        std::vector<std::size_t> newShape,
        std::size_t logicalOffset,
        std::size_t count) const;
    [[nodiscard]] ArrayExpr transposed() const;
    [[nodiscard]] ArrayExpr replacedExpressions(
        std::span<const std::size_t> indices,
        std::vector<Expr> values) const;

    [[nodiscard]] bool operator==(const ArrayExpr& rhs) const;

private:
    struct Storage;

    std::shared_ptr<const Storage> storage_;
    std::vector<std::size_t> strides_;
    std::size_t offset_ = 0;
    std::size_t elementCount_ = 0;
    ArrayStorageKind storageKind_ = ArrayStorageKind::Generic;
    bool hasStoredExpressions_ = false;

    ArrayExpr(
        std::vector<std::size_t> shape,
        std::shared_ptr<const Storage> storage,
        std::vector<std::size_t> strides,
        std::size_t offset,
        std::size_t elementCount,
        ArrayStorageKind storageKind,
        bool hasStoredExpressions);

    [[nodiscard]] std::size_t physicalIndex(std::size_t logicalIndex) const;

    friend class ArrayBuilder;
};

// Array構築時のpromotionを固定ページ内へ限定する。末尾でsymbolic値が現れても既存packed pageは再構築しない。
class ArrayBuilder final {
public:
    ArrayBuilder();
    ~ArrayBuilder();
    ArrayBuilder(ArrayBuilder&&) noexcept;
    ArrayBuilder& operator=(ArrayBuilder&&) noexcept;
    ArrayBuilder(const ArrayBuilder&) = delete;
    ArrayBuilder& operator=(const ArrayBuilder&) = delete;

    void reserve(std::size_t elementCount);
    void append(Expr value);
    void append(numeric::BigInt value);
    void append(numeric::Rational value);
    void append(numeric::Number value);
    void append(numeric::DecimalApproximation value);
    void append(numeric::ComplexDecimalApproximation value);
    void appendArray(const ArrayExpr& array);

    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] ArrayExpr finish(std::vector<std::size_t> shape);

private:
    struct State;
    std::unique_ptr<State> state_;
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
    // Rootのcanonical構造には含めない内部算術cache。
    std::shared_ptr<const symbolic::AlgebraicNumber> algebraicValue;

    [[nodiscard]] bool operator==(const CallExpr& rhs) const {
        return head == rhs.head && arguments == rhs.arguments;
    }
};

} // namespace mmcal::expression
