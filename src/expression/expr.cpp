// 式ASTと値表現
#include "expr.hpp"

#include "solver/solution_set.hpp"
#include "array_utils.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <variant>

namespace mmcal::expression {
namespace {

constexpr std::size_t arrayPageCapacity = 1024;

using PageStorage = std::variant<
    std::vector<numeric::BigInt>,
    std::vector<numeric::Rational>,
    std::vector<numeric::Number>,
    std::vector<numeric::DecimalApproximation>,
    std::vector<numeric::ComplexDecimalApproximation>,
    std::vector<Expr>>;

struct ArrayPage final {
    explicit ArrayPage(PageStorage pageStorage)
        : storage(std::move(pageStorage)) {}

    PageStorage storage;

    [[nodiscard]] std::size_t size() const noexcept {
        return std::visit([](const auto& values) { return values.size(); }, storage);
    }

    [[nodiscard]] ArrayStorageKind kind() const noexcept {
        if (std::holds_alternative<std::vector<numeric::BigInt>>(storage))
            return ArrayStorageKind::Integer;
        if (std::holds_alternative<std::vector<numeric::Rational>>(storage))
            return ArrayStorageKind::Rational;
        if (std::holds_alternative<std::vector<numeric::Number>>(storage))
            return ArrayStorageKind::Number;
        if (std::holds_alternative<std::vector<numeric::DecimalApproximation>>(storage))
            return ArrayStorageKind::DecimalApproximation;
        if (std::holds_alternative<std::vector<numeric::ComplexDecimalApproximation>>(storage))
            return ArrayStorageKind::ComplexDecimalApproximation;
        return ArrayStorageKind::Generic;
    }
};

[[nodiscard]] ArrayStorageKind combinedStorageKind(
    ArrayStorageKind lhs,
    ArrayStorageKind rhs) noexcept {
    if (lhs == rhs)
        return lhs;

    const auto exactRank = [](ArrayStorageKind kind) -> int {
        switch (kind) {
        case ArrayStorageKind::Integer: return 0;
        case ArrayStorageKind::Rational: return 1;
        case ArrayStorageKind::Number: return 2;
        default: return -1;
        }
    };

    const int lhsRank = exactRank(lhs);
    const int rhsRank = exactRank(rhs);
    if (lhsRank >= 0 && rhsRank >= 0) {
        switch (std::max(lhsRank, rhsRank)) {
        case 0: return ArrayStorageKind::Integer;
        case 1: return ArrayStorageKind::Rational;
        default: return ArrayStorageKind::Number;
        }
    }
    return ArrayStorageKind::Generic;
}

[[nodiscard]] std::vector<std::size_t> contiguousStrides(
    std::span<const std::size_t> shape) {
    std::vector<std::size_t> strides(shape.size(), 1);
    std::size_t stride = 1;
    for (std::size_t i = shape.size(); i-- > 0;) {
        strides[i] = stride;
        if (shape[i] != 0)
            stride *= shape[i];
    }
    return strides;
}

[[nodiscard]] std::shared_ptr<const ArrayPage> canonicalExprPage(
    std::vector<Expr> values) {
    if (values.empty())
        return std::make_shared<const ArrayPage>(PageStorage{std::move(values)});

    const bool allNumbers = std::all_of(values.begin(), values.end(),
        [](const Expr& value) { return value.isNumber(); });
    if (allNumbers) {
        const bool allReal = std::all_of(values.begin(), values.end(),
            [](const Expr& value) { return value.asNumber().isReal(); });
        if (allReal) {
            const bool allInteger = std::all_of(values.begin(), values.end(),
                [](const Expr& value) { return value.asNumber().asReal().isInteger(); });
            if (allInteger) {
                std::vector<numeric::BigInt> packed;
                packed.reserve(values.size());
                for (const Expr& value : values)
                    packed.push_back(value.asNumber().asReal().asInteger());
                return std::make_shared<const ArrayPage>(PageStorage{std::move(packed)});
            }

            std::vector<numeric::Rational> packed;
            packed.reserve(values.size());
            for (const Expr& value : values)
                packed.push_back(value.asNumber().asReal().toRational());
            return std::make_shared<const ArrayPage>(PageStorage{std::move(packed)});
        }

        std::vector<numeric::Number> packed;
        packed.reserve(values.size());
        for (const Expr& value : values)
            packed.push_back(value.asNumber());
        return std::make_shared<const ArrayPage>(PageStorage{std::move(packed)});
    }

    if (std::all_of(values.begin(), values.end(),
        [](const Expr& value) { return value.isDecimalApproximation(); })) {
        std::vector<numeric::DecimalApproximation> packed;
        packed.reserve(values.size());
        for (const Expr& value : values)
            packed.push_back(value.asDecimalApproximation());
        return std::make_shared<const ArrayPage>(PageStorage{std::move(packed)});
    }

    if (std::all_of(values.begin(), values.end(),
        [](const Expr& value) { return value.isComplexDecimalApproximation(); })) {
        std::vector<numeric::ComplexDecimalApproximation> packed;
        packed.reserve(values.size());
        for (const Expr& value : values)
            packed.push_back(value.asComplexDecimalApproximation());
        return std::make_shared<const ArrayPage>(PageStorage{std::move(packed)});
    }

    return std::make_shared<const ArrayPage>(PageStorage{std::move(values)});
}

} // namespace

struct ArrayExpr::Storage final {
    std::vector<std::shared_ptr<const ArrayPage>> pages;
    std::size_t size = 0;
};

struct ArrayBuilder::State final {
    using CurrentStorage = std::variant<
        std::monostate,
        std::vector<numeric::BigInt>,
        std::vector<numeric::Rational>,
        std::vector<numeric::Number>,
        std::vector<numeric::DecimalApproximation>,
        std::vector<numeric::ComplexDecimalApproximation>,
        std::vector<Expr>>;

    std::vector<std::shared_ptr<const ArrayPage>> pages;
    CurrentStorage current;
    std::size_t size = 0;
    bool hasKind = false;
    ArrayStorageKind storageKind = ArrayStorageKind::Generic;
    bool hasStoredExpressions = false;

    [[nodiscard]] std::size_t currentSize() const noexcept {
        return std::visit([](const auto& values) -> std::size_t {
            using T = std::decay_t<decltype(values)>;
            if constexpr (std::is_same_v<T, std::monostate>)
                return 0;
            else
                return values.size();
        }, current);
    }

    [[nodiscard]] ArrayStorageKind currentKind() const noexcept {
        if (std::holds_alternative<std::vector<numeric::BigInt>>(current))
            return ArrayStorageKind::Integer;
        if (std::holds_alternative<std::vector<numeric::Rational>>(current))
            return ArrayStorageKind::Rational;
        if (std::holds_alternative<std::vector<numeric::Number>>(current))
            return ArrayStorageKind::Number;
        if (std::holds_alternative<std::vector<numeric::DecimalApproximation>>(current))
            return ArrayStorageKind::DecimalApproximation;
        if (std::holds_alternative<std::vector<numeric::ComplexDecimalApproximation>>(current))
            return ArrayStorageKind::ComplexDecimalApproximation;
        if (std::holds_alternative<std::vector<Expr>>(current))
            return ArrayStorageKind::Generic;
        return ArrayStorageKind::Generic;
    }

    void recordKind(ArrayStorageKind kind) noexcept {
        storageKind = hasKind ? combinedStorageKind(storageKind, kind) : kind;
        hasKind = true;
        if (kind == ArrayStorageKind::Generic)
            hasStoredExpressions = true;
    }

    void flush() {
        if (currentSize() == 0) {
            current = std::monostate{};
            return;
        }

        std::shared_ptr<const ArrayPage> page = std::visit([](auto&& values)
            -> std::shared_ptr<const ArrayPage> {
            using T = std::decay_t<decltype(values)>;
            if constexpr (std::is_same_v<T, std::monostate>)
                throw std::logic_error("ArrayBuilder has no current page");
            else {
                using Vector = T;
                return std::make_shared<const ArrayPage>(PageStorage{Vector{std::move(values)}});
            }
        }, std::move(current));

        recordKind(page->kind());
        pages.push_back(std::move(page));
        current = std::monostate{};
    }

    template <class Vector>
    void start(Vector values = {}) {
        values.reserve(arrayPageCapacity);
        current = std::move(values);
    }

    void flushIfFull() {
        if (currentSize() == arrayPageCapacity)
            flush();
    }

    void promoteExact(ArrayStorageKind target) {
        const ArrayStorageKind currentStorageKind = currentKind();
        if (currentStorageKind == target)
            return;

        if (currentStorageKind == ArrayStorageKind::Integer
            && target == ArrayStorageKind::Rational) {
            auto integers = std::move(std::get<std::vector<numeric::BigInt>>(current));
            std::vector<numeric::Rational> rationals;
            rationals.reserve(arrayPageCapacity);
            for (numeric::BigInt& value : integers)
                rationals.emplace_back(std::move(value));
            current = std::move(rationals);
            return;
        }

        if (currentStorageKind == ArrayStorageKind::Integer
            && target == ArrayStorageKind::Number) {
            auto integers = std::move(std::get<std::vector<numeric::BigInt>>(current));
            std::vector<numeric::Number> numbers;
            numbers.reserve(arrayPageCapacity);
            for (numeric::BigInt& value : integers)
                numbers.emplace_back(std::move(value));
            current = std::move(numbers);
            return;
        }

        if (currentStorageKind == ArrayStorageKind::Rational
            && target == ArrayStorageKind::Number) {
            auto rationals = std::move(std::get<std::vector<numeric::Rational>>(current));
            std::vector<numeric::Number> numbers;
            numbers.reserve(arrayPageCapacity);
            for (numeric::Rational& value : rationals)
                numbers.emplace_back(std::move(value));
            current = std::move(numbers);
            return;
        }

        throw std::logic_error("Invalid exact Array page promotion");
    }

    void promoteCurrentToGeneric() {
        if (currentSize() == 0) {
            start(std::vector<Expr>{});
            return;
        }
        if (currentKind() == ArrayStorageKind::Generic)
            return;

        std::vector<Expr> expressions = std::visit([](auto&& values) {
            using T = std::decay_t<decltype(values)>;
            std::vector<Expr> result;
            if constexpr (!std::is_same_v<T, std::monostate>) {
                result.reserve(arrayPageCapacity);
                if constexpr (std::is_same_v<T, std::vector<numeric::BigInt>>) {
                    for (auto& value : values)
                        result.emplace_back(numeric::Number{std::move(value)});
                }
                else if constexpr (std::is_same_v<T, std::vector<numeric::Rational>>) {
                    for (auto& value : values)
                        result.emplace_back(numeric::Number{std::move(value)});
                }
                else if constexpr (std::is_same_v<T, std::vector<numeric::Number>>) {
                    for (auto& value : values)
                        result.emplace_back(std::move(value));
                }
                else if constexpr (std::is_same_v<T, std::vector<numeric::DecimalApproximation>>) {
                    for (auto& value : values)
                        result.emplace_back(std::move(value));
                }
                else if constexpr (std::is_same_v<T, std::vector<numeric::ComplexDecimalApproximation>>) {
                    for (auto& value : values)
                        result.emplace_back(std::move(value));
                }
                else if constexpr (std::is_same_v<T, std::vector<Expr>>) {
                    result = std::move(values);
                }
            }
            return result;
        }, std::move(current));
        current = std::move(expressions);
    }

    void appendInteger(numeric::BigInt value) {
        if (currentSize() == 0)
            start(std::vector<numeric::BigInt>{});

        switch (currentKind()) {
        case ArrayStorageKind::Integer:
            std::get<std::vector<numeric::BigInt>>(current).push_back(std::move(value));
            break;
        case ArrayStorageKind::Rational:
            std::get<std::vector<numeric::Rational>>(current).emplace_back(std::move(value));
            break;
        case ArrayStorageKind::Number:
            std::get<std::vector<numeric::Number>>(current).emplace_back(std::move(value));
            break;
        case ArrayStorageKind::Generic:
            std::get<std::vector<Expr>>(current).emplace_back(numeric::Number{std::move(value)});
            break;
        default:
            promoteCurrentToGeneric();
            std::get<std::vector<Expr>>(current).emplace_back(numeric::Number{std::move(value)});
            break;
        }
        ++size;
        flushIfFull();
    }

    void appendRational(numeric::Rational value) {
        if (value.isInteger()) {
            appendInteger(value.numerator());
            return;
        }
        if (currentSize() == 0)
            start(std::vector<numeric::Rational>{});

        switch (currentKind()) {
        case ArrayStorageKind::Integer:
            promoteExact(ArrayStorageKind::Rational);
            [[fallthrough]];
        case ArrayStorageKind::Rational:
            std::get<std::vector<numeric::Rational>>(current).push_back(std::move(value));
            break;
        case ArrayStorageKind::Number:
            std::get<std::vector<numeric::Number>>(current).emplace_back(std::move(value));
            break;
        case ArrayStorageKind::Generic:
            std::get<std::vector<Expr>>(current).emplace_back(numeric::Number{std::move(value)});
            break;
        default:
            promoteCurrentToGeneric();
            std::get<std::vector<Expr>>(current).emplace_back(numeric::Number{std::move(value)});
            break;
        }
        ++size;
        flushIfFull();
    }

    void appendNumber(numeric::Number value) {
        if (value.isReal()) {
            if (value.asReal().isInteger()) {
                appendInteger(value.asReal().asInteger());
                return;
            }
            appendRational(value.asReal().toRational());
            return;
        }
        if (currentSize() == 0)
            start(std::vector<numeric::Number>{});

        switch (currentKind()) {
        case ArrayStorageKind::Integer:
        case ArrayStorageKind::Rational:
            promoteExact(ArrayStorageKind::Number);
            [[fallthrough]];
        case ArrayStorageKind::Number:
            std::get<std::vector<numeric::Number>>(current).push_back(std::move(value));
            break;
        case ArrayStorageKind::Generic:
            std::get<std::vector<Expr>>(current).emplace_back(std::move(value));
            break;
        default:
            promoteCurrentToGeneric();
            std::get<std::vector<Expr>>(current).emplace_back(std::move(value));
            break;
        }
        ++size;
        flushIfFull();
    }

    void appendDecimal(numeric::DecimalApproximation value) {
        if (currentSize() == 0)
            start(std::vector<numeric::DecimalApproximation>{});
        if (currentKind() == ArrayStorageKind::DecimalApproximation)
            std::get<std::vector<numeric::DecimalApproximation>>(current).push_back(std::move(value));
        else if (currentKind() == ArrayStorageKind::Generic)
            std::get<std::vector<Expr>>(current).emplace_back(std::move(value));
        else {
            promoteCurrentToGeneric();
            std::get<std::vector<Expr>>(current).emplace_back(std::move(value));
        }
        ++size;
        flushIfFull();
    }

    void appendComplexDecimal(numeric::ComplexDecimalApproximation value) {
        if (currentSize() == 0)
            start(std::vector<numeric::ComplexDecimalApproximation>{});
        if (currentKind() == ArrayStorageKind::ComplexDecimalApproximation)
            std::get<std::vector<numeric::ComplexDecimalApproximation>>(current).push_back(std::move(value));
        else if (currentKind() == ArrayStorageKind::Generic)
            std::get<std::vector<Expr>>(current).emplace_back(std::move(value));
        else {
            promoteCurrentToGeneric();
            std::get<std::vector<Expr>>(current).emplace_back(std::move(value));
        }
        ++size;
        flushIfFull();
    }

    void appendExpr(Expr value) {
        // Generic pageに入った後はpage境界までExprのまま保持し，交互入力で小pageが乱立するのを避ける。
        if (currentSize() != 0 && currentKind() == ArrayStorageKind::Generic) {
            std::get<std::vector<Expr>>(current).push_back(std::move(value));
            ++size;
            flushIfFull();
            return;
        }

        if (value.isNumber()) {
            appendNumber(value.asNumber());
            return;
        }
        if (value.isDecimalApproximation()) {
            appendDecimal(value.asDecimalApproximation());
            return;
        }
        if (value.isComplexDecimalApproximation()) {
            appendComplexDecimal(value.asComplexDecimalApproximation());
            return;
        }

        promoteCurrentToGeneric();
        std::get<std::vector<Expr>>(current).push_back(std::move(value));
        ++size;
        flushIfFull();
    }

    void appendSharedPage(std::shared_ptr<const ArrayPage> page) {
        if (currentSize() != 0)
            throw std::logic_error("Shared Array page requires an aligned builder");
        if (page->size() != arrayPageCapacity)
            throw std::logic_error("Only complete Array pages may be shared");
        recordKind(page->kind());
        pages.push_back(std::move(page));
        size += arrayPageCapacity;
    }
};

struct Expr::Node {
    explicit Node(ExprKind nodeKind) noexcept
        : kind(nodeKind) {}

    ExprKind kind;
};

template <ExprKind Kind, class Value>
struct Expr::TypedNode final : Node {
    explicit TypedNode(Value nodeValue)
        : Node(Kind), value(std::move(nodeValue)) {}

    Value value;
};

Expr::Expr(numeric::Number number)
    : node_(std::make_shared<TypedNode<ExprKind::Number, numeric::Number>>(
        std::move(number))) {}

Expr::Expr(numeric::DecimalApproximation decimalApproximation)
    : node_(std::make_shared<TypedNode<
        ExprKind::DecimalApproximation, numeric::DecimalApproximation>>(
        std::move(decimalApproximation))) {}

Expr::Expr(numeric::ComplexDecimalApproximation complexDecimalApproximation)
    : node_(std::make_shared<TypedNode<
        ExprKind::ComplexDecimalApproximation, numeric::ComplexDecimalApproximation>>(
        std::move(complexDecimalApproximation))) {}

Expr::Expr(bool boolean)
    : node_(std::make_shared<TypedNode<ExprKind::Boolean, bool>>(boolean)) {}

Expr::Expr(std::string string)
    : node_(std::make_shared<TypedNode<ExprKind::String, std::string>>(
        std::move(string))) {}

Expr::Expr(Symbol symbol)
    : node_(std::make_shared<TypedNode<ExprKind::Symbol, Symbol>>(
        std::move(symbol))) {}

Expr::Expr(std::shared_ptr<const Node> node)
    : node_(std::move(node)) {}

Expr Expr::solutionSet(solver::SolutionSet value) {
    return Expr{std::make_shared<TypedNode<
        ExprKind::SolutionSet, std::shared_ptr<const solver::SolutionSet>>>(
        std::make_shared<const solver::SolutionSet>(std::move(value)))};
}

Expr Expr::array(
    std::vector<std::size_t> shape,
    std::vector<Expr> elements) {
    if (arrayElementCount(shape) != elements.size())
        throw std::invalid_argument("Array shape does not match the element count");

    ArrayBuilder builder;
    if (!elements.empty()) {
        for (const Expr& element : elements)
            if (element.isList())
                throw std::invalid_argument("Dense Array cannot contain a non-rectangular brace value");

        const bool nested = elements.front().isArray();
        for (const Expr& element : elements)
            if (element.isArray() != nested)
                throw std::invalid_argument("Array elements must have a uniform rank and shape");

        if (nested) {
            const std::vector<std::size_t> childShape = elements.front().asArray().shape;
            std::size_t flattenedCount = 0;
            for (const Expr& element : elements) {
                const ArrayExpr& child = element.asArray();
                if (child.shape != childShape)
                    throw std::invalid_argument("Array child shapes must be identical");
                if (flattenedCount > std::numeric_limits<std::size_t>::max() - child.size())
                    throw std::length_error("Array element count exceeds the size_t range");
                flattenedCount += child.size();
            }

            shape.insert(shape.end(), childShape.begin(), childShape.end());
            builder.reserve(flattenedCount);
            for (const Expr& element : elements)
                builder.appendArray(element.asArray());
        }
        else {
            builder.reserve(elements.size());
            for (Expr& element : elements)
                builder.append(std::move(element));
        }
    }

    ArrayExpr array = builder.finish(std::move(shape));
    return Expr{std::make_shared<TypedNode<ExprKind::Array, ArrayExpr>>(
        std::move(array))};
}

Expr Expr::array(ArrayExpr array) {
    if (arrayElementCount(array.shape) != array.size())
        throw std::invalid_argument("Array shape does not match the element count");
    return Expr{std::make_shared<TypedNode<ExprKind::Array, ArrayExpr>>(
        std::move(array))};
}

Expr Expr::integerArray(
    std::vector<std::size_t> shape,
    std::vector<numeric::BigInt> elements) {
    return Expr::array(ArrayExpr{std::move(shape), std::move(elements)});
}

Expr Expr::rationalArray(
    std::vector<std::size_t> shape,
    std::vector<numeric::Rational> elements) {
    return Expr::array(ArrayExpr{std::move(shape), std::move(elements)});
}

Expr Expr::numberArray(
    std::vector<std::size_t> shape,
    std::vector<numeric::Number> elements) {
    return Expr::array(ArrayExpr{std::move(shape), std::move(elements)});
}

Expr Expr::decimalArray(
    std::vector<std::size_t> shape,
    std::vector<numeric::DecimalApproximation> elements) {
    return Expr::array(ArrayExpr{std::move(shape), std::move(elements)});
}

Expr Expr::complexDecimalArray(
    std::vector<std::size_t> shape,
    std::vector<numeric::ComplexDecimalApproximation> elements) {
    return Expr::array(ArrayExpr{std::move(shape), std::move(elements)});
}

Expr Expr::list(std::vector<Expr> elements) {
    ListExpr list{std::move(elements)};
    return Expr{std::make_shared<TypedNode<ExprKind::List, ListExpr>>(
        std::move(list))};
}

Expr Expr::call(
    Symbol head,
    std::vector<Expr> arguments,
    std::shared_ptr<const symbolic::AlgebraicNumber> algebraicValue) {
    CallExpr call{std::move(head), std::move(arguments), std::move(algebraicValue)};
    return Expr{std::make_shared<TypedNode<ExprKind::Call, CallExpr>>(
        std::move(call))};
}

Expr Expr::rebuildCall(
    const CallExpr& source,
    std::vector<Expr> arguments) {
    auto algebraicValue = source.arguments == arguments
        ? source.algebraicValue : nullptr;
    return call(source.head, std::move(arguments), std::move(algebraicValue));
}

ExprKind Expr::kind() const noexcept {
    return node_->kind;
}

bool Expr::isNumber() const noexcept { return kind() == ExprKind::Number; }
bool Expr::isDecimalApproximation() const noexcept { return kind() == ExprKind::DecimalApproximation; }
bool Expr::isComplexDecimalApproximation() const noexcept { return kind() == ExprKind::ComplexDecimalApproximation; }
bool Expr::isBoolean() const noexcept { return kind() == ExprKind::Boolean; }
bool Expr::isString() const noexcept { return kind() == ExprKind::String; }
bool Expr::isSymbol() const noexcept { return kind() == ExprKind::Symbol; }
bool Expr::isArray() const noexcept { return kind() == ExprKind::Array; }
bool Expr::isList() const noexcept { return kind() == ExprKind::List; }
bool Expr::isCall() const noexcept { return kind() == ExprKind::Call; }
bool Expr::isSolutionSet() const noexcept { return kind() == ExprKind::SolutionSet; }

const numeric::Number& Expr::asNumber() const {
    if (!isNumber())
        throw std::logic_error("Expr does not contain a number");
    return static_cast<const TypedNode<ExprKind::Number, numeric::Number>&>(*node_).value;
}

const numeric::DecimalApproximation& Expr::asDecimalApproximation() const {
    if (!isDecimalApproximation())
        throw std::logic_error("Expr does not contain a decimal approximation");
    return static_cast<const TypedNode<
        ExprKind::DecimalApproximation, numeric::DecimalApproximation>&>(*node_).value;
}

const numeric::ComplexDecimalApproximation& Expr::asComplexDecimalApproximation() const {
    if (!isComplexDecimalApproximation())
        throw std::logic_error("Expr does not contain a complex decimal approximation");
    return static_cast<const TypedNode<
        ExprKind::ComplexDecimalApproximation, numeric::ComplexDecimalApproximation>&>(*node_).value;
}

bool Expr::asBoolean() const {
    if (!isBoolean())
        throw std::logic_error("Expr does not contain a boolean");
    return static_cast<const TypedNode<ExprKind::Boolean, bool>&>(*node_).value;
}

const std::string& Expr::asString() const {
    if (!isString())
        throw std::logic_error("Expr does not contain a string");
    return static_cast<const TypedNode<ExprKind::String, std::string>&>(*node_).value;
}

const Symbol& Expr::asSymbol() const {
    if (!isSymbol())
        throw std::logic_error("Expr does not contain a symbol");
    return static_cast<const TypedNode<ExprKind::Symbol, Symbol>&>(*node_).value;
}

const ArrayExpr& Expr::asArray() const {
    if (!isArray())
        throw std::logic_error("Expr does not contain an array");
    return static_cast<const TypedNode<ExprKind::Array, ArrayExpr>&>(*node_).value;
}

const ListExpr& Expr::asList() const {
    if (!isList())
        throw std::logic_error("Expr does not contain a list");
    return static_cast<const TypedNode<ExprKind::List, ListExpr>&>(*node_).value;
}

const CallExpr& Expr::asCall() const {
    if (!isCall())
        throw std::logic_error("Expr does not contain a function call");
    return static_cast<const TypedNode<ExprKind::Call, CallExpr>&>(*node_).value;
}

const solver::SolutionSet& Expr::asSolutionSet() const {
    if (!isSolutionSet())
        throw std::logic_error("Expr does not contain a solution set");
    return *static_cast<const TypedNode<
        ExprKind::SolutionSet, std::shared_ptr<const solver::SolutionSet>>&>(*node_).value;
}

const void* Expr::identity() const noexcept { return node_.get(); }

bool Expr::operator==(const Expr& rhs) const {
    if (node_ == rhs.node_)
        return true;
    if (kind() != rhs.kind())
        return false;

    switch (kind()) {
    case ExprKind::Number: return asNumber() == rhs.asNumber();
    case ExprKind::DecimalApproximation: return asDecimalApproximation() == rhs.asDecimalApproximation();
    case ExprKind::ComplexDecimalApproximation: return asComplexDecimalApproximation() == rhs.asComplexDecimalApproximation();
    case ExprKind::Boolean: return asBoolean() == rhs.asBoolean();
    case ExprKind::String: return asString() == rhs.asString();
    case ExprKind::Symbol: return asSymbol() == rhs.asSymbol();
    case ExprKind::Array: return asArray() == rhs.asArray();
    case ExprKind::List: return asList() == rhs.asList();
    case ExprKind::Call: return asCall() == rhs.asCall();
    case ExprKind::SolutionSet: return asSolutionSet() == rhs.asSolutionSet();
    }
    return false;
}

ArrayBuilder::ArrayBuilder()
    : state_(std::make_unique<State>()) {}

ArrayBuilder::~ArrayBuilder() = default;
ArrayBuilder::ArrayBuilder(ArrayBuilder&&) noexcept = default;
ArrayBuilder& ArrayBuilder::operator=(ArrayBuilder&&) noexcept = default;

void ArrayBuilder::reserve(std::size_t elementCount) {
    state_->pages.reserve((elementCount + arrayPageCapacity - 1) / arrayPageCapacity);
}

void ArrayBuilder::append(Expr value) { state_->appendExpr(std::move(value)); }
void ArrayBuilder::append(numeric::BigInt value) { state_->appendInteger(std::move(value)); }
void ArrayBuilder::append(numeric::Rational value) { state_->appendRational(std::move(value)); }
void ArrayBuilder::append(numeric::Number value) { state_->appendNumber(std::move(value)); }
void ArrayBuilder::append(numeric::DecimalApproximation value) { state_->appendDecimal(std::move(value)); }
void ArrayBuilder::append(numeric::ComplexDecimalApproximation value) { state_->appendComplexDecimal(std::move(value)); }

void ArrayBuilder::appendArray(const ArrayExpr& array) {
    std::size_t logical = 0;
    while (logical < array.size()) {
        if (state_->currentSize() == 0 && array.isContiguous()) {
            const std::size_t physical = array.offset_ + logical;
            const std::size_t remaining = array.size() - logical;
            if (physical % arrayPageCapacity == 0 && remaining >= arrayPageCapacity) {
                const auto& page = array.storage_->pages[physical / arrayPageCapacity];
                if (page->size() == arrayPageCapacity) {
                    state_->appendSharedPage(page);
                    logical += arrayPageCapacity;
                    continue;
                }
            }
        }

        switch (array.storedKindAt(logical)) {
        case ArrayStorageKind::Integer:
            append(array.integerAt(logical));
            break;
        case ArrayStorageKind::Rational:
            append(array.rationalAt(logical));
            break;
        case ArrayStorageKind::Number:
            append(array.numberAt(logical));
            break;
        case ArrayStorageKind::DecimalApproximation:
            append(array.decimalAt(logical));
            break;
        case ArrayStorageKind::ComplexDecimalApproximation:
            append(array.complexDecimalAt(logical));
            break;
        case ArrayStorageKind::Generic:
            append(array.expressionAt(logical));
            break;
        }
        ++logical;
    }
}

std::size_t ArrayBuilder::size() const noexcept { return state_->size; }

ArrayExpr ArrayBuilder::finish(std::vector<std::size_t> shape) {
    if (arrayElementCount(shape) != state_->size)
        throw std::invalid_argument("Array shape does not match the element count");
    state_->flush();

    auto storage = std::make_shared<ArrayExpr::Storage>();
    storage->pages = std::move(state_->pages);
    storage->size = state_->size;

    const ArrayStorageKind kind = state_->hasKind
        ? state_->storageKind
        : ArrayStorageKind::Generic;
    const bool hasExpressions = state_->hasStoredExpressions;
    const std::size_t count = state_->size;
    std::vector<std::size_t> strides = contiguousStrides(shape);
    return ArrayExpr{
        std::move(shape),
        std::move(storage),
        std::move(strides),
        0,
        count,
        kind,
        hasExpressions};
}

ArrayExpr::ArrayExpr(std::vector<std::size_t> arrayShape, std::vector<Expr> elements) {
    ArrayBuilder builder;
    builder.reserve(elements.size());
    for (Expr& element : elements)
        builder.append(std::move(element));
    *this = builder.finish(std::move(arrayShape));
}

ArrayExpr::ArrayExpr(
    std::vector<std::size_t> arrayShape,
    std::vector<numeric::BigInt> elements) {
    ArrayBuilder builder;
    builder.reserve(elements.size());
    for (auto& element : elements)
        builder.append(std::move(element));
    *this = builder.finish(std::move(arrayShape));
}

ArrayExpr::ArrayExpr(
    std::vector<std::size_t> arrayShape,
    std::vector<numeric::Rational> elements) {
    ArrayBuilder builder;
    builder.reserve(elements.size());
    for (auto& element : elements)
        builder.append(std::move(element));
    *this = builder.finish(std::move(arrayShape));
}

ArrayExpr::ArrayExpr(
    std::vector<std::size_t> arrayShape,
    std::vector<numeric::Number> elements) {
    ArrayBuilder builder;
    builder.reserve(elements.size());
    for (auto& element : elements)
        builder.append(std::move(element));
    *this = builder.finish(std::move(arrayShape));
}

ArrayExpr::ArrayExpr(
    std::vector<std::size_t> arrayShape,
    std::vector<numeric::DecimalApproximation> elements) {
    ArrayBuilder builder;
    builder.reserve(elements.size());
    for (auto& element : elements)
        builder.append(std::move(element));
    *this = builder.finish(std::move(arrayShape));
}

ArrayExpr::ArrayExpr(
    std::vector<std::size_t> arrayShape,
    std::vector<numeric::ComplexDecimalApproximation> elements) {
    ArrayBuilder builder;
    builder.reserve(elements.size());
    for (auto& element : elements)
        builder.append(std::move(element));
    *this = builder.finish(std::move(arrayShape));
}

ArrayExpr::ArrayExpr(
    std::vector<std::size_t> arrayShape,
    std::shared_ptr<const Storage> backingStorage,
    std::vector<std::size_t> arrayStrides,
    std::size_t offset,
    std::size_t elementCount,
    ArrayStorageKind storageKind,
    bool hasStoredExpressions)
    : shape(std::move(arrayShape)),
      storage_(std::move(backingStorage)),
      strides_(std::move(arrayStrides)),
      offset_(offset),
      elementCount_(elementCount),
      storageKind_(storageKind),
      hasStoredExpressions_(hasStoredExpressions) {}

std::size_t ArrayExpr::rank() const noexcept { return shape.size(); }
std::size_t ArrayExpr::size() const noexcept { return elementCount_; }
bool ArrayExpr::empty() const noexcept { return elementCount_ == 0; }

std::size_t ArrayExpr::extent(std::size_t dimension) const {
    if (dimension >= shape.size())
        throw std::out_of_range("Array dimension is out of range");
    return shape[dimension];
}

std::size_t ArrayExpr::flatIndex(std::span<const std::size_t> indices) const {
    if (indices.size() != shape.size())
        throw std::invalid_argument("Array index rank does not match the array rank");

    std::size_t index = 0;
    for (std::size_t dimension = 0; dimension < shape.size(); ++dimension) {
        if (indices[dimension] >= shape[dimension])
            throw std::out_of_range("Array index is out of range");
        index = index * shape[dimension] + indices[dimension];
    }
    return index;
}

bool ArrayExpr::isVector() const noexcept { return rank() == 1; }
bool ArrayExpr::isMatrix() const noexcept { return rank() == 2; }

bool ArrayExpr::isContiguous() const noexcept {
    if (empty())
        return true;
    std::size_t expected = 1;
    for (std::size_t i = shape.size(); i-- > 0;) {
        if (shape[i] > 1 && strides_[i] != expected)
            return false;
        expected *= shape[i];
    }
    return true;
}

ArrayStorageKind ArrayExpr::storageKind() const noexcept { return storageKind_; }
bool ArrayExpr::hasExactNumberStorage() const noexcept {
    return storageKind_ == ArrayStorageKind::Integer
        || storageKind_ == ArrayStorageKind::Rational
        || storageKind_ == ArrayStorageKind::Number;
}

bool ArrayExpr::hasExactRealStorage() const noexcept {
    return storageKind_ == ArrayStorageKind::Integer
        || storageKind_ == ArrayStorageKind::Rational;
}

bool ArrayExpr::hasStoredExpressions() const noexcept { return hasStoredExpressions_; }

std::size_t ArrayExpr::physicalIndex(std::size_t logicalIndex) const {
    if (logicalIndex >= elementCount_)
        throw std::out_of_range("Array flat index is out of range");
    if (isContiguous())
        return offset_ + logicalIndex;

    if (shape.size() == 2) {
        const std::size_t columns = shape[1];
        const std::size_t row = columns == 0 ? 0 : logicalIndex / columns;
        const std::size_t column = columns == 0 ? 0 : logicalIndex % columns;
        return offset_ + row * strides_[0] + column * strides_[1];
    }

    std::size_t remaining = logicalIndex;
    std::size_t physical = offset_;
    for (std::size_t i = shape.size(); i-- > 0;) {
        const std::size_t coordinate = shape[i] == 0 ? 0 : remaining % shape[i];
        if (shape[i] != 0)
            remaining /= shape[i];
        physical += coordinate * strides_[i];
    }
    return physical;
}

ArrayStorageKind ArrayExpr::storedKindAt(std::size_t index) const {
    const std::size_t physical = physicalIndex(index);
    return storage_->pages[physical / arrayPageCapacity]->kind();
}

const numeric::BigInt& ArrayExpr::integerAt(std::size_t index) const {
    const std::size_t physical = physicalIndex(index);
    const auto& page = storage_->pages[physical / arrayPageCapacity];
    if (page->kind() != ArrayStorageKind::Integer)
        throw std::logic_error("Array element is not stored as an integer");
    return std::get<std::vector<numeric::BigInt>>(page->storage)[physical % arrayPageCapacity];
}

const numeric::Rational& ArrayExpr::rationalAt(std::size_t index) const {
    const std::size_t physical = physicalIndex(index);
    const auto& page = storage_->pages[physical / arrayPageCapacity];
    if (page->kind() != ArrayStorageKind::Rational)
        throw std::logic_error("Array element is not stored as a rational");
    return std::get<std::vector<numeric::Rational>>(page->storage)[physical % arrayPageCapacity];
}

const numeric::Number& ArrayExpr::numberAt(std::size_t index) const {
    const std::size_t physical = physicalIndex(index);
    const auto& page = storage_->pages[physical / arrayPageCapacity];
    if (page->kind() != ArrayStorageKind::Number)
        throw std::logic_error("Array element is not stored as a Number");
    return std::get<std::vector<numeric::Number>>(page->storage)[physical % arrayPageCapacity];
}

const numeric::DecimalApproximation& ArrayExpr::decimalAt(std::size_t index) const {
    const std::size_t physical = physicalIndex(index);
    const auto& page = storage_->pages[physical / arrayPageCapacity];
    if (page->kind() != ArrayStorageKind::DecimalApproximation)
        throw std::logic_error("Array element is not stored as a decimal approximation");
    return std::get<std::vector<numeric::DecimalApproximation>>(page->storage)[physical % arrayPageCapacity];
}

const numeric::ComplexDecimalApproximation& ArrayExpr::complexDecimalAt(
    std::size_t index) const {
    const std::size_t physical = physicalIndex(index);
    const auto& page = storage_->pages[physical / arrayPageCapacity];
    if (page->kind() != ArrayStorageKind::ComplexDecimalApproximation)
        throw std::logic_error("Array element is not stored as a complex decimal approximation");
    return std::get<std::vector<numeric::ComplexDecimalApproximation>>(page->storage)
        [physical % arrayPageCapacity];
}

const Expr& ArrayExpr::expressionAt(std::size_t index) const {
    const std::size_t physical = physicalIndex(index);
    const auto& page = storage_->pages[physical / arrayPageCapacity];
    if (page->kind() != ArrayStorageKind::Generic)
        throw std::logic_error("Array element is not stored as an expression");
    return std::get<std::vector<Expr>>(page->storage)[physical % arrayPageCapacity];
}

numeric::Number ArrayExpr::exactNumber(std::size_t index) const {
    switch (storedKindAt(index)) {
    case ArrayStorageKind::Integer:
        return numeric::Number{integerAt(index)};
    case ArrayStorageKind::Rational:
        return numeric::Number{rationalAt(index)};
    case ArrayStorageKind::Number:
        return numberAt(index);
    default:
        throw std::logic_error("Array element is not an exact Number");
    }
}

Expr ArrayExpr::element(std::size_t index) const {
    switch (storedKindAt(index)) {
    case ArrayStorageKind::Integer:
        return Expr{numeric::Number{integerAt(index)}};
    case ArrayStorageKind::Rational:
        return Expr{numeric::Number{rationalAt(index)}};
    case ArrayStorageKind::Number:
        return Expr{numberAt(index)};
    case ArrayStorageKind::DecimalApproximation:
        return Expr{decimalAt(index)};
    case ArrayStorageKind::ComplexDecimalApproximation:
        return Expr{complexDecimalAt(index)};
    case ArrayStorageKind::Generic:
        return expressionAt(index);
    }
    throw std::logic_error("Unknown Array storage kind");
}

std::vector<ArrayExpressionEntry> ArrayExpr::expressionEntries() const {
    std::vector<ArrayExpressionEntry> result;
    if (!hasStoredExpressions_)
        return result;

    for (std::size_t i = 0; i < size(); ++i)
        if (storedKindAt(i) == ArrayStorageKind::Generic)
            result.push_back(ArrayExpressionEntry{i, expressionAt(i)});
    return result;
}

std::vector<Expr> ArrayExpr::storedExpressions() const {
    std::vector<Expr> result;
    if (!hasStoredExpressions_)
        return result;
    for (std::size_t i = 0; i < size(); ++i)
        if (storedKindAt(i) == ArrayStorageKind::Generic)
            result.push_back(expressionAt(i));
    return result;
}

std::vector<Expr> ArrayExpr::materialize() const {
    std::vector<Expr> output;
    output.reserve(size());
    appendMaterialized(output);
    return output;
}

void ArrayExpr::appendMaterialized(std::vector<Expr>& output) const {
    for (std::size_t i = 0; i < size(); ++i)
        output.push_back(element(i));
}

ArrayExpr ArrayExpr::reshaped(std::vector<std::size_t> newShape) const {
    if (arrayElementCount(newShape) != size())
        throw std::invalid_argument("Array shape does not match the element count");
    if (isContiguous()) {
        std::vector<std::size_t> newStrides = contiguousStrides(newShape);
        return ArrayExpr{
            std::move(newShape), storage_, std::move(newStrides), offset_, size(),
            storageKind_, hasStoredExpressions_};
    }

    ArrayBuilder builder;
    builder.reserve(size());
    builder.appendArray(*this);
    return builder.finish(std::move(newShape));
}

ArrayExpr ArrayExpr::sliced(
    std::vector<std::size_t> newShape,
    std::size_t logicalOffset,
    std::size_t count) const {
    if (logicalOffset > size() || count > size() - logicalOffset)
        throw std::out_of_range("Array slice is out of range");
    if (arrayElementCount(newShape) != count)
        throw std::invalid_argument("Array slice shape does not match the element count");
    if (newShape.size() > shape.size())
        throw std::invalid_argument("Array slice rank exceeds the source rank");

    const std::size_t newOffset = count == 0 ? offset_ : physicalIndex(logicalOffset);
    const std::size_t suffix = shape.size() - newShape.size();
    std::vector<std::size_t> newStrides(
        strides_.begin() + static_cast<std::ptrdiff_t>(suffix), strides_.end());
    return ArrayExpr{
        std::move(newShape), storage_, std::move(newStrides), newOffset, count,
        storageKind_, hasStoredExpressions_};
}

ArrayExpr ArrayExpr::transposed() const {
    if (rank() != 2)
        throw std::invalid_argument("Array transpose requires rank 2");
    std::vector<std::size_t> newShape{shape[1], shape[0]};
    std::vector<std::size_t> newStrides{strides_[1], strides_[0]};
    return ArrayExpr{
        std::move(newShape), storage_, std::move(newStrides), offset_, size(),
        storageKind_, hasStoredExpressions_};
}

ArrayExpr ArrayExpr::replacedExpressions(
    std::span<const std::size_t> indices,
    std::vector<Expr> values) const {
    if (indices.size() != values.size())
        throw std::invalid_argument("Array replacement count does not match");
    if (indices.empty())
        return *this;

    auto storage = std::make_shared<Storage>(*storage_);
    std::vector<std::size_t> affectedPages;
    affectedPages.reserve(indices.size());
    for (const std::size_t index : indices) {
        const std::size_t physical = physicalIndex(index);
        const std::size_t pageIndex = physical / arrayPageCapacity;
        if (storage->pages[pageIndex]->kind() != ArrayStorageKind::Generic)
            throw std::logic_error("Array replacement targets a packed numeric element");
        if (std::find(affectedPages.begin(), affectedPages.end(), pageIndex) == affectedPages.end())
            affectedPages.push_back(pageIndex);
    }

    for (const std::size_t pageIndex : affectedPages) {
        const auto& sourcePage = storage->pages[pageIndex];
        std::vector<Expr> pageValues = std::get<std::vector<Expr>>(sourcePage->storage);
        for (std::size_t i = 0; i < indices.size(); ++i) {
            const std::size_t physical = physicalIndex(indices[i]);
            if (physical / arrayPageCapacity == pageIndex)
                pageValues[physical % arrayPageCapacity] = std::move(values[i]);
        }
        storage->pages[pageIndex] = canonicalExprPage(std::move(pageValues));
    }

    bool hasKind = false;
    ArrayStorageKind kind = ArrayStorageKind::Generic;
    bool hasExpressions = false;
    if (isContiguous()) {
        std::size_t physical = offset_;
        const std::size_t end = offset_ + size();
        while (physical < end) {
            const auto& page = storage->pages[physical / arrayPageCapacity];
            kind = hasKind ? combinedStorageKind(kind, page->kind()) : page->kind();
            hasKind = true;
            hasExpressions = hasExpressions || page->kind() == ArrayStorageKind::Generic;
            const std::size_t nextPage = (physical / arrayPageCapacity + 1) * arrayPageCapacity;
            physical = std::min(end, nextPage);
        }
    }
    else {
        for (std::size_t i = 0; i < size(); ++i) {
            const std::size_t physical = physicalIndex(i);
            const auto& page = storage->pages[physical / arrayPageCapacity];
            kind = hasKind ? combinedStorageKind(kind, page->kind()) : page->kind();
            hasKind = true;
            hasExpressions = hasExpressions || page->kind() == ArrayStorageKind::Generic;
        }
    }

    return ArrayExpr{
        shape, std::move(storage), strides_, offset_, size(),
        hasKind ? kind : ArrayStorageKind::Generic, hasExpressions};
}

bool ArrayExpr::operator==(const ArrayExpr& rhs) const {
    if (shape != rhs.shape || size() != rhs.size())
        return false;
    if (storage_ == rhs.storage_ && offset_ == rhs.offset_ && strides_ == rhs.strides_)
        return true;
    if (hasExactNumberStorage() && rhs.hasExactNumberStorage()) {
        for (std::size_t i = 0; i < size(); ++i)
            if (!(exactNumber(i) == rhs.exactNumber(i)))
                return false;
        return true;
    }
    for (std::size_t i = 0; i < size(); ++i)
        if (!(element(i) == rhs.element(i)))
            return false;
    return true;
}

} // namespace mmcal::expression
