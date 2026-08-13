// 式ASTと値表現
#include "expr.hpp"

#include "solver/solution_set.hpp"
#include "array_utils.hpp"

#include <limits>
#include <stdexcept>
#include <utility>
#include <variant>

namespace mmcal::expression {

struct Expr::Node final {
    using Value = std::variant<
        numeric::Number,
        numeric::DecimalApproximation,
        numeric::ComplexDecimalApproximation,
        bool,
        std::string,
        Symbol,
        ArrayExpr,
        ListExpr,
        CallExpr,
        std::shared_ptr<const solver::SolutionSet>>;

    explicit Node(Value value)
        : value(std::move(value)) {}

    Value value;
};

Expr::Expr(numeric::Number number)
    : node_(std::make_shared<Node>(Node::Value{std::move(number)})) {}

Expr::Expr(numeric::DecimalApproximation decimalApproximation)
    : node_(std::make_shared<Node>(Node::Value{std::move(decimalApproximation)})) {}

Expr::Expr(numeric::ComplexDecimalApproximation complexDecimalApproximation)
    : node_(std::make_shared<Node>(Node::Value{std::move(complexDecimalApproximation)})) {}

Expr::Expr(bool boolean)
    : node_(std::make_shared<Node>(Node::Value{boolean})) {}

Expr::Expr(std::string string)
    : node_(std::make_shared<Node>(Node::Value{std::move(string)})) {}

Expr::Expr(Symbol symbol)
    : node_(std::make_shared<Node>(Node::Value{std::move(symbol)})) {}

Expr::Expr(std::shared_ptr<const Node> node)
    : node_(std::move(node)) {}

Expr Expr::solutionSet(solver::SolutionSet value) {
    return Expr{std::make_shared<Node>(Node::Value{
        std::make_shared<const solver::SolutionSet>(std::move(value))})};
}

Expr Expr::array(
    std::vector<std::size_t> shape,
    std::vector<Expr> elements) {
    if (arrayElementCount(shape) != elements.size())
        throw std::invalid_argument("Array shape does not match the element count");
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
                if (flattenedCount > std::numeric_limits<std::size_t>::max() - child.elements.size())
                    throw std::length_error("Array element count exceeds the size_t range");
                flattenedCount += child.elements.size();
            }

            shape.insert(shape.end(), childShape.begin(), childShape.end());
            std::vector<Expr> flattened;
            flattened.reserve(flattenedCount);
            for (const Expr& element : elements) {
                const ArrayExpr& child = element.asArray();
                flattened.insert(flattened.end(), child.elements.begin(), child.elements.end());
            }
            elements = std::move(flattened);
        }
    }

    if (arrayElementCount(shape) != elements.size())
        throw std::invalid_argument("Normalized array shape does not match the element count");
    ArrayExpr array{std::move(shape), std::move(elements)};
    return Expr{std::make_shared<Node>(Node::Value{std::move(array)})};
}

Expr Expr::list(std::vector<Expr> elements) {
    ListExpr list{std::move(elements)};
    return Expr{std::make_shared<Node>(Node::Value{std::move(list)})};
}

Expr Expr::call(Symbol head, std::vector<Expr> arguments) {
    CallExpr call{std::move(head), std::move(arguments)};
    return Expr{std::make_shared<Node>(Node::Value{std::move(call)})};
}

ExprKind Expr::kind() const noexcept {
    return static_cast<ExprKind>(node_->value.index());
}

bool Expr::isNumber() const noexcept {
    return kind() == ExprKind::Number;
}

bool Expr::isDecimalApproximation() const noexcept {
    return kind() == ExprKind::DecimalApproximation;
}

bool Expr::isComplexDecimalApproximation() const noexcept {
    return kind() == ExprKind::ComplexDecimalApproximation;
}

bool Expr::isBoolean() const noexcept {
    return kind() == ExprKind::Boolean;
}

bool Expr::isString() const noexcept {
    return kind() == ExprKind::String;
}

bool Expr::isSymbol() const noexcept {
    return kind() == ExprKind::Symbol;
}

bool Expr::isArray() const noexcept {
    return kind() == ExprKind::Array;
}

bool Expr::isList() const noexcept {
    return kind() == ExprKind::List;
}

bool Expr::isCall() const noexcept {
    return kind() == ExprKind::Call;
}

bool Expr::isSolutionSet() const noexcept {
    return kind() == ExprKind::SolutionSet;
}

const numeric::Number& Expr::asNumber() const {
    if (!isNumber())
        throw std::logic_error("Expr does not contain a number");

    return std::get<numeric::Number>(node_->value);
}

const numeric::DecimalApproximation& Expr::asDecimalApproximation() const {
    if (!isDecimalApproximation())
        throw std::logic_error("Expr does not contain a decimal approximation");

    return std::get<numeric::DecimalApproximation>(node_->value);
}

const numeric::ComplexDecimalApproximation& Expr::asComplexDecimalApproximation() const {
    if (!isComplexDecimalApproximation())
        throw std::logic_error("Expr does not contain a complex decimal approximation");

    return std::get<numeric::ComplexDecimalApproximation>(node_->value);
}

bool Expr::asBoolean() const {
    if (!isBoolean())
        throw std::logic_error("Expr does not contain a boolean");

    return std::get<bool>(node_->value);
}

const std::string& Expr::asString() const {
    if (!isString())
        throw std::logic_error("Expr does not contain a string");

    return std::get<std::string>(node_->value);
}

const Symbol& Expr::asSymbol() const {
    if (!isSymbol())
        throw std::logic_error("Expr does not contain a symbol");

    return std::get<Symbol>(node_->value);
}

const ArrayExpr& Expr::asArray() const {
    if (!isArray())
        throw std::logic_error("Expr does not contain an array");

    return std::get<ArrayExpr>(node_->value);
}

const ListExpr& Expr::asList() const {
    if (!isList())
        throw std::logic_error("Expr does not contain a list");

    return std::get<ListExpr>(node_->value);
}

const CallExpr& Expr::asCall() const {
    if (!isCall())
        throw std::logic_error("Expr does not contain a function call");

    return std::get<CallExpr>(node_->value);
}

const solver::SolutionSet& Expr::asSolutionSet() const {
    if (!isSolutionSet())
        throw std::logic_error("Expr does not contain a solution set");

    return *std::get<std::shared_ptr<const solver::SolutionSet>>(node_->value);
}

const void* Expr::identity() const noexcept {
    return node_.get();
}

bool Expr::operator==(const Expr& rhs) const {
    if (node_ == rhs.node_)
        return true;
    if (kind() != rhs.kind())
        return false;
    if (isSolutionSet())
        return asSolutionSet() == rhs.asSolutionSet();
    return node_->value == rhs.node_->value;
}

std::size_t ArrayExpr::rank() const noexcept {
    return shape.size();
}

std::size_t ArrayExpr::size() const noexcept {
    return elements.size();
}

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

bool ArrayExpr::isVector() const noexcept {
    return rank() == 1;
}

bool ArrayExpr::isMatrix() const noexcept {
    return rank() == 2;
}

} // namespace mmcal::expression
