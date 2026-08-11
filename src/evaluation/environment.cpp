// ユーザー変数環境
#include "environment.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace mmcal::evaluation {

void Environment::set(expression::Symbol symbol, expression::Expr value) {
    const symbols::SymbolId id = symbol.id();
    bindings_.insert_or_assign(id, Binding{std::move(symbol), std::move(value)});
}

void Environment::assign(expression::Symbol symbol, expression::Expr value) {
    for (auto scope = localScopes_.rbegin(); scope != localScopes_.rend(); ++scope) {
        const auto binding = scope->find(symbol.id());
        if (binding == scope->end())
            continue;

        binding->second.value = std::move(value);
        return;
    }

    set(std::move(symbol), std::move(value));
}

void Environment::setLocal(expression::Symbol symbol, expression::Expr value) {
    if (localScopes_.empty())
        throw std::logic_error("No local scope is active");

    const symbols::SymbolId id = symbol.id();
    localScopes_.back().insert_or_assign(id, Binding{std::move(symbol), std::move(value)});
}

const expression::Expr* Environment::find(const expression::Symbol& symbol) const noexcept {
    for (auto scope = localScopes_.rbegin(); scope != localScopes_.rend(); ++scope) {
        const auto binding = scope->find(symbol.id());
        if (binding != scope->end())
            return &binding->second.value;
    }

    const auto binding = bindings_.find(symbol.id());
    return binding == bindings_.end() ? nullptr : &binding->second.value;
}

bool Environment::contains(const expression::Symbol& symbol) const noexcept {
    return find(symbol) != nullptr;
}

bool Environment::containsLocal(const expression::Symbol& symbol) const noexcept {
    for (auto scope = localScopes_.rbegin(); scope != localScopes_.rend(); ++scope)
        if (scope->contains(symbol.id()))
            return true;

    return false;
}

bool Environment::erase(const expression::Symbol& symbol) {
    return bindings_.erase(symbol.id()) != 0;
}

std::vector<std::pair<expression::Symbol, expression::Expr>> Environment::definitions() const {
    std::vector<std::pair<expression::Symbol, expression::Expr>> result;
    result.reserve(bindings_.size());
    for (const auto& [id, binding] : bindings_) {
        static_cast<void>(id);
        result.emplace_back(binding.symbol, binding.value);
    }

    std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first.name() < rhs.first.name();
    });
    return result;
}

void Environment::pushScope() {
    localScopes_.emplace_back();
}

void Environment::popScope() {
    if (localScopes_.empty())
        throw std::logic_error("No local scope is active");

    localScopes_.pop_back();
}

std::size_t Environment::localDepth() const noexcept {
    return localScopes_.size();
}

void Environment::clear() noexcept {
    bindings_.clear();
    localScopes_.clear();
}

std::size_t Environment::size() const noexcept {
    return bindings_.size();
}

} // namespace mmcal::evaluation
