// ユーザー定義函数
#include "user_function_registry.hpp"

#include <algorithm>
#include <stdexcept>
#include <utility>

namespace mmcal::evaluation {

UserFunctionDefinition::UserFunctionDefinition(
    expression::Symbol name,
    std::vector<expression::Symbol> parameters,
    expression::Expr body,
    expression::OriginMap origins)
    : name(std::move(name)),
      parameters(std::move(parameters)),
      body(std::move(body)),
      origins(std::move(origins)) {}

std::size_t UserFunctionDefinition::arity() const noexcept {
    return parameters.size();
}

void UserFunctionRegistry::define(UserFunctionDefinition definition) {
    if (!definition.name.valid())
        throw std::invalid_argument("Function name cannot be empty");

    const expression::Symbol name = definition.name;
    const std::size_t arity = definition.arity();
    auto [iterator, inserted] = definitions_.try_emplace(
        name.id(),
        FunctionEntry{name, {}});
    static_cast<void>(inserted);
    iterator->second.overloads.insert_or_assign(arity, std::move(definition));
}

const UserFunctionDefinition* UserFunctionRegistry::find(
    const expression::Symbol& name,
    std::size_t arity) const noexcept {
    const auto function = definitions_.find(name.id());
    if (function == definitions_.end())
        return nullptr;

    const auto overload = function->second.overloads.find(arity);
    return overload == function->second.overloads.end() ? nullptr : &overload->second;
}

bool UserFunctionRegistry::contains(const expression::Symbol& name) const noexcept {
    return definitions_.contains(name.id());
}

std::vector<std::size_t> UserFunctionRegistry::arities(const expression::Symbol& name) const {
    std::vector<std::size_t> result;
    const auto function = definitions_.find(name.id());
    if (function == definitions_.end())
        return result;

    result.reserve(function->second.overloads.size());
    for (const auto& [arity, definition] : function->second.overloads) {
        static_cast<void>(definition);
        result.push_back(arity);
    }
    return result;
}

std::set<std::string> UserFunctionRegistry::names() const {
    std::set<std::string> result;
    for (const auto& [id, entry] : definitions_) {
        static_cast<void>(id);
        result.insert(entry.name.name());
    }
    return result;
}

std::vector<UserFunctionDefinition> UserFunctionRegistry::definitions() const {
    std::vector<UserFunctionDefinition> result;
    result.reserve(size());
    for (const auto& [id, entry] : definitions_) {
        static_cast<void>(id);
        for (const auto& [arity, definition] : entry.overloads) {
            static_cast<void>(arity);
            result.push_back(definition);
        }
    }

    std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
        if (lhs.name.name() != rhs.name.name())
            return lhs.name.name() < rhs.name.name();
        return lhs.arity() < rhs.arity();
    });
    return result;
}

bool UserFunctionRegistry::erase(const expression::Symbol& name) {
    return definitions_.erase(name.id()) != 0;
}

void UserFunctionRegistry::clear() noexcept {
    definitions_.clear();
}

std::size_t UserFunctionRegistry::size() const noexcept {
    std::size_t result = 0;
    for (const auto& [id, entry] : definitions_) {
        static_cast<void>(id);
        result += entry.overloads.size();
    }
    return result;
}

} // namespace mmcal::evaluation
