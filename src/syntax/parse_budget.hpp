#pragma once

#include "error/error_message.hpp"
#include "source/source_span.hpp"

#include <cstddef>
#include <string>
#include <utility>

namespace mmcal::syntax {

// Lexer/ParserからLowererへ危険な深さ・幅の構文木を渡さないための共通上限。
// 上限超過は数学的DomainErrorではなく，入力資源制限として報告する。
struct ParseLimits final {
    std::size_t maxTokens = 1'000'000;
    std::size_t maxNodes = 250'000;
    std::size_t maxRecursionDepth = 1024;
    std::size_t maxOperatorChain = 256;
    std::size_t maxLiteralDigits = 100'000;
    std::size_t maxCallArguments = 100'000;
    std::size_t maxArrayElements = 250'000;
};

class ParseBudget final {
public:
    class DepthGuard final {
    public:
        DepthGuard(const DepthGuard&) = delete;
        DepthGuard& operator=(const DepthGuard&) = delete;

        DepthGuard(DepthGuard&& other) noexcept
            : owner_(std::exchange(other.owner_, nullptr)) {}

        DepthGuard& operator=(DepthGuard&&) = delete;

        ~DepthGuard() {
            if (owner_)
                owner_->leaveDepth();
        }

    private:
        friend class ParseBudget;

        explicit DepthGuard(ParseBudget& owner) noexcept
            : owner_(&owner) {}

        ParseBudget* owner_ = nullptr;
    };

    explicit ParseBudget(ParseLimits limits = {}) noexcept
        : limits_(limits) {}

    [[nodiscard]] const ParseLimits& limits() const noexcept {
        return limits_;
    }

    [[nodiscard]] DepthGuard enter(source::SourceSpan span) {
        if (depth_ >= limits_.maxRecursionDepth)
            exceeded("Parser nesting depth", limits_.maxRecursionDepth, span);
        ++depth_;
        return DepthGuard{*this};
    }

    void consumeNode(source::SourceSpan span) {
        if (nodeCount_ >= limits_.maxNodes)
            exceeded("Parser AST node count", limits_.maxNodes, span);
        ++nodeCount_;
    }

    void checkTokenCount(std::size_t count, source::SourceSpan span) const {
        if (count > limits_.maxTokens)
            exceeded("Lexer token count", limits_.maxTokens, span);
    }

    void checkOperatorChain(std::size_t count, source::SourceSpan span) const {
        if (count > limits_.maxOperatorChain)
            exceeded("Parser operator-chain length", limits_.maxOperatorChain, span);
    }

    void checkLiteralDigits(std::size_t count, source::SourceSpan span) const {
        if (count > limits_.maxLiteralDigits)
            exceeded("Numeric literal digit count", limits_.maxLiteralDigits, span);
    }

    void checkCallArguments(std::size_t count, source::SourceSpan span) const {
        if (count > limits_.maxCallArguments)
            exceeded("Function argument count", limits_.maxCallArguments, span);
    }

    void checkArrayElements(std::size_t count, source::SourceSpan span) const {
        if (count > limits_.maxArrayElements)
            exceeded("Array element count", limits_.maxArrayElements, span);
    }

    [[nodiscard]] std::size_t nodeCount() const noexcept {
        return nodeCount_;
    }

private:
    ParseLimits limits_;
    std::size_t nodeCount_ = 0;
    std::size_t depth_ = 0;

    void leaveDepth() noexcept {
        if (depth_ != 0)
            --depth_;
    }

    [[noreturn]] static void exceeded(
        std::string name,
        std::size_t limit,
        source::SourceSpan span) {
        error::throwCalcError(
            error::CalcErrorType::ResourceLimit,
            std::move(name) + " exceeded (limit " + std::to_string(limit) + ')',
            span);
    }
};

} // namespace mmcal::syntax
