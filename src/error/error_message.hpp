#pragma once

#include "source/source_document.hpp"
#include "source/source_span.hpp"

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace mmcal::error {

enum class CalcErrorType {
    Syntax,
    Domain,
    Type,
    Overflow,
    Name,
    Evaluation,
    Internal
};

struct ErrorTraceFrame final {
    std::string label;
    source::SourceReference source;
};

class CalcError final : public std::runtime_error {
public:
    CalcError(CalcErrorType type, std::string message);
    CalcError(CalcErrorType type, std::string message, source::SourceSpan span);
    CalcError(
        CalcErrorType type,
        std::string message,
        source::SourceReference source,
        std::string sourceLabel = "At");

    [[nodiscard]] CalcErrorType type() const noexcept;
    [[nodiscard]] const std::optional<source::SourceSpan>& span() const noexcept;
    [[nodiscard]] const std::shared_ptr<const source::SourceDocument>& document() const noexcept;
    [[nodiscard]] const std::string& sourceLabel() const noexcept;
    [[nodiscard]] const std::vector<ErrorTraceFrame>& trace() const noexcept;

    void attachSourceIfMissing(
        source::SourceReference source,
        std::string_view label = "At");
    void attachDocumentIfMissing(std::shared_ptr<const source::SourceDocument> document);
    void setSourceLabel(std::string label);
    void addTrace(std::string label, source::SourceReference source);

private:
    CalcErrorType type_;
    std::optional<source::SourceSpan> span_;
    std::shared_ptr<const source::SourceDocument> document_;
    std::string sourceLabel_ = "At";
    std::vector<ErrorTraceFrame> trace_;
};

[[nodiscard]] std::string_view calcErrorTypeName(CalcErrorType type) noexcept;
[[nodiscard]] std::string errorMessage(
    const CalcError& error,
    std::string_view sourceText = {});

[[noreturn]] void throwCalcError(CalcErrorType type, std::string message);
[[noreturn]] void throwCalcError(
    CalcErrorType type,
    std::string message,
    source::SourceSpan span);
[[noreturn]] void throwCalcError(
    CalcErrorType type,
    std::string message,
    source::SourceReference source,
    std::string sourceLabel = "At");

} // namespace mmcal::error
