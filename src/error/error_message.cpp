#include "error_message.hpp"

#include <algorithm>
#include <sstream>
#include <utility>

namespace mmcal::error {
namespace {

[[nodiscard]] std::string sourceExcerpt(
    std::string_view sourceText,
    source::SourceSpan span) {
    if (sourceText.empty() || span.begin.offset > sourceText.size())
        return {};

    const std::size_t searchOffset = span.begin.offset == 0 ? 0 : span.begin.offset - 1;
    const std::size_t lineStart = sourceText.rfind('\n', searchOffset);
    const std::size_t start = lineStart == std::string_view::npos ? 0 : lineStart + 1;
    const std::size_t lineEnd = sourceText.find('\n', span.begin.offset);
    const std::size_t end = lineEnd == std::string_view::npos ? sourceText.size() : lineEnd;
    const std::string_view line = sourceText.substr(start, end - start);

    const std::size_t caretOffset = std::min(span.begin.offset - start, line.size());
    std::size_t caretLength = 1;
    if (span.end.offset > span.begin.offset)
        caretLength = std::min(span.end.offset - span.begin.offset, line.size() - caretOffset);

    std::string result;
    result.reserve(line.size() + caretOffset + caretLength + 4);
    result.append(line);
    result.push_back('\n');
    result.append(caretOffset, ' ');
    result.append(std::max<std::size_t>(caretLength, 1), '^');
    return result;
}

void appendSourceReference(
    std::ostringstream& output,
    std::string_view label,
    const source::SourceReference& reference) {
    if (!reference.document)
        return;

    output << label << " In [" << reference.document->inputNumber() << "], "
           << reference.span.begin.line << ':' << reference.span.begin.column << ":\n";

    const std::string excerpt = sourceExcerpt(reference.document->text(), reference.span);
    if (!excerpt.empty())
        output << excerpt;
}

[[nodiscard]] bool sameSource(
    const ErrorTraceFrame& frame,
    std::string_view label,
    const source::SourceReference& source) noexcept {
    return frame.label == label
        && frame.source.document == source.document
        && frame.source.span == source.span;
}

} // namespace

CalcError::CalcError(CalcErrorType type, std::string message)
    : std::runtime_error(std::move(message)), type_(type) {}

CalcError::CalcError(
    CalcErrorType type,
    std::string message,
    source::SourceSpan span)
    : std::runtime_error(std::move(message)), type_(type), span_(span) {}

CalcError::CalcError(
    CalcErrorType type,
    std::string message,
    source::SourceReference source,
    std::string sourceLabel)
    : std::runtime_error(std::move(message)),
      type_(type),
      span_(source.span),
      document_(std::move(source.document)),
      sourceLabel_(std::move(sourceLabel)) {}

CalcErrorType CalcError::type() const noexcept {
    return type_;
}

const std::optional<source::SourceSpan>& CalcError::span() const noexcept {
    return span_;
}

const std::shared_ptr<const source::SourceDocument>& CalcError::document() const noexcept {
    return document_;
}

const std::string& CalcError::sourceLabel() const noexcept {
    return sourceLabel_;
}

const std::vector<ErrorTraceFrame>& CalcError::trace() const noexcept {
    return trace_;
}

void CalcError::attachSourceIfMissing(
    source::SourceReference source,
    std::string_view label) {
    if (!span_)
        span_ = source.span;
    if (!document_)
        document_ = std::move(source.document);
    if (sourceLabel_ == "At" && label != "At")
        sourceLabel_ = label;
}

void CalcError::attachDocumentIfMissing(
    std::shared_ptr<const source::SourceDocument> document) {
    if (!document_)
        document_ = std::move(document);
}

void CalcError::setSourceLabel(std::string label) {
    sourceLabel_ = std::move(label);
}

void CalcError::addTrace(std::string label, source::SourceReference source) {
    if (!source.document)
        return;
    if (!trace_.empty() && sameSource(trace_.back(), label, source))
        return;

    trace_.push_back(ErrorTraceFrame{std::move(label), std::move(source)});
}

std::string_view calcErrorTypeName(CalcErrorType type) noexcept {
    switch (type) {
    case CalcErrorType::Syntax:
        return "SyntaxError";
    case CalcErrorType::Domain:
        return "DomainError";
    case CalcErrorType::Type:
        return "TypeError";
    case CalcErrorType::Overflow:
        return "OverflowError";
    case CalcErrorType::Name:
        return "NameError";
    case CalcErrorType::Evaluation:
        return "EvaluationError";
    case CalcErrorType::Internal:
        return "InternalError";
    }

    return "CalcError";
}

std::string errorMessage(const CalcError& error, std::string_view sourceText) {
    std::ostringstream output;

    if (error.document()) {
        output << calcErrorTypeName(error.type()) << ": " << error.what();
        if (error.span()) {
            output << '\n';
            appendSourceReference(
                output,
                error.sourceLabel(),
                source::SourceReference{error.document(), *error.span()});
        }
    }
    else {
        output << calcErrorTypeName(error.type());
        if (error.span())
            output << " at " << error.span()->begin.line << ':' << error.span()->begin.column;
        output << ": " << error.what();

        if (error.span()) {
            const std::string excerpt = sourceExcerpt(sourceText, *error.span());
            if (!excerpt.empty())
                output << '\n' << excerpt;
        }
    }

    for (const ErrorTraceFrame& frame : error.trace()) {
        output << '\n';
        appendSourceReference(output, frame.label, frame.source);
    }

    return output.str();
}

[[noreturn]] void throwCalcError(CalcErrorType type, std::string message) {
    throw CalcError{type, std::move(message)};
}

[[noreturn]] void throwCalcError(
    CalcErrorType type,
    std::string message,
    source::SourceSpan span) {
    throw CalcError{type, std::move(message), span};
}

[[noreturn]] void throwCalcError(
    CalcErrorType type,
    std::string message,
    source::SourceReference source,
    std::string sourceLabel) {
    throw CalcError{
        type,
        std::move(message),
        std::move(source),
        std::move(sourceLabel)};
}

} // namespace mmcal::error
