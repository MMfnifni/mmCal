// CalcErrorの表示文字列生成の回帰テスト
#include "error_message_tests.hpp"

#include "error/error_message.hpp"
#include "test_framework.hpp"

#include <memory>
#include <string>

namespace mmcal::tests {

void runErrorMessageTests(TestRunner& tests) {
    using error::CalcError;
    using error::CalcErrorType;
    using source::SourcePosition;
    using source::SourceSpan;

    const CalcError syntaxError{
        CalcErrorType::Syntax,
        "Expected ')'",
        SourceSpan{SourcePosition{4, 1, 5}, SourcePosition{5, 1, 6}}};
    const std::string formatted = error::errorMessage(syntaxError, "1 + (2");

    tests.expect(
        formatted.find("SyntaxError at 1:5: Expected ')'") != std::string::npos,
        "error message includes type and position");
    tests.expect(
        formatted.find("1 + (2\n    ^") != std::string::npos,
        "error message includes source excerpt");
    const auto sourceText = std::make_shared<const std::string>("1 / 0");
    const auto document = std::make_shared<const source::SourceDocument>(7, sourceText);
    CalcError trackedError{
        CalcErrorType::Domain,
        "Division by zero",
        source::SourceReference{
            document,
            SourceSpan{SourcePosition{4, 1, 5}, SourcePosition{5, 1, 6}}}};
    trackedError.addTrace(
        "Called from",
        source::SourceReference{
            document,
            SourceSpan{SourcePosition{0, 1, 1}, SourcePosition{5, 1, 6}}});
    const std::string tracked = error::errorMessage(trackedError);
    tests.expect(tracked.find("At In [7], 1:5:") != std::string::npos,
        "tracked error message includes input number");
    tests.expect(tracked.find("Called from In [7], 1:1:") != std::string::npos,
        "tracked error message includes trace frame");

    tests.expect(
        error::calcErrorTypeName(CalcErrorType::Domain) == "DomainError",
        "domain error type name");
}

} // namespace mmcal::tests
