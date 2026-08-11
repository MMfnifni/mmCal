#pragma once

#include <cstddef>
#include <iostream>
#include <string_view>
#include <utility>

namespace mmcal::tests {

class TestRunner final {
public:
    void expect(bool condition, std::string_view name) {
        if (condition) {
            ++passed_;
            std::cout << "[PASS] " << name << '\n';
            return;
        }

        ++failed_;
        std::cout << "[FAIL] " << name << '\n';
    }

    template <class Actual, class Expected>
    void expectEqual(
        const Actual& actual,
        const Expected& expected,
        std::string_view name) {
        const bool equal = actual == expected;
        expect(equal, name);

        if (!equal) {
            std::cout << "       actual:   " << actual << '\n';
            std::cout << "       expected: " << expected << '\n';
        }
    }

    template <class Exception, class Function>
    void expectThrows(Function&& function, std::string_view name) {
        bool expectedExceptionCaught = false;

        try {
            std::forward<Function>(function)();
        }
        catch (const Exception&) {
            expectedExceptionCaught = true;
        }
        catch (...) {
            // 別種の例外が送出された場合も、このテストは失敗とする。
        }

        expect(expectedExceptionCaught, name);
    }

    [[nodiscard]] int result() const {
        std::cout << "\npassed: " << passed_
                  << ", failed: " << failed_ << '\n';
        return failed_ == 0 ? 0 : 1;
    }

private:
    std::size_t passed_ = 0;
    std::size_t failed_ = 0;
};

} // namespace mmcal::tests
