// 保証付き計算の精度不足エラー
#pragma once

#include <stdexcept>
#include <string>

namespace mmcal::approximation {

// 「数学的に未定義」ではなく、「現在の作業精度では必要な符号・非零性・branchを
// まだ証明できない」ことを表す内部例外。Nはこれを捕捉して作業precisionを増やす。
class PrecisionInsufficient final : public std::runtime_error {
public:
    explicit PrecisionInsufficient(std::string message)
        : std::runtime_error(std::move(message)) {}
};

} // namespace mmcal::approximation
