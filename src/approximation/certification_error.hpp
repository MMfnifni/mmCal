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

// 数学的な値は存在するが，現在のcertified backendが必要な値域/branchを未実装であることを表す。
// DomainErrorと区別し，Nは式を保持してunsupported warningへ落とす。
class CertifiedBackendUnsupported final : public std::domain_error {
public:
    explicit CertifiedBackendUnsupported(std::string message)
        : std::domain_error(std::move(message)) {}
};

} // namespace mmcal::approximation
