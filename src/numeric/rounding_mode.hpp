#pragma once

namespace mmcal::numeric {

// BigFloatおよび区間演算で共通に使う丸め方向。
// RealIntervalでは TowardNegative / TowardPositive を使って、真値を必ず含む下端・上端を構成する。
enum class RoundingMode {
    NearestEven,
    TowardZero,
    TowardPositive,
    TowardNegative
};

} // namespace mmcal::numeric
