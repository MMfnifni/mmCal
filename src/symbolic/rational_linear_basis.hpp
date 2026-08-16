#pragma once

#include "numeric/rational.hpp"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace mmcal::symbolic::detail {

// Q上のcolumn basisを逐次row-echelon化する小規模exact線形基盤。
// Krylov列のfirst dependency検出と，確立済みbasisへの座標変換を同じ消去状態で行う。
class RationalLinearBasis final {
public:
    explicit RationalLinearBasis(std::size_t dimension)
        : dimension_(dimension) {}

    // vectorを次のcolumnとして追加する。独立ならnullopt，従属なら
    // sum_i relation[i] * column_i == 0 を満たすrelationを返す。
    [[nodiscard]] std::optional<std::vector<numeric::Rational>> append(
        std::span<const numeric::Rational> vector) {
        if (vector.size() != dimension_)
            throw std::logic_error("Rational linear basis dimension mismatch");

        std::vector<numeric::Rational> work(vector.begin(), vector.end());
        std::vector<numeric::Rational> representation(appended_ + 1);
        representation.back() = numeric::Rational{numeric::BigInt{1}};

        reduce(work, representation, false);

        std::size_t pivot = 0;
        while (pivot < dimension_ && work[pivot].isZero())
            ++pivot;

        ++appended_;
        if (pivot == dimension_)
            return representation;

        const numeric::Rational pivotValue = work[pivot];
        for (std::size_t row = pivot; row < dimension_; ++row)
            work[row] /= pivotValue;
        for (numeric::Rational& coefficient : representation)
            coefficient /= pivotValue;

        basis_.push_back(BasisVector{
            pivot, std::move(work), std::move(representation)});
        return std::nullopt;
    }

    // 現basis列の一意な線形結合としてtargetを表す。spanされていなければnullopt。
    [[nodiscard]] std::optional<std::vector<numeric::Rational>> coordinates(
        std::span<const numeric::Rational> target) const {
        if (target.size() != dimension_)
            throw std::logic_error("Rational linear basis target dimension mismatch");
        if (basis_.size() != appended_)
            return std::nullopt;

        std::vector<numeric::Rational> work(target.begin(), target.end());
        std::vector<numeric::Rational> result(appended_);
        reduce(work, result, true);

        if (std::any_of(work.begin(), work.end(),
                [](const numeric::Rational& value) { return !value.isZero(); }))
            return std::nullopt;
        return result;
    }

    [[nodiscard]] std::size_t rank() const noexcept { return basis_.size(); }
    [[nodiscard]] std::size_t columnCount() const noexcept { return appended_; }

private:
    struct BasisVector final {
        std::size_t pivot = 0;
        std::vector<numeric::Rational> values;
        // このechelon vectorを追加済み元columnsで表した係数。
        std::vector<numeric::Rational> representation;
    };

    std::size_t dimension_ = 0;
    std::size_t appended_ = 0;
    std::vector<BasisVector> basis_;

    void reduce(
        std::vector<numeric::Rational>& work,
        std::vector<numeric::Rational>& representation,
        bool accumulateCoordinates) const {
        for (const BasisVector& basis : basis_) {
            const numeric::Rational factor = work[basis.pivot];
            if (factor.isZero())
                continue;

            for (std::size_t row = basis.pivot; row < dimension_; ++row)
                work[row] -= factor * basis.values[row];

            const std::size_t count = basis.representation.size();
            if (accumulateCoordinates) {
                for (std::size_t i = 0; i < count; ++i)
                    representation[i] += factor * basis.representation[i];
            }
            else {
                for (std::size_t i = 0; i < count; ++i)
                    representation[i] -= factor * basis.representation[i];
            }
        }
    }
};

} // namespace mmcal::symbolic::detail
