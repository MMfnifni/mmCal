// 小規模整数lattice用のexact row-LLL。
#include "exact_lll.hpp"

#include "numeric/rational.hpp"
#include "numeric/rational_rounding.hpp"

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace mmcal::linear_algebra {
namespace {

using numeric::BigInt;
using numeric::Rational;

struct GramSchmidtData final {
    std::vector<std::vector<Rational>> mu;
    std::vector<Rational> normSquared;
};

enum class GramSchmidtStatus {
    Ready,
    RankDeficient,
    BitLimitExceeded
};

[[nodiscard]] bool observe(
    const BigInt& value,
    const ExactLllOptions& options,
    std::size_t& maximumObservedBits) noexcept {
    const std::size_t bits = value.bitLength();
    maximumObservedBits = std::max(maximumObservedBits, bits);
    return bits <= options.maximumIntermediateBits;
}

[[nodiscard]] bool observe(
    const Rational& value,
    const ExactLllOptions& options,
    std::size_t& maximumObservedBits) noexcept {
    return observe(value.numerator(), options, maximumObservedBits)
        && observe(value.denominator(), options, maximumObservedBits);
}

[[nodiscard]] bool observe(
    const IntegerLatticeBasis& basis,
    const ExactLllOptions& options,
    std::size_t& maximumObservedBits) noexcept {
    for (const IntegerLatticeRow& row : basis)
        for (const BigInt& value : row)
            if (!observe(value, options, maximumObservedBits))
                return false;
    return true;
}

[[nodiscard]] bool integerDot(
    const IntegerLatticeRow& lhs,
    const IntegerLatticeRow& rhs,
    const ExactLllOptions& options,
    std::size_t& maximumObservedBits,
    BigInt& result) {
    result = BigInt{};
    for (std::size_t column = 0; column < lhs.size(); ++column) {
        BigInt term = lhs[column] * rhs[column];
        if (!observe(term, options, maximumObservedBits))
            return false;
        result += term;
        if (!observe(result, options, maximumObservedBits))
            return false;
    }
    return true;
}

[[nodiscard]] GramSchmidtStatus exactGramSchmidt(
    const IntegerLatticeBasis& basis,
    const ExactLllOptions& options,
    std::size_t& maximumObservedBits,
    GramSchmidtData& result) {
    const std::size_t rank = basis.size();
    result.mu.assign(rank, std::vector<Rational>(rank));
    result.normSquared.assign(rank, Rational{});

    for (std::size_t row = 0; row < rank; ++row) {
        BigInt normInteger;
        if (!integerDot(
                basis[row], basis[row], options,
                maximumObservedBits, normInteger))
            return GramSchmidtStatus::BitLimitExceeded;
        Rational norm{std::move(normInteger)};

        for (std::size_t previous = 0; previous < row; ++previous) {
            BigInt dotInteger;
            if (!integerDot(
                    basis[row], basis[previous], options,
                    maximumObservedBits, dotInteger))
                return GramSchmidtStatus::BitLimitExceeded;
            Rational numerator{std::move(dotInteger)};
            for (std::size_t projection = 0;
                 projection < previous; ++projection) {
                Rational term = result.mu[row][projection]
                    * result.mu[previous][projection]
                    * result.normSquared[projection];
                if (!observe(term, options, maximumObservedBits))
                    return GramSchmidtStatus::BitLimitExceeded;
                numerator -= term;
                if (!observe(numerator, options, maximumObservedBits))
                    return GramSchmidtStatus::BitLimitExceeded;
            }
            if (result.normSquared[previous].isZero())
                return GramSchmidtStatus::RankDeficient;
            result.mu[row][previous] =
                numerator / result.normSquared[previous];
            if (!observe(
                    result.mu[row][previous], options,
                    maximumObservedBits))
                return GramSchmidtStatus::BitLimitExceeded;

            Rational projectionNorm = result.mu[row][previous]
                * result.mu[row][previous]
                * result.normSquared[previous];
            if (!observe(projectionNorm, options, maximumObservedBits))
                return GramSchmidtStatus::BitLimitExceeded;
            norm -= projectionNorm;
            if (!observe(norm, options, maximumObservedBits))
                return GramSchmidtStatus::BitLimitExceeded;
        }

        if (norm.isZero() || norm.numerator().isNegative())
            return GramSchmidtStatus::RankDeficient;
        result.normSquared[row] = std::move(norm);
    }
    return GramSchmidtStatus::Ready;
}

[[nodiscard]] IntegerLatticeBasis identityBasis(std::size_t size) {
    IntegerLatticeBasis result(
        size, IntegerLatticeRow(size, BigInt{}));
    for (std::size_t row = 0; row < size; ++row)
        result[row][row] = BigInt{1};
    return result;
}

[[nodiscard]] bool subtractRowMultiple(
    const IntegerLatticeRow& target,
    const IntegerLatticeRow& source,
    const BigInt& multiple,
    const ExactLllOptions& options,
    std::size_t& maximumObservedBits,
    IntegerLatticeRow& result) {
    result.resize(target.size());
    for (std::size_t column = 0; column < target.size(); ++column) {
        BigInt product = multiple * source[column];
        if (!observe(product, options, maximumObservedBits))
            return false;
        result[column] = target[column] - product;
        if (!observe(result[column], options, maximumObservedBits))
            return false;
    }
    return true;
}

[[nodiscard]] ExactLllStatus gramFailureStatus(
    GramSchmidtStatus status) noexcept {
    return status == GramSchmidtStatus::RankDeficient
        ? ExactLllStatus::RankDeficient
        : ExactLllStatus::IntermediateBitLimitExceeded;
}

} // namespace

ExactLllResult exactLllReduceRows(
    IntegerLatticeBasis basis,
    ExactLllOptions options) {
    ExactLllResult result;
    result.basis = std::move(basis);
    if (result.basis.empty()) {
        result.status = ExactLllStatus::Reduced;
        return result;
    }

    const std::size_t rank = result.basis.size();
    const std::size_t columns = result.basis.front().size();
    if (columns == 0
        || std::any_of(
            result.basis.begin(), result.basis.end(),
            [&](const IntegerLatticeRow& row) {
                return row.size() != columns;
            })) {
        result.status = ExactLllStatus::DimensionMismatch;
        return result;
    }
    // well-formedな入力では，dimension budgetで開始前に止まる場合も
    // basis == transformation * inputという公開contractを保つ。
    result.transformation = identityBasis(rank);
    if (rank > options.maximumRank) {
        result.status = ExactLllStatus::RankLimitExceeded;
        return result;
    }
    if (columns > options.maximumColumns) {
        result.status = ExactLllStatus::ColumnLimitExceeded;
        return result;
    }
    if (rank > columns) {
        result.status = ExactLllStatus::RankDeficient;
        return result;
    }

    if (!observe(result.basis, options, result.maximumObservedBits)
        || !observe(
            result.transformation, options,
            result.maximumObservedBits)) {
        result.status = ExactLllStatus::IntermediateBitLimitExceeded;
        return result;
    }

    GramSchmidtData gram;
    GramSchmidtStatus gramStatus = exactGramSchmidt(
        result.basis, options, result.maximumObservedBits, gram);
    if (gramStatus != GramSchmidtStatus::Ready) {
        result.status = gramFailureStatus(gramStatus);
        return result;
    }

    const Rational delta{BigInt{3}, BigInt{4}};
    std::size_t row = 1;
    while (row < rank) {
        for (std::size_t previous = row; previous-- > 0;) {
            const BigInt multiple = numeric::roundToNearestEvenInteger(
                gram.mu[row][previous]);
            if (multiple.isZero())
                continue;
            if (result.sizeReductions >= options.maximumSizeReductions) {
                result.status = ExactLllStatus::SizeReductionLimitExceeded;
                return result;
            }

            IntegerLatticeRow nextBasisRow;
            IntegerLatticeRow nextTransformRow;
            if (!subtractRowMultiple(
                    result.basis[row], result.basis[previous], multiple,
                    options, result.maximumObservedBits, nextBasisRow)
                || !subtractRowMultiple(
                    result.transformation[row],
                    result.transformation[previous], multiple,
                    options, result.maximumObservedBits,
                    nextTransformRow)) {
                result.status = ExactLllStatus::IntermediateBitLimitExceeded;
                return result;
            }
            result.basis[row] = std::move(nextBasisRow);
            result.transformation[row] = std::move(nextTransformRow);
            ++result.sizeReductions;

            const Rational exactMultiple{multiple};
            for (std::size_t projection = 0;
                 projection < previous; ++projection) {
                Rational next = gram.mu[row][projection]
                    - exactMultiple * gram.mu[previous][projection];
                if (!observe(next, options, result.maximumObservedBits)) {
                    result.status =
                        ExactLllStatus::IntermediateBitLimitExceeded;
                    return result;
                }
                gram.mu[row][projection] = std::move(next);
            }
            gram.mu[row][previous] -= exactMultiple;
            if (!observe(
                    gram.mu[row][previous], options,
                    result.maximumObservedBits)) {
                result.status = ExactLllStatus::IntermediateBitLimitExceeded;
                return result;
            }
        }

        const Rational muSquared = gram.mu[row][row - 1]
            * gram.mu[row][row - 1];
        if (!observe(muSquared, options, result.maximumObservedBits)) {
            result.status = ExactLllStatus::IntermediateBitLimitExceeded;
            return result;
        }
        const Rational right =
            (delta - muSquared) * gram.normSquared[row - 1];
        if (!observe(right, options, result.maximumObservedBits)) {
            result.status = ExactLllStatus::IntermediateBitLimitExceeded;
            return result;
        }
        if (gram.normSquared[row] >= right) {
            ++row;
            continue;
        }

        if (result.swaps >= options.maximumSwaps) {
            result.status = ExactLllStatus::SwapLimitExceeded;
            return result;
        }

        const Rational adjacentMu = gram.mu[row][row - 1];
        const Rational previousNorm = gram.normSquared[row - 1];
        const Rational currentNorm = gram.normSquared[row];
        Rational nextPreviousNorm = currentNorm
            + adjacentMu * adjacentMu * previousNorm;
        Rational nextAdjacentMu =
            adjacentMu * previousNorm / nextPreviousNorm;
        Rational nextCurrentNorm =
            previousNorm * currentNorm / nextPreviousNorm;
        if (!observe(nextPreviousNorm, options, result.maximumObservedBits)
            || !observe(nextAdjacentMu, options, result.maximumObservedBits)
            || !observe(nextCurrentNorm, options, result.maximumObservedBits)) {
            result.status = ExactLllStatus::IntermediateBitLimitExceeded;
            return result;
        }

        std::vector<std::pair<Rational, Rational>> followingUpdates;
        followingUpdates.reserve(rank - row - 1);
        for (std::size_t following = row + 1;
             following < rank; ++following) {
            Rational nextCurrentMu = gram.mu[following][row - 1]
                - adjacentMu * gram.mu[following][row];
            Rational nextPreviousMu = gram.mu[following][row]
                + nextAdjacentMu * nextCurrentMu;
            if (!observe(nextCurrentMu, options, result.maximumObservedBits)
                || !observe(
                    nextPreviousMu, options,
                    result.maximumObservedBits)) {
                result.status = ExactLllStatus::IntermediateBitLimitExceeded;
                return result;
            }
            followingUpdates.emplace_back(
                std::move(nextPreviousMu), std::move(nextCurrentMu));
        }

        std::swap(result.basis[row], result.basis[row - 1]);
        std::swap(
            result.transformation[row],
            result.transformation[row - 1]);
        ++result.swaps;
        for (std::size_t projection = 0;
             projection + 1 < row; ++projection)
            std::swap(
                gram.mu[row][projection],
                gram.mu[row - 1][projection]);
        gram.normSquared[row - 1] = std::move(nextPreviousNorm);
        gram.normSquared[row] = std::move(nextCurrentNorm);
        gram.mu[row][row - 1] = std::move(nextAdjacentMu);
        for (std::size_t following = row + 1;
             following < rank; ++following) {
            auto& [nextPreviousMu, nextCurrentMu] =
                followingUpdates[following - row - 1];
            gram.mu[following][row - 1] = std::move(nextPreviousMu);
            gram.mu[following][row] = std::move(nextCurrentMu);
        }
        if (row > 1)
            --row;
    }

    result.status = ExactLllStatus::Reduced;
    return result;
}

} // namespace mmcal::linear_algebra
