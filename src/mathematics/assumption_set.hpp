#pragma once

#include "predicate.hpp"

#include <span>
#include <vector>

namespace mmcal::mathematics {

// Solverの条件分岐と、将来のAssuming/Refineが共有する前提集合。
// ここ自体は推論器ではなく、前提を失わず保持する小さな値型に徹する。
class AssumptionSet final {
public:
    AssumptionSet() = default;
    explicit AssumptionSet(std::vector<Predicate> predicates);

    void add(Predicate predicate);
    [[nodiscard]] bool contains(const Predicate& predicate) const;
    [[nodiscard]] bool empty() const noexcept;
    [[nodiscard]] std::size_t size() const noexcept;
    [[nodiscard]] std::span<const Predicate> predicates() const noexcept;

    [[nodiscard]] bool operator==(const AssumptionSet&) const = default;

private:
    std::vector<Predicate> predicates_;
};

} // namespace mmcal::mathematics
