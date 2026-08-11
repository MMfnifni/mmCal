// 仮定集合
#include "assumption_set.hpp"

#include <algorithm>
#include <utility>

namespace mmcal::mathematics {

AssumptionSet::AssumptionSet(std::vector<Predicate> predicates) {
    for (Predicate& predicate : predicates)
        add(std::move(predicate));
}

void AssumptionSet::add(Predicate predicate) {
    if (!contains(predicate))
        predicates_.push_back(std::move(predicate));
}

bool AssumptionSet::contains(const Predicate& predicate) const {
    return std::find(predicates_.begin(), predicates_.end(), predicate) != predicates_.end();
}

bool AssumptionSet::empty() const noexcept {
    return predicates_.empty();
}

std::size_t AssumptionSet::size() const noexcept {
    return predicates_.size();
}

std::span<const Predicate> AssumptionSet::predicates() const noexcept {
    return predicates_;
}

} // namespace mmcal::mathematics
