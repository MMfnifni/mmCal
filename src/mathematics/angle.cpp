// 角度単位とセッション既定角度
#include "angle.hpp"

namespace mmcal::mathematics {

std::optional<AngleUnit> AngleSemantics::parseUnit(std::string_view name) noexcept {
    if (name == "Deg" || name == "deg")
        return AngleUnit::Degree;
    if (name == "Rad" || name == "rad")
        return AngleUnit::Radian;
    if (name == "Grad" || name == "grad")
        return AngleUnit::Gradian;
    return std::nullopt;
}

std::string_view AngleSemantics::canonicalName(AngleUnit unit) noexcept {
    switch (unit) {
    case AngleUnit::Degree:
        return "Deg";
    case AngleUnit::Radian:
        return "Rad";
    case AngleUnit::Gradian:
        return "Grad";
    }

    return "Deg";
}

const AngleSemantics& defaultAngleSemantics() {
    static const AngleSemantics semantics{AngleUnit::Radian};
    return semantics;
}

} // namespace mmcal::mathematics
