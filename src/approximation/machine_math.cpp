// 機械精度補助演算
#include "machine_math.hpp"

#include <cmath>
#include <stdexcept>

namespace mmcal::approximation {

double MachineMath::evaluate(ElementaryFunction function, double value) {
    switch (function) {
    case ElementaryFunction::Sin:
        return std::sin(value);
    case ElementaryFunction::Cos:
        return std::cos(value);
    case ElementaryFunction::Tan:
        return std::tan(value);
    case ElementaryFunction::Asin:
        return std::asin(value);
    case ElementaryFunction::Acos:
        return std::acos(value);
    case ElementaryFunction::Atan:
        return std::atan(value);
    case ElementaryFunction::Exp:
        return std::exp(value);
    case ElementaryFunction::Log:
        return std::log(value);
    case ElementaryFunction::Log10:
        return std::log10(value);
    case ElementaryFunction::Sqrt:
        return std::sqrt(value);
    }

    throw std::logic_error("Unknown elementary function");
}

double MachineMath::atan2(double y, double x) {
    return std::atan2(y, x);
}

double MachineMath::pow(double base, double exponent) {
    return std::pow(base, exponent);
}

} // namespace mmcal::approximation
