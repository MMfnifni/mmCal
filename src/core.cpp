#include "core.hpp"

namespace mm::cal {

 Value roundValue(const Value &v, const SystemConfig &cfg) {
  if (isInvalid(v)) return v;

  auto roundDouble = [&](double d) { return std::round(d / cnst_precision_inv) * cnst_precision_inv; };

  if (std::holds_alternative<double>(v)) { return roundDouble(std::get<double>(v)); }

  if (std::holds_alternative<std::complex<double>>(v)) {
   auto c = std::get<std::complex<double>>(v);
   return std::complex<double>(roundDouble(c.real()), roundDouble(c.imag()));
  }

  return InvalidValue{};
 }

} // namespace mm::cal