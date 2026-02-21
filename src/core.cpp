#include "core.hpp"

namespace mm::cal {

 Value roundValue(const Value &v, const SystemConfig &cfg) {
  auto roundDouble = [&](double d) { return std::round(d / cnst_precision_inv) * cnst_precision_inv; };

  return std::visit(Overloaded{// double
                               [&](double d) -> Value { return roundDouble(d); },

                               // complex
                               [&](std::complex<double> c) -> Value { return std::complex<double>(roundDouble(c.real()), roundDouble(c.imag())); },

                               // multi (vectorized)
                               [&](std::shared_ptr<MultiValue> m) -> Value {
                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(m->elems_.size());

                                for (const auto &e : m->elems_)
                                 result->elems_.push_back(roundValue(e, cfg));

                                return result;
                               },

                               // string / invalid
                               [&](auto &&x) -> Value {
                                if constexpr (std::is_same_v<std::decay_t<decltype(x)>, std::monostate>) {
                                 return Value();
                                } else if constexpr (std::is_same_v<std::decay_t<decltype(x)>, std::string>) {
                                 return Value(x);
                                } else {
                                 return Value(x); // double / complex / MultiPtr だけ
                                }
                               }},
                    v.storage());
 }

 void Value::checkFiniteImpl(const MultiPtr &mv, size_t pos) {
  if (!mv) return;
  for (const auto &e : mv->elems_)
   Value::checkFinite(e, pos);
 }

} // namespace mm::cal