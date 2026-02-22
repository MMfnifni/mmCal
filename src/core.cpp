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

 bool Value::hasNestedMulti(size_t pos) const {
  const auto &mv = asMultiRef(pos);
  for (const auto &v : mv.elems()) {
   if (v.isMulti()) return true;
  }
  return false;
 }

 // multi util
 inline const MultiValue &Value::asMultiRef(size_t pos) const { return *asMulti(pos); }

 std::size_t Value::multiSize(size_t pos) const { return asMultiRef(pos).size(); }

 bool Value::multiEmpty(size_t pos) const { return asMultiRef(pos).size() == 0; }

 const Value &Value::multiAt(std::size_t i, size_t pos) const {
  const auto &mv = asMultiRef(pos);
  if (i >= mv.size()) throw CalcError(CalcErrorType::OutOfRange, "multivalue index out of range", pos);
  return mv[i];
 }

 bool Value::isEmptyMulti() const noexcept {
  if (!isMulti()) return false;
  return std::get<MultiPtr>(data_)->size() == 0;
 }

 void Value::checkFiniteImpl(const MultiPtr &mv, size_t pos) {
  if (!mv) return;
  for (const auto &e : mv->elems_)
   Value::checkFinite(e, pos);
 }

 template <class F> Value MultiValue::map(F &&f) const {
  std::vector<Value> res;
  res.reserve(elems_.size());

  for (const auto &v : elems_)
   res.push_back(f(v));

  return Value(std::make_shared<MultiValue>(std::move(res)));
 }

 bool MultiValue::hasNested() const noexcept {
  for (const auto &v : elems_)
   if (v.isMulti()) return true;
  return false;
 }
 std::vector<std::vector<double>> MultiValue::toMatrix() const {
  std::vector<std::vector<double>> matrix;
  for (const auto &row : elems_) {
   if (!row.isMulti()) { throw CalcError(CalcErrorType::TypeError, "Invalid matrix row", 0); }
   std::vector<double> row_vals;
   const auto &multi_row = row.asMultiRef(0);
   for (size_t i = 0; i < multi_row.size(); ++i) {
    row_vals.push_back(multi_row[i].asScalar(0));
   }
   matrix.push_back(std::move(row_vals));
  }
  return matrix;
 }

} // namespace mm::cal