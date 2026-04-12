#pragma once

#include "core.hpp"
#include <complex>
#include <string>

namespace mm::cal {

 inline constexpr int INDENT_WIDTH = 2;

 inline double normalizeZero(double x) { return (x == 0.0 ? 0.0 : x); } // -0 → 0

 inline bool nearlyZero(double x, int precision) { return x < cnst_precision_inv && x > -cnst_precision_inv; }

 std::string formatReal(double x, const SystemConfig &cfg);
 std::string formatComplex(const std::complex<double> &c, const SystemConfig &cfg);

 void appendValue(const Value &v, const SystemConfig &cfg, std::string &out);

 std::string formatResult(const Value &v, const SystemConfig &cfg);
 std::string formatResultMulti(const Value &v, const SystemConfig &cfg, int indent = 0);

 inline void setupStream(std::ostringstream &oss, const SystemConfig &cfg) {
  if (cfg.precision >= 0) {
   oss << std::fixed << std::setprecision(cfg.precision);
  } else {
   oss << std::setprecision(15);
  }
 }

} // namespace mm::cal