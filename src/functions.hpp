#pragma once

#include "core.hpp"
#include "math_util.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace mm::cal {

 /* ============================
   地獄の無限関数定義編
   ============================ */

 void initFunctions(SystemConfig &cfg);

 Value evaluateFunction(const std::string &name, const std::vector<Value> &args, FunctionContext &ctx);

} // namespace mm::cal
