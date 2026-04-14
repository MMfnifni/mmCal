#pragma once

#include "core.hpp"
#include "evaluate.hpp"
#include "math_util.hpp"
#include <algorithm>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <functional>
#include <numeric>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

namespace mm::cal {

 static std::filesystem::path path_from_utf8_string(const std::string &utf8) {
#ifdef _WIN32
  return std::filesystem::path(std::u8string(reinterpret_cast<const char8_t *>(utf8.data()), reinterpret_cast<const char8_t *>(utf8.data() + utf8.size())));
#else
  return std::filesystem::path(utf8);
#endif
 }

 void registerBasicMath(SystemConfig &cfg);
 void registerExpLog(SystemConfig &cfg);
 void registerTrig(SystemConfig &cfg);
 void registerHyperbolic(SystemConfig &cfg);
 void registerStatistics(SystemConfig &cfg);
 void registerSignal(SystemConfig &cfg);
 void registerGeoVec(SystemConfig &cfg);
 void registerVector(SystemConfig &cfg);
 void registerComplex(SystemConfig &cfg);
 void registerRandom(SystemConfig &cfg);
 void registerAreaVol(SystemConfig &cfg);
 void registerEngineering(SystemConfig &cfg);
 void registerMoldInjection(SystemConfig &cfg);
 void registerFinance(SystemConfig &cfg);
 void registerCalculus(SystemConfig &cfg);
 void registerOthers(SystemConfig &cfg);

 void initFunctions(SystemConfig &cfg);
 Value evaluateFunction(const std::string &name, const std::vector<Value> &args, FunctionContext &ctx);

 // ---- help ----
 const std::unordered_map<std::string, std::string> &functionHelpMap();
 std::string getFunctionHelp(const std::string &name);
 std::string getFunctionHelpIndex();

} // namespace mm::cal
