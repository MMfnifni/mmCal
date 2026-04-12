#pragma once

#include "core.hpp"
#include <string>
#include <vector>

namespace mm::cal {

 std::string makeAssignDisplayMessage(const std::string &name, const Value &newValue, const SystemConfig &cfg, bool wasRedefined, const Value *oldValue);

 std::string makeFunctionDefDisplayMessage(const std::string &name, const std::vector<std::string> &params);

} // namespace mm::cal