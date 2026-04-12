#pragma once
#include "core.hpp"
#include <string>

namespace mm::cal {
 std::string dataExplainStructure(const Value &value);
 std::string dataExplain(Value value, const SystemConfig &cfg);
} // namespace mm::cal