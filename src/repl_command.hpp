#pragma once

#include "core.hpp"
#include <string>
#include <vector>

namespace mm::cal {

 EvalResult evalCommandLine(const std::string &line, SystemConfig &cfg, std::vector<InputEntry> &history, EvaluationContext &ectx);

} // namespace mm::cal