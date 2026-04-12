#pragma once

#include "core.hpp"
#include <string>
#include <vector>

namespace mm::cal {

 bool unsetVariable(EvaluationContext &ectx, const std::string &name);
 bool undefFunction(EvaluationContext &ectx, const std::string &name);

 void clearHistory(std::vector<InputEntry> &history);

 // cfg は維持, base は維持, variables / userFunctions / sideEffects / history は消す
 void resetSession(EvaluationContext &ectx, std::vector<InputEntry> &history);

} // namespace mm::cal