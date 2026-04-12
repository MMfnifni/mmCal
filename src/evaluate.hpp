#pragma once

#include "core.hpp"
#include <string>

namespace mm::cal {

 /* ============================
    評価マン
    ============================ */
 struct AssignNode;
 struct FunctionDefNode;

 EvalResult evaluate(const std::string &src, EvaluationContext &ectx);
 EvalResult evalLine(const std::string &line, SystemConfig &cfg, std::vector<InputEntry> &history, EvaluationContext &ectx);

} // namespace mm::cal
