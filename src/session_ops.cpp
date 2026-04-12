#include "session_ops.hpp"

namespace mm::cal {

 bool unsetVariable(EvaluationContext &ectx, const std::string &name) { return ectx.session.globals.erase(name) > 0; }

 bool undefFunction(EvaluationContext &ectx, const std::string &name) { return ectx.session.userFunctions.erase(name) > 0; }

 // void clearHistory(std::vector<InputEntry> &history) { history.clear(); }

 void resetSession(EvaluationContext &ectx, std::vector<InputEntry> &history) {
  ectx.session.globals.clear();
  ectx.session.userFunctions.clear();
  ectx.sideEffects.clear();
  history.clear();
 }

} // namespace mm::cal