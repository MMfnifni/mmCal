#include "display_messages.hpp"
#include "formatter.hpp"

namespace mm::cal {

 std::string makeAssignDisplayMessage(const std::string &name, const Value &newValue, const SystemConfig &cfg, bool wasRedefined, const Value *oldValue) {

  if (wasRedefined && oldValue != nullptr) { return "[redefined: " + name + " := " + formatResult(newValue, cfg) + " (was " + formatResult(*oldValue, cfg) + ")]"; }

  return "[defined: " + name + " := " + formatResult(newValue, cfg) + "]";
 }

 std::string makeFunctionDefDisplayMessage(const std::string &name, const std::vector<std::string> &params) {

  std::string sig = name + "(";
  for (size_t i = 0; i < params.size(); ++i) {
   if (i > 0) sig += ", ";
   sig += params[i];
  }
  sig += ")";

  return "[defined: " + sig + "]";
 }

} // namespace mm::cal