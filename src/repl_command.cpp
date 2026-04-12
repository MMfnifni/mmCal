#include "repl_command.hpp"
#include "evaluate.hpp"
#include "formatter.hpp"
#include "session_ops.hpp"
#include <algorithm>
#include <sstream>
#include <string>
#include <vector>

namespace mm::cal {

 static std::vector<std::string> splitWords(const std::string &line) {
  std::vector<std::string> words;
  std::istringstream iss(line);
  std::string w;
  while (iss >> w)
   words.push_back(w);
  return words;
 }

 static std::string joinNamesSorted(const std::vector<std::string> &items, const std::string &indent = "  ") {
  if (items.empty()) return indent + "(none)\n";

  std::string out;
  for (const auto &s : items) {
   out += indent;
   out += s;
   out += "\n";
  }
  return out;
 }

 static std::string formatFunctionSignature(const std::string &name, const UserFunction &fn) {
  std::string sig = name;
  sig += "(";
  for (size_t i = 0; i < fn.params.size(); ++i) {
   if (i > 0) sig += ", ";
   sig += fn.params[i];
  }
  sig += ")";
  return sig;
 }

 EvalResult evalCommandLine(const std::string &line, SystemConfig &cfg, std::vector<InputEntry> &history, EvaluationContext &ectx) {
  (void)cfg;
  EvalResult res{};

  auto words = splitWords(line);
  if (words.empty()) { throw CalcError(CalcErrorType::SyntaxError, "empty command", 0); }

  const std::string &cmd = words[0];

  // -----------------------
  // :defs
  // -----------------------
  if (cmd == ":defs") {
   if (words.size() != 1) { throw CalcError(CalcErrorType::SyntaxError, "Usage: :defs", 0); }

   std::vector<std::string> varLines;
   varLines.reserve(ectx.session.globals.size());

   for (const auto &kv : ectx.session.globals) {
    const std::string &name = kv.first;
    const Value &val = kv.second;

    varLines.push_back(name + " := " + formatResultMulti(val, cfg));
   }

   std::sort(varLines.begin(), varLines.end());

   std::vector<std::string> funcSigs;
   funcSigs.reserve(ectx.session.userFunctions.size());
   for (const auto &[name, fn] : ectx.session.userFunctions) {
    funcSigs.push_back(formatFunctionSignature(name, fn));
   }
   std::sort(funcSigs.begin(), funcSigs.end());

   std::string out;
   out += "definitions =>\n";
   out += "Variables:\n";
   out += joinNamesSorted(varLines);
   out += "Functions:\n";
   out += joinNamesSorted(funcSigs);

   res.displayOverride = out;
   return res;
  }

  // -----------------------
  // :unset <var>
  // -----------------------
  if (cmd == ":unset") {
   if (words.size() != 2) { throw CalcError(CalcErrorType::SyntaxError, "Usage: :unset <variable]", 0); }

   const std::string &name = words[1];

   if (!unsetVariable(ectx, name)) { throw CalcError(CalcErrorType::UnknownIdentifier, "unknown variable: " + name, 0); }

   res.displayOverride = "[unset: " + name + "]";
   return res;
  }

  // -----------------------
  // :undef <func>
  // -----------------------
  if (cmd == ":undef") {
   if (words.size() != 2) { throw CalcError(CalcErrorType::SyntaxError, "Usage: :undef <function]", 0); }

   const std::string &name = words[1];

   if (!undefFunction(ectx, name)) { throw CalcError(CalcErrorType::UnknownIdentifier, "unknown function: " + name, 0); }

   res.displayOverride = "[undefined: " + name + "]";
   return res;
  }

  // -----------------------
  // :exit
  // -----------------------
  if (cmd == ":exit") {
   if (words.size() != 1) { throw CalcError(CalcErrorType::SyntaxError, "Usage: :exit", 0); }
   throw ControlRequest(ControlRequest::Kind::Exit);
  }

  // -----------------------
  // :clear
  // -----------------------
  if (cmd == ":clear") {
   if (words.size() != 1) { throw CalcError(CalcErrorType::SyntaxError, "Usage: :clear", 0); }
   throw ControlRequest(ControlRequest::Kind::Clear);
  }

  throw CalcError(CalcErrorType::InvalidArgument, "unknown command: " + cmd, 0);
 }

} // namespace mm::cal