#include "evaluate.hpp"
#include "display_messages.hpp"
#include "formatter.hpp"
#include "functions.hpp"
#include "lexer_parser.hpp"
#include "repl_command.hpp"
#include <iomanip>
#include <sstream>

namespace mm::cal {

 /* ============================
    評価マン(ToDo: Semanticの責務を兼ねている)
    ============================ */

 EvalResult evaluate(const std::string &src, EvaluationContext &ectx) {
  Parser p(ectx.session.cfg, src);
  auto ast = p.parse();

  if (p.cur.type != TokenType::End) { throw CalcError(CalcErrorType::SyntaxError, "Syntax Error: not closed", p.cur.pos); }

  Value v = ast->eval(ectx);
  return {v};
 }

 EvalResult evalLine(const std::string &line, SystemConfig &cfg, std::vector<InputEntry> &history, EvaluationContext &ectx) {
  if (!line.empty() && line.front() == ':') return evalCommandLine(line, cfg, history, ectx);

  EvalResult res{};

  ectx.sideEffects.clear();

  Parser p(cfg, line);
  auto ast = p.parse();

  if (p.cur.type != TokenType::End) { throw CalcError(CalcErrorType::SyntaxError, "Syntax Error: not closed", p.cur.pos); }
  bool assignWasRedefinition = false;
  Value oldAssignedValue;

  if (auto *n = dynamic_cast<AssignNode *>(ast.get())) {
   if (ectx.callStack.empty()) {
    auto it = ectx.session.globals.find(n->name);
    if (it != ectx.session.globals.end()) {
     assignWasRedefinition = true;
     oldAssignedValue = it->second;
    }
   } else {
    auto it = ectx.locals.find(n->name);
    if (it != ectx.locals.end()) {
     assignWasRedefinition = true;
     oldAssignedValue = it->second;
    }
   }
  }

  Value v = ast->eval(ectx);

  if (auto *n = dynamic_cast<AssignNode *>(ast.get())) {
   res.value = v;
   res.displayOverride = makeAssignDisplayMessage(n->name, v, cfg, assignWasRedefinition, assignWasRedefinition ? &oldAssignedValue : nullptr);
   return res;
  }

  if (auto *n = dynamic_cast<FunctionDefNode *>(ast.get())) {
   res.value = v;
   res.displayOverride = makeFunctionDefDisplayMessage(n->name, n->params);
   return res;
  }

  res.value = v;
  return res;
 }

} // namespace mm::cal