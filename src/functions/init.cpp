#include "ast.hpp"
#include "functions/functions.hpp"

namespace mm::cal {

 void initFunctions(SystemConfig &cfg) {
  cfg.functions.reserve(1024);
  registerBasicMath(cfg);
  registerExpLog(cfg);
  registerTrig(cfg);
  registerHyperbolic(cfg);
  registerGeoVec(cfg);
  registerVector(cfg);
  registerSignal(cfg);
  registerStatistics(cfg);
  registerComplex(cfg);
  registerRandom(cfg);
  registerAreaVol(cfg);
  registerEngineering(cfg);
  registerMoldInjection(cfg);
  registerFinance(cfg);
  registerCalculus(cfg);
  registerOthers(cfg);
 }

 Value evaluateFunction(const std::string &name, const std::vector<Value> &args, FunctionContext &ctx) {
  // 1. user-defined function
  if (auto uit = ctx.session.userFunctions.find(name); uit != ctx.session.userFunctions.end()) {

   const UserFunction &fn = uit->second;
   if (args.size() != fn.params.size()) { throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: wrong number of arguments for user-defined function: " + name, ctx.pos); }

   EvaluationContext local(ctx.session);
   for (size_t i = 0; i < fn.params.size(); ++i) {
    local.locals[fn.params[i]] = args[i];
   }
   return fn.body.get()->eval(local);
  }

  // 2. built-in
  auto it = ctx.session.cfg.functions.find(name);
  if (it == ctx.session.cfg.functions.end()) { throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), ctx.pos); }

  auto &f = it->second;
  if (!f.validArgc(static_cast<int>(args.size()))) { throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: wrong number of arguments for function: " + name, ctx.pos); }

  return f.f(const_cast<std::vector<Value> &>(args), ctx);
 }

} // namespace mm::cal
