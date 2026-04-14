#include "functions/functions.hpp"

namespace mm::cal {

 void initFunctions(SystemConfig &cfg) {
  cfg.functions.reserve(400);
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
  registerOthers(cfg);
 }

 Value evaluateFunction(const std::string &name, const std::vector<Value> &args, FunctionContext &ctx) {
  auto it = ctx.session.cfg.functions.find(name);
  if (it == ctx.session.cfg.functions.end()) throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), ctx.pos);

  const FunctionDef &f = it->second;

  if (!f.validArgc(static_cast<int>(args.size()))) throw CalcError(CalcErrorType::SyntaxError, errorMessage(CalcErrorType::SyntaxError), ctx.pos);

  Value result = f.f(args, ctx);
  result.validateFinite(ctx.pos);
  return result;
 }

} // namespace mm::cal
