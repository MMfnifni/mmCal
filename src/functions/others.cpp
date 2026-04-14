#include "constants.hpp"
#include "explain_formatter.hpp"
#include "formatter.hpp"
#include "functions/functions.hpp"
#include "unit_conv.hpp"

namespace mm::cal {

 void registerOthers(SystemConfig &cfg) {
  cfg.functions["cnst"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            // v[0] が string であることを要求
                            const std::string &name = v[0].asString(ctx.pos);
                            auto it = constants_dic.find(name);
                            if (it == constants_dic.end()) { throw CalcError(CalcErrorType::UnknownIdentifier, "unknown constant: " + name, ctx.pos); }
                            return it->second;
                           }};
  cfg.functions["fib"] = {1, 1, [](auto &v, auto &ctx) -> Value { return (double)fibULL(requireInt(v[0], ctx.pos), ctx.pos); }};
  cfg.functions["isprime"] = {1, 1, [](auto &v, auto &ctx) -> Value { return isPrimeLL(requireInt(v[0], ctx.pos)) ? 1.0 : 0.0; }};
  cfg.functions["nextprime"] = {1, 1, [](auto &v, auto &ctx) {
                                 long long n = requireInt(v[0], ctx.pos) + 1;
                                 while (!isPrimeLL(n))
                                  ++n;
                                 return (double)n;
                                }};
  cfg.functions["prevprime"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                                 int n = (int)std::floor(v[0].asScalar(ctx.pos));
                                 while (--n >= 2)
                                  if (isPrimeLL(n)) return (double)n;
                                 throwDomain(ctx.pos);
                                }};
  cfg.functions["if"] = {3, 3, [](auto &v, auto &ctx) -> Value { return (v[0].asScalar(ctx.pos) != 0.0 ? v[1] : v[2]); }};

  cfg.functions["convert"] = {2, 3, [](auto &v, auto &ctx) -> Value {
                               double value = v[0].asScalar(ctx.pos);
                               if (v.size() == 2) {
                                // 引数2つ: SI単位系に変換
                                std::string unit = v[1].asString(ctx.pos);
                                // SI単位系に変換
                                if (unit == "m") return Value(value);
                                else if (unit == "kg") return Value(value);
                                else if (unit == "s") return Value(value);
                                else if (unit == "K") return Value(value);
                                else if (unit == "A") return Value(value);
                                else if (unit == "mol") return Value(value);
                                else if (unit == "cd") return Value(value);
                                else {
                                 // 他の単位もSI単位系に変換
                                 for (const auto &[quantity, units] : unitTable) {
                                  auto it = units.find(unit);
                                  if (it != units.end()) {
                                   // SI単位系に変換
                                   double si_value = value * it->second;
                                   // SI単位系の値を返す
                                   return Value(si_value);
                                  }
                                 }
                                 throw CalcError(CalcErrorType::UnknownIdentifier, "unknown unit: " + unit, ctx.pos);
                                }
                               } else {
                                // 引数3つ: 通常の変換
                                std::string from_unit = v[1].asString(ctx.pos);
                                std::string to_unit = v[2].asString(ctx.pos);

                                // 温度特殊処理(羊のせいでえらい目にあった)
                                if ((from_unit == "C" || from_unit == "F" || from_unit == "K") && (to_unit == "C" || to_unit == "F" || to_unit == "K")) {

                                 // ケルビンに変換
                                 double kelvin;
                                 if (from_unit == "C") {
                                  kelvin = value + 273.15;
                                 } else if (from_unit == "F") {
                                  kelvin = (value - 32) * 5 / 9 + 273.15;
                                 } else {
                                  kelvin = value;
                                 }

                                 // ターゲットに変換
                                 if (to_unit == "C") {
                                  return Value(kelvin - 273.15);
                                 } else if (to_unit == "F") {
                                  return Value((kelvin - 273.15) * 9 / 5 + 32);
                                 } else {
                                  return Value(kelvin);
                                 }
                                }

                                // 他の単位の変換
                                for (const auto &[quantity, units] : unitTable) {
                                 auto from_it = units.find(from_unit);
                                 auto to_it = units.find(to_unit);

                                 if (from_it != units.end() && to_it != units.end()) {
                                  // SI単位系に変換してから変換
                                  double si_value = value * from_it->second;
                                  double result = si_value / to_it->second;
                                  return Value(result);
                                 }
                                }

                                throw CalcError(CalcErrorType::UnknownIdentifier, "unknown unit combination: " + from_unit + " to " + to_unit, ctx.pos);
                               }
                              }};

  // ---- history / state ----
  cfg.functions["In"] = {1, 1, [](auto &v, FunctionContext &ctx) -> Value {
                          int i = (int)requireInt(v[0], ctx.pos);
                          if (i <= 0 || i > (int)ctx.session.history.size()) { throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos); }

                          EvaluationContext ectx{ctx.session};
                          return evaluate(ctx.session.history[i - 1].expr, ectx).value;
                         }};
  cfg.functions["Out"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                           int i = (int)requireInt(v[0], ctx.pos);
                           if (i <= 0 || i > (int)ctx.session.history.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos);
                           return ctx.session.history[i - 1].value;
                          }};
  cfg.functions["Prev"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            int k = (int)requireInt(v[0], ctx.pos);
                            int idx = ctx.session.base - k + 1;
                            if (idx <= 0 || idx > (int)ctx.session.history.size()) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), ctx.pos);
                            return ctx.session.history[idx - 1].value;
                           }};
  cfg.functions["Clear"] = {0, 0, [](auto &, auto &) -> Value { throw ControlRequest{ControlRequest::Kind::Clear}; }};
  cfg.functions["Exit"] = {0, 0, [](auto &, auto &) -> Value { throw ControlRequest{ControlRequest::Kind::Exit}; }};

  cfg.functions["filein"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                              namespace fs = std::filesystem;

                              const std::string &utf8 = v[0].asString(ctx.pos);
                              fs::path path = fs::path{std::u8string_view(reinterpret_cast<const char8_t *>(utf8.data()), utf8.size())};

                              std::ifstream ifs;
                              ifs.exceptions(std::ios::failbit | std::ios::badbit);

                              try {
                               ifs.open(path, std::ios::in);
                              } catch (...) { throw CalcError(CalcErrorType::IOError, "IOError: cannot open file", ctx.pos); }

                              std::string content;

                              try {
                               ifs.seekg(0, std::ios::end);
                               std::size_t size = static_cast<std::size_t>(ifs.tellg());
                               ifs.seekg(0, std::ios::beg);

                               content.resize(size);
                               ifs.read(content.data(), size);
                              } catch (...) { throw CalcError(CalcErrorType::IOError, "IOError: read failure", ctx.pos); }

                              try {
                               std::vector<InputEntry> localHistory = ctx.session.history;
                               int localBase = ctx.session.base;
                               std::unordered_map<std::string, Value> localGlobals;
                               std::unordered_map<std::string, UserFunction> localUserFunctions;

                               SessionState session{ctx.session.cfg, localHistory, localBase, localGlobals, localUserFunctions};
                               EvaluationContext ectx{session};

                               EvalResult result = evaluate(content, ectx);
                               return result.value;
                              } catch (...) { throw CalcError(CalcErrorType::IOError, "IOError: file content cannot be evaluated", ctx.pos); }
                             }};
  cfg.functions["fileout"] = {2, 2, [](auto &v, auto &ctx) -> Value {
                               std::ostringstream oss;
                               serializeValue(v[0], oss, ctx.session.cfg);
                               std::string content = std::move(oss).str();

                               ctx.sideEffects.push_back({
                                   SideEffect::Kind::FileWrite,
                                   v[1].asString(ctx.pos), // path
                                   std::move(content)      // serialized content
                               });

                               return v[0];
                              }};

  cfg.functions["clip"] = {1, 1, [](auto &v, auto &ctx) -> Value {
                            std::string s = formatResult(v[0], ctx.session.cfg);

                            if (s.size() >= 10 * 1024 * 1024) calcWarn("clip: data size exceeds 10MB");
                            if (s.size() >= 100 * 1024 * 1024) { throw CalcError(CalcErrorType::DomainError, "DomainError: clip: data size exceeds 100MB safety limit", ctx.pos); }

                            ctx.sideEffects.push_back({SideEffect::Kind::ClipboardCopy, std::move(s), {}});

                            return v[0];
                           }};

  cfg.functions["silent"] = {1, 1, [](auto &v, FunctionContext &ctx) -> Value {
                              ctx.sideEffects.push_back({SideEffect::Kind::SuppressDisplay, "", ""});
                              return v[0];
                             }};
  cfg.functions["explain"] = {1, 1, [](auto &v, FunctionContext &ctx) -> Value {
                               ctx.sideEffects.push_back({SideEffect::Kind::Explain, dataExplain(v[0], ctx.session.cfg), {}});
                               return v[0];
                              }};
 }

} // namespace mm::cal
