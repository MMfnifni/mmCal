#include "evaluate.hpp"

namespace mm::cal {

 /* ============================
    評価マン
    ============================ */

 Value evaluate(const std::string &src, SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) {
  Parser p(cfg, src);
  std::unique_ptr<Parser::ASTNode> ast = p.parseCompare();
  if (p.cur.type != TokenType::End) throw CalcError(CalcErrorType::Syntax, errorMessage(CalcErrorType::Syntax), p.cur.pos);
  return ast->eval(cfg, hist, base);
 }

 CalcResult calcEval(const std::string &expr, SystemConfig &cfg, std::vector<InputEntry> &history, int base) {
  try {
   Value v = evaluate(expr, cfg, history, base);

   // 表示用に丸め
   Value vr = roundValue(v, cfg);
   std::string out = formatResult(vr, cfg);

   history.push_back({expr, vr});

   return {true, out, 0};

  } catch (const CalcError &e) { return {false, e.what(), e.pos}; } catch (const std::exception &e) {
   return {false, e.what(), 0};

  } catch (...) { return {false, "unknown error", 0}; }
 }

 EvalResult evalLine(const std::string &line, SystemConfig &cfg, RuntimeState &rt, std::vector<InputEntry> &history) {
  EvalResult res{};

  try {
   Value v = evaluate(line, cfg, history, static_cast<int>(history.size()));
   res.kind = EvalKind::Value;
   res.value = v;
   return res;
  } catch (const ClearRequest &) {
   rt.shouldClear = true;
   res.kind = EvalKind::Clear;
   return res;
  } catch (const ExitRequest &) {
   rt.shouldExit = true;
   std::cout << "bye...nara\n";
   res.kind = EvalKind::Exit;
   return res;
  }
 }

 /* ============================
   フォーマッタ
   ============================ */

 std::string formatReal(double x, const SystemConfig &cfg) {
  std::ostringstream oss;

  if (cfg.trimTrailingZeros) {
   oss << std::fixed << std::setprecision(cfg.precision);
  } else {
   oss << std::setprecision(cfg.precision);
  }

  oss << x;
  std::string s = oss.str();

  // fixed + trimTrailingZeros のとき末尾処理
  if (cfg.trimTrailingZeros && s.find('.') != std::string::npos) {
   while (!s.empty() && s.back() == '0')
    s.pop_back();
   if (!s.empty() && s.back() == '.') s.pop_back();
  }

  return s;
 }

 std::string formatComplex(const std::complex<double> &c, const SystemConfig &cfg) {
  double re = c.real();
  double im = c.imag();

  double eps = std::pow(10.0, -(cnst_precision));
  bool re0 = std::abs(re) < eps;
  bool im0 = std::abs(im) < eps;

  // 純実数
  if (im0) return formatReal(re, cfg);

  std::ostringstream oss;

  // 実部
  if (!re0) oss << formatReal(re, cfg);

  // 虚部の符号
  if (!re0 && im > 0) oss << "+";

  // 虚部の大きさ
  if (std::abs(im - 1.0) < eps) oss << "I";
  else if (std::abs(im + 1.0) < eps) oss << "-I";
  else oss << formatReal(im, cfg) << "I";

  return oss.str();
 }

 // std::string formatResult(const Value &v, const SystemConfig &cfg) {
 //  if (isReal(v)) return formatReal(std::get<double>(v), cfg);
 //  return formatComplex(std::get<std::complex<double>>(v), cfg);
 // }

 std::string formatResult(const Value &v, const SystemConfig &cfg) {
  return std::visit(
      [&](auto &&x) -> std::string {
       using T = std::decay_t<decltype(x)>;

       if constexpr (std::is_same_v<T, double>) {
        if (std::isinf(x)) return x > 0 ? "inf" : "-inf";
        if (std::isnan(x)) return "nan";
        return formatReal(x, cfg);

       } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        double re = x.real(), im = x.imag();

        if (std::isinf(re) || std::isinf(im)) {
         std::string s;
         if (std::isinf(re)) s += (re > 0 ? "inf" : "-inf");
         if (std::isinf(im)) {
          if (!s.empty() && im > 0) s += "+";
          s += (im > 0 ? "inf" : "-inf");
          s += "i";
         }
         return s;
        }

        return formatComplex(x, cfg);

       } else if constexpr (std::is_same_v<T, InvalidValue>) {
        return "Invalid";

       } else {
        static_assert(always_false<T>, "Unhandled Value type");
       }
      },
      v);
 }

 /* ============================
    ユーティリティの兄貴
    ============================ */
 Value Parser::ASTNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const { return Value(); }
 Value div(const Value &a, const Value &b, size_t pos) {
  if (isDouble(a) && isDouble(b)) {
   double r = asDouble(b, pos);
   if (r == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
   return asDouble(a) / r;
  }
  return asComplex(a) / asComplex(b);
 }

 /* ============================
   DLL
   ============================ */

 struct mmcal_context {
   SystemConfig cfg;
   std::vector<InputEntry> history;
   int base = 10;

   std::string last_error;
 };

 static int write_out(char *out, int cap, const std::string &s) {
  if (!out || cap <= 0) return -1;
  // cap-1 までコピーして末尾\0
  int n = (int)s.size();
  if (n >= cap) n = cap - 1;
  std::memcpy(out, s.data(), n);
  out[n] = '\0';
  return n;
 }

 extern "C" {

  MMCAL_API mmcal_context *mmcal_create(void) {
   try {
    auto *ctx = new mmcal_context();
    initFunctions(ctx->cfg); // あなたの関数登録
    return ctx;
   } catch (...) { return nullptr; }
  }

  MMCAL_API void mmcal_destroy(mmcal_context *ctx) { delete ctx; }

  MMCAL_API void mmcal_set_base(mmcal_context *ctx, int base) {
   if (!ctx) return;
   ctx->base = base;
  }

  MMCAL_API void mmcal_set_precision(mmcal_context *ctx, int precision) {
   if (!ctx) return;
   ctx->cfg.precision = precision;
  }

  MMCAL_API int mmcal_eval(mmcal_context *ctx, const char *expr, char *out, int out_cap) {
   if (!ctx) return 1;
   if (!expr) {
    ctx->last_error = "expr is null";
    write_out(out, out_cap, "ERROR: expr is null");
    return 2;
   }

   try {
    Value v = evaluate(expr, ctx->cfg, ctx->history, ctx->base);
    Value vr = roundValue(v, ctx->cfg);
    std::string s = formatResult(vr, ctx->cfg);

    ctx->history.push_back({expr, vr});
    ctx->last_error.clear();

    write_out(out, out_cap, s);
    return 0;

   } catch (const CalcError &e) {
    ctx->last_error = e.what();
    write_out(out, out_cap, std::string("ERROR: ") + e.what());
    return 10;

   } catch (const std::exception &e) {
    ctx->last_error = e.what();
    write_out(out, out_cap, std::string("ERROR: ") + e.what());
    return 11;

   } catch (...) {
    ctx->last_error = "unknown error";
    write_out(out, out_cap, "ERROR: unknown error");
    return 12;
   }
  }

  MMCAL_API int mmcal_last_error(mmcal_context *ctx, char *out, int out_cap) {
   if (!ctx) return 1;
   write_out(out, out_cap, ctx->last_error);
   return 0;
  }

 } // extern "C"
} // namespace mm::cal