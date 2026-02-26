#include "evaluate.hpp"

namespace mm::cal {

 static constexpr int INDENT_WIDTH = 2;

 /* ============================
    評価マン
    ============================ */

 Value evaluate(const std::string &src, SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) {
  Parser p(cfg, src);
  std::unique_ptr<Parser::ASTNode> ast = p.parseCompare();
  if (p.cur.type != TokenType::End) throw CalcError(CalcErrorType::SyntaxError, errorMessage(CalcErrorType::SyntaxError), p.cur.pos);
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

 inline static double normalizeZero(double x) {
  return (x == 0.0 ? 0.0 : x); // -0 → 0
 }

 inline static double epsilon(int precision) { return std::pow(10.0, -precision); }

 inline static bool nearlyZero(double x, int precision) { return std::abs(x) < epsilon(precision); }

 inline static void trimZeros(std::string &s) {
  if (auto pos = s.find('.'); pos != std::string::npos) {
   while (!s.empty() && s.back() == '0')
    s.pop_back();
   if (!s.empty() && s.back() == '.') s.pop_back();
  }
 }

 std::string formatReal(double x, const SystemConfig &cfg) {
  if (std::isnan(x)) return "Indeterminate";

  if (std::isinf(x)) return x > 0 ? "Infinity" : "-Infinity";

  x = normalizeZero(x);

  std::array<char, 64> buffer;
  auto *begin = buffer.data();
  auto *end = buffer.data() + buffer.size();

  std::chars_format fmt = cfg.trimTrailingZeros ? std::chars_format::fixed : std::chars_format::general;

  auto result = std::to_chars(begin, end, x, fmt, cfg.precision);

  if (result.ec != std::errc{}) return "0"; // fallback（ほぼ起きない）

  std::string s(begin, result.ptr);

  if (cfg.trimTrailingZeros) trimZeros(s);

  return s;
 }

 std::string formatComplex(const std::complex<double> &c, const SystemConfig &cfg) {
  double re = normalizeZero(c.real());
  double im = normalizeZero(c.imag());

  const bool re0 = nearlyZero(re, cfg.precision);
  const bool im0 = nearlyZero(im, cfg.precision);

  if (im0) return formatReal(re, cfg);

  std::string out;

  // ---- 実部 ----
  if (!re0) out += formatReal(re, cfg);

  // ---- 虚部 ----
  auto appendImag = [&](double value) {
   const bool neg = value < 0;
   const double absIm = std::abs(value);

   if (!out.empty()) out += (neg ? "-" : "+");
   else if (neg) out += "-";

   // ±1I のとき数値を出さない
   if (!nearlyZero(absIm - 1.0, cfg.precision)) out += formatReal(absIm, cfg);

   out += "I";
  };

  appendImag(im);

  return out;
 }

 std::string formatMultiValue(const MultiValue &mv, const SystemConfig &cfg, int indent) {
  if (mv.elems_.empty()) return "{}";

  // ----- 1D判定 -----
  bool isVector = std::all_of(mv.elems_.begin(), mv.elems_.end(), [](const Value &v) { return !v.isMulti(); });
  // ----- n×1 行列判定（縦ベクトル） -----
  bool isColumnVector = true;

  if (!mv.elems_.empty()) {
   for (const auto &e : mv.elems_) {
    if (!e.isMulti()) {
     isColumnVector = false;
     break;
    }

    const auto &inner = *e.asMulti(0);
    if (inner.elems_.size() != 1) {
     isColumnVector = false;
     break;
    }
   }
  } else {
   isColumnVector = false;
  }

  if (isColumnVector) {
   std::string out = "{";
   for (size_t i = 0; i < mv.elems_.size(); ++i) {
    if (i > 0) out += ", ";
    out += formatResultMulti(mv.elems_[i], cfg, indent);
   }
   out += "}";
   return out;
  }
  if (isVector) {
   std::string out = "{";
   for (size_t i = 0; i < mv.elems_.size(); ++i) {
    if (i > 0) out += ", ";
    out += formatResultMulti(mv.elems_[i], cfg, indent);
   }
   out += "}";
   return out;
  }

  // ----- すべてが「1行行列」なら横圧縮 -----
  bool allSingleRowMatrix = true;

  for (const auto &e : mv.elems_) {

   if (!e.isMulti()) {
    allSingleRowMatrix = false;
    break;
   }

   const auto &inner = *e.asMulti(0);

   // 行列であること
   if (inner.elems_.empty()) {
    allSingleRowMatrix = false;
    break;
   }

   // 1行のみか？
   bool isSingleRow = true;
   for (const auto &row : inner.elems_) {
    if (row.isMulti()) {
     isSingleRow = false;
     break;
    }
   }

   if (!isSingleRow) {
    allSingleRowMatrix = false;
    break;
   }
  }

  if (allSingleRowMatrix) {
   std::string out = "{";
   for (size_t i = 0; i < mv.elems_.size(); ++i) {
    if (i > 0) out += ", ";
    out += formatResultMulti(mv.elems_[i], cfg, indent);
   }
   out += "}";
   return out;
  }

  // ----- 2D以上 -----
  std::string ind(indent, ' ');
  std::string out = "{\n";

  for (size_t i = 0; i < mv.elems_.size(); ++i) {

   out += std::string(indent + INDENT_WIDTH, ' ');
   out += formatResultMulti(mv.elems_[i], cfg, indent + INDENT_WIDTH);

   if (i + 1 < mv.elems_.size()) out += ",";

   out += "\n";
  }

  out += ind + "}";

  return out;
 }

 // std::string formatResult(const Value &v, const SystemConfig &cfg) {
 //  if (isReal(v)) return formatReal(std::get<double>(v), cfg);
 //  return formatComplex(std::get<std::complex<double>>(v), cfg);
 // }
 static inline bool containsMultiValue(const MultiValue &mv) {
  for (const auto &e : mv.elems_) {
   if (std::holds_alternative<std::shared_ptr<MultiValue>>(e.storage())) { return true; }
  }
  return false;
 }
 static inline std::string indentStr(int indent) { return std::string(indent * 2, ' '); }

 static inline int depthOf(const Value &v) {
  if (!std::holds_alternative<std::shared_ptr<MultiValue>>(v.storage())) return 0;

  const auto &mv = *std::get<std::shared_ptr<MultiValue>>(v.storage());

  int maxChild = 0;
  for (const auto &e : mv.elems_) {
   maxChild = std::max(maxChild, depthOf(e));
  }
  return 1 + maxChild;
 }

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

       } else if constexpr (std::is_same_v<T, std::string>) {
        return x;

       } else if constexpr (std::is_same_v<T, InvalidValue>) {
        return "Invalid";

       } else if constexpr (std::is_same_v<T, std::shared_ptr<MultiValue>>) {

        std::ostringstream oss;
        oss << "{";

        for (size_t i = 0; i < x->elems_.size(); ++i) {
         if (i > 0) oss << ", ";
         oss << formatResult(x->elems_[i], cfg);
        }
        oss << "}";
        return oss.str();

       } else if constexpr (std::is_same_v<T, std::monostate>) {
        return "Invalid";
       } else {
        static_assert(always_false<T>, "Unhandled Value type");
       }
      },
      v.storage());
 }

 std::string formatResultMulti(const Value &v, const SystemConfig &cfg, int indent) {
  return std::visit(
      [&](auto &&x) -> std::string {
       using T = std::decay_t<decltype(x)>;

       if constexpr (std::is_same_v<T, std::shared_ptr<MultiValue>>) {
        return formatMultiValue(*x, cfg, indent);
       } else if constexpr (std::is_same_v<T, double>) {
        return formatReal(x, cfg);
       } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        return formatComplex(x, cfg);
       } else if constexpr (std::is_same_v<T, std::string>) {
        return x;
       } else if constexpr (std::is_same_v<T, InvalidValue>) {
        return "Invalid";
       } else {
        throw std::runtime_error("Unknown Value type");
       }
      },
      v.storage());
 }

 /* ============================
    ユーティリティの兄貴
    ============================ */
 Value Parser::ASTNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const { return Value(); }

 /* ============================
   DLL
   ============================ */

 enum mmcal_status {
  MMCAL_OK = 0,
  MMCAL_E_NULL_CTX = -1,
  MMCAL_E_NULL_EXPR = -2,
  MMCAL_E_BUFFER_TOO_SMALL = -3,
  MMCAL_E_CALC_ERROR = -10,
  MMCAL_E_STD_EXCEPTION = -11,
  MMCAL_E_UNKNOWN = -12,
  MMCAL_S_CLEARED = -13,
  MMCAL_S_EXIT = -14,
 };

 struct mmcal_context {
   SystemConfig cfg;
   std::vector<InputEntry> history;
   int base = 10;
   int last_error_pos = -1;

   std::string last_error;
 };

 static int write_out(char *out, int cap, const std::string &s) {
  const int need = (int)s.size() + 1; // '\0' 含む

  if (!out || cap <= 0) return need; // バッファ無しなら、必要サイズだけ返す
  if (cap < need) return need;       // 足りないなら何も書かない（呼び出し側が out_need を見て再確保する）

  std::memcpy(out, s.data(), s.size());
  out[s.size()] = '\0';
  return need;
 }

 extern "C" {

  MMCAL_API mmcal_context *mmcal_create(void) {
   try {
    auto *ctx = new mmcal_context();
    initFunctions(ctx->cfg);
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

  MMCAL_API int mmcal_last_error_pos(mmcal_context *ctx) {
   if (!ctx) return -1;
   return ctx->last_error_pos;
  }

  MMCAL_API int mmcal_last_error(mmcal_context *ctx, char *out, int out_cap, int *out_need) {
   if (!ctx) return MMCAL_E_NULL_CTX;

   const int need = (int)ctx->last_error.size() + 1;
   if (out_need) *out_need = need;

   write_out(out, out_cap, ctx->last_error);
   return MMCAL_OK;
  }

  MMCAL_API int mmcal_eval_ex(mmcal_context *ctx, const char *expr, char *out, int out_cap, int *out_need, char *err, int err_cap, int *err_need, int *err_pos) {
   if (!ctx) return MMCAL_E_NULL_CTX;

   if (out_need) *out_need = 0;
   if (err_need) *err_need = 0;
   if (err_pos) *err_pos = -1;

   ctx->last_error.clear();
   ctx->last_error_pos = -1;

   if (!expr) {
    ctx->last_error = "expr is null";
    ctx->last_error_pos = -1;

    const int need = (int)ctx->last_error.size() + 1;
    if (err_need) *err_need = need;

    write_out(err, err_cap, ctx->last_error);
    return MMCAL_E_NULL_EXPR;
   }

   try {
    Value v = evaluate(expr, ctx->cfg, ctx->history, ctx->base);
    Value vr = roundValue(v, ctx->cfg);
    std::string s = formatResult(vr, ctx->cfg);

    ctx->history.push_back({expr, vr});

    const int need = (int)s.size() + 1;
    if (out_need) *out_need = need;

    // out が無い/足りないならコピーせず、ステータスで返す
    if (!out || out_cap <= 0) return MMCAL_OK;
    if (out_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;

    write_out(out, out_cap, s);
    return MMCAL_OK;

   } catch (const CalcError &e) {
    ctx->last_error = e.what();
    ctx->last_error_pos = (int)e.pos;

    const int need = (int)ctx->last_error.size() + 1;
    if (err_need) *err_need = need;
    if (err_pos) *err_pos = (int)e.pos;

    if (!err || err_cap <= 0) return MMCAL_E_CALC_ERROR;
    if (err_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;

    write_out(err, err_cap, ctx->last_error);
    return MMCAL_E_CALC_ERROR;

   } catch (const ClearRequest &) {
    ctx->history.clear();
    ctx->last_error.clear();
    ctx->last_error_pos = -1;
    return MMCAL_S_CLEARED;
   } catch (const ExitRequest &) { return MMCAL_S_EXIT; } catch (const std::exception &e) {
    ctx->last_error = e.what();
    ctx->last_error_pos = -1;

    const int need = (int)ctx->last_error.size() + 1;
    if (err_need) *err_need = need;
    if (err_pos) *err_pos = -1;

    if (!err || err_cap <= 0) return MMCAL_E_STD_EXCEPTION;
    if (err_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;

    write_out(err, err_cap, ctx->last_error);
    return MMCAL_E_STD_EXCEPTION;

   } catch (...) {
    ctx->last_error = "unknown error";
    ctx->last_error_pos = -1;

    const int need = (int)ctx->last_error.size() + 1;
    if (err_need) *err_need = need;
    if (err_pos) *err_pos = -1;

    if (!err || err_cap <= 0) return MMCAL_E_UNKNOWN;
    if (err_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;

    write_out(err, err_cap, ctx->last_error);
    return MMCAL_E_UNKNOWN;
   }
  }

  // expr の構文だけ解析してエラー位置を返す
  MMCAL_API int mmcal_syntax_check(mmcal_context *ctx, const char *expr, char *err, int err_cap, int *err_pos) {
   if (!ctx) return MMCAL_E_NULL_CTX;
   if (!expr) return MMCAL_E_NULL_EXPR;

   ctx->last_error.clear();
   ctx->last_error_pos = -1;

   try {
    Parser p(ctx->cfg, expr);
    std::unique_ptr<Parser::ASTNode> ast = p.parseCompare();

    // 未消化トークンがあれば構文エラー
    if (p.cur.type != TokenType::End) {
     ctx->last_error = "Syntax error";
     ctx->last_error_pos = static_cast<int>(p.cur.pos);
    }

   } catch (const CalcError &e) {
    ctx->last_error = e.what();
    ctx->last_error_pos = (int)e.pos;
   } catch (const std::exception &e) {
    ctx->last_error = e.what();
    ctx->last_error_pos = -1;
   } catch (...) {
    ctx->last_error = "unknown error";
    ctx->last_error_pos = -1;
   }

   if (err_pos) *err_pos = ctx->last_error_pos;
   if (err && err_cap > 0) write_out(err, err_cap, ctx->last_error);

   return ctx->last_error_pos >= 0 ? MMCAL_E_CALC_ERROR : MMCAL_OK;
  }

 } // extern "C"
} // namespace mm::cal