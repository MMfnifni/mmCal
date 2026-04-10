#include "evaluate.hpp"

namespace mm::cal {

 static constexpr int INDENT_WIDTH = 2;

 /* ============================
    評価マン(ToDo: Semanticの責務を兼ねている)
    ============================ */

 EvalResult evaluate(const std::string &src, SystemConfig &cfg, const std::vector<InputEntry> &hist, int base, EvaluationContext &ectx) {
  Parser p(cfg, src);
  auto ast = p.parse();

  if (p.cur.type != TokenType::End) throw CalcError(CalcErrorType::SyntaxError, "Syntax Error: not closed", p.cur.pos);
  Value v = ast->eval(ectx);
  return {v};
 }

 CalcResult calcEval(const std::string &expr, SystemConfig &cfg, std::vector<InputEntry> &history, int base, EvaluationContext &ectx) {
  try {
   EvalResult v = evaluate(expr, cfg, history, base, ectx);

   // 表示用に丸め
   Value vr = roundValue(v.value, cfg);
   std::string out = formatResult(vr, cfg);

   history.push_back({expr, vr});

   return {true, out, 0};

  } catch (const CalcError &e) { return {false, e.what(), e.pos}; } catch (const std::exception &e) {
   return {false, e.what(), 0};

  } catch (...) { return {false, "unknown error", 0}; }
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
   auto it = ectx.variables.find(n->name);
   if (it != ectx.variables.end()) {
    assignWasRedefinition = true;
    oldAssignedValue = it->second;
   }
  }

  Value v = ast->eval(ectx);

  if (auto *n = dynamic_cast<AssignNode *>(ast.get())) {
   std::string msg;

   if (assignWasRedefinition) {
    msg = "<redefined: " + n->name + " := " + formatResult(v, cfg) + " (was " + formatResult(oldAssignedValue, cfg) + ")>";
   } else {
    msg = "<defined: " + n->name + " := " + formatResult(v, cfg) + ">";
   }

   res.value = Value(msg);
   return res;
  }

  if (auto *n = dynamic_cast<FunctionDefNode *>(ast.get())) {
   std::string sig = n->name + "(";
   for (size_t i = 0; i < n->params.size(); ++i) {
    if (i > 0) sig += ", ";
    sig += n->params[i];
   }
   sig += ")";
   res.value = Value("<defined: " + sig + ">");
   return res;
  }

  res.value = v;
  return res;
 }

 // Value evalAssign(AssignNode *n, EvaluationContext &ctx) {
 //  if (ctx.userFunctions.contains(n->name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: name already used as function", 0);
 //  Value v = n->expr->eval(ctx);
 //  ctx.variables[n->name] = v;
 //  return v;
 // }

 // Value evalFunctionDef(FunctionDefNode *n, EvaluationContext &ctx) {
 //  if (ctx.variables.contains(n->name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: name already used as variable", 0);
 //  if (ctx.cfg.functions.contains(n->name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: cannot override builtin", 0);
 //  if (ctx.userFunctions.contains(n->name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: function already defined", 0);

 // UserFunction fn;
 // fn.params = n->params;
 // fn.body = cloneAST(n->body.get()); // 重要

 // ctx.userFunctions[n->name] = std::move(fn);

 // return Value();
 //}

 /* ============================
   フォーマッタ
   ============================ */

 inline static double normalizeZero(double x) {
  return (x == 0.0 ? 0.0 : x); // -0 → 0
 }

 inline static bool nearlyZero(double x, int precision) { return x < cnst_precision_inv && x > -cnst_precision_inv; }

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

  std::array<char, 512> buffer;
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

  // ---- 1パス構造判定 ----
  bool isVector = true;
  bool isColumnVector = true;
  bool allSingleRowMatrix = true;

  for (const auto &e : mv.elems_) {

   if (!e.isMulti()) {
    isColumnVector = false;
    allSingleRowMatrix = false;
   } else {
    const auto &inner = *e.asMulti(0);

    if (inner.elems_.size() != 1) isColumnVector = false;

    bool singleRow = true;
    for (const auto &row : inner.elems_) {
     if (row.isMulti()) {
      singleRow = false;
      break;
     }
    }
    if (!singleRow) allSingleRowMatrix = false;
   }

   if (e.isMulti()) isVector = false;
  }

  // ---- 横1行表示パターン ----
  if (isVector || isColumnVector || allSingleRowMatrix) {
   std::string out;
   out.reserve(mv.elems_.size() * 32);

   out += "{";
   for (size_t i = 0; i < mv.elems_.size(); ++i) {
    if (i > 0) out += ", ";
    out += formatResultMulti(mv.elems_[i], cfg, indent);
   }
   out += "}";

   return out;
  }

  // ---- 2D以上ブロック表示 ----
  std::string out;
  out.reserve(mv.elems_.size() * 64);

  std::string ind(indent, ' ');

  out += "{\n";

  for (size_t i = 0; i < mv.elems_.size(); ++i) {

   out += std::string(indent + INDENT_WIDTH, ' ');
   out += formatResultMulti(mv.elems_[i], cfg, indent + INDENT_WIDTH);

   if (i + 1 < mv.elems_.size()) out += ",";

   out += "\n";
  }

  out += ind + "}";

  return out;
 }

 inline void appendReal(double x, const SystemConfig &cfg, std::string &out) {
  if (std::isnan(x)) {
   out += "Indeterminate";
   return;
  }

  if (std::isinf(x)) {
   out += (x > 0 ? "Infinity" : "-Infinity");
   return;
  }

  std::array<char, 64> buffer;

  auto result = std::to_chars(buffer.data(), buffer.data() + buffer.size(), x, std::chars_format::general, cfg.precision);

  out.append(buffer.data(), result.ptr);
 }

 inline void appendComplex(const std::complex<double> &c, const SystemConfig &cfg, std::string &out) {
  double re = (c.real() == 0.0 ? 0.0 : c.real());
  double im = (c.imag() == 0.0 ? 0.0 : c.imag());

  const bool re0 = (re < cnst_precision_inv && re > -cnst_precision_inv);
  const bool im0 = (im < cnst_precision_inv && im > -cnst_precision_inv);

  // 虚部ゼロ → 実数扱い
  if (im0) {
   appendReal(re, cfg, out);
   return;
  }

  // ---- 実部 ----
  if (!re0) { appendReal(re, cfg, out); }

  // ---- 虚部 ----
  const bool neg = im < 0.0;
  const double absIm = neg ? -im : im;

  if (!re0) {
   out += (neg ? "-" : "+");
  } else if (neg) {
   out += "-";
  }

  // ±1I のときは数値を出さない
  if (!(absIm < 1.0 + cnst_precision_inv && absIm > 1.0 - cnst_precision_inv)) { appendReal(absIm, cfg, out); }

  out += "I";
 }

 inline void appendMulti(const MultiValue &mv, const SystemConfig &cfg, std::string &out) {
  out.push_back('{');

  for (size_t i = 0; i < mv.elems_.size(); ++i) {
   if (i > 0) {
    out.push_back(',');
    out.push_back(' ');
   }

   appendValue(mv.elems_[i], cfg, out);
  }

  out.push_back('}');
 }

 void appendValue(const Value &v, const SystemConfig &cfg, std::string &out) {
  std::visit(
      [&](auto &&x) {
       using T = std::decay_t<decltype(x)>;

       if constexpr (std::is_same_v<T, double>) {
        appendReal(x, cfg, out);
       } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        appendComplex(x, cfg, out);
       } else if constexpr (std::is_same_v<T, std::shared_ptr<MultiValue>>) {
        appendMulti(*x, cfg, out);
       } else if constexpr (std::is_same_v<T, std::string>) {
        out += x;
       } else if constexpr (std::is_same_v<T, std::monostate>) {
        // 関数定義などの「値なし」は何も表示しない
       } else {
        out += "Invalid";
       }
      },
      v.storage());
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

 std::string formatResultMulti(const Value &v, const SystemConfig &cfg, int indent) {
  std::string r, chk;
  r = std::visit(
      [&](auto &&x) -> std::string {
       using T = std::decay_t<decltype(x)>;

       if constexpr (std::is_same_v<T, std::monostate>) {
        return "";
       } else if constexpr (std::is_same_v<T, std::shared_ptr<MultiValue>>) {
        return formatMultiValue(*x, cfg, indent);
       } else if constexpr (std::is_same_v<T, double>) {
        return formatReal(x, cfg);
       } else if constexpr (std::is_same_v<T, std::complex<double>>) {
        return formatComplex(x, cfg);
       } else if constexpr (std::is_same_v<T, std::string>) {
        return x;
       } else {
        throw CalcError(CalcErrorType::RuntimeError, "Unknown Value type", 0);
       }
      },
      v.storage());
  return r;
 }

 std::string formatResult(const Value &v, const SystemConfig &cfg) {
  std::string out;
  out.reserve(256);

  appendValue(v, cfg, out);

  return out;
 }

 std::string dataExplainStructure(const Value &value) {
  std::function<std::string(const Value &)> dump = [&](const Value &v) -> std::string {
   if (v.isMulti()) {
    auto mv = v.asMulti(0);

    std::string out = "{";

    for (size_t i = 0; i < mv->size(); ++i) {
     if (i > 0) out += ",";

     out += dump((*mv)[i]);
    }

    out += "}";
    return out;
   }

   if (v.isScalar()) return "R";
   if (v.isComplex()) return "C";
   if (v.isString()) return "S";

   return "I";
  };

  return "Structural Overview: " + dump(value);
 }

 std::string dataExplain(Value value, const SystemConfig &cfg) {
  std::ostringstream oss;

  oss << "\n[Data Explain]\n\n";

  std::string displayed = formatResult(value, cfg);

  oss << "Type = ";
  if (value.isScalar()) oss << "Real Number\n";
  else if (value.isComplex()) oss << "Complex Number\n";
  else if (value.isString()) oss << "String\n";
  else if (value.isMulti()) oss << "MultiValue\n";
  else oss << "Invalid\n";

  oss << "Displayed = " << displayed << "\n";

  // 数値系のみ誤差解析
  if (value.isNumeric()) {
   std::complex<double> z = value.toComplex(0);

   double re = z.real();
   double im = z.imag();

   re = normalizeZero(re);
   im = normalizeZero(im);

   oss << "Raw = ";

   if (im == 0.0) oss << formatReal(re, cfg);
   else oss << formatComplex(z, cfg);

   oss << "\n";

   double threshold = cnst_precision_inv;

   double errMag = std::abs(z);

   oss << "Magnitude = " << formatReal(errMag, cfg) << "\n";

   oss << "Noise threshold = " << formatReal(threshold, cfg) << "\n";

   // if (errMag < threshold) oss << "Status = Numerical noise level\n";
   // else oss << "Status = Significant magnitude\n";

  } else if (value.isMulti()) {
   auto mv = value.asMulti(0);
   oss << "MultiValue = true\n";

   // 再帰的に構造表示
   std::function<void(const MultiValue &, int)> dumpStructure = [&](const MultiValue &m, int depth) {
    oss << std::string(depth * 2, ' ') << "MultiValue(depth=" << depth << ") "
        << "size=" << m.size() << "\n";

    size_t index = 0;
    for (const auto &elem : m.elems_) {

     oss << std::string(depth * 2 + 2, ' ') << "[" << index++ << "] = ";

     if (elem.isMulti()) {
      oss << "\n";
      dumpStructure(*elem.asMulti(0), depth + 1);
     } else if (elem.isScalar()) {
      oss << formatReal(elem.asScalar(0), cfg) << "\n";
     } else if (elem.isComplex()) {
      oss << formatComplex(elem.asComplex(0), cfg) << "\n";
     } else if (elem.isString()) {
      oss << elem.asString(0) << "\n";
     } else {
      oss << "Invalid\n";
     }
    }
   };

   dumpStructure(*mv, 0);
  }
  oss << dataExplainStructure(value) << "\n\n[Data Explain END]\n\n";
  return oss.str();
 }

 /* ============================
    ユーティリティの兄貴
    ============================ */
 Value evalImpl(EvaluationContext &ectx) { return Value(); }

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
    EvaluationContext ectx{ctx->cfg, ctx->history, ctx->base};
    EvalResult v = evaluate(expr, ctx->cfg, ctx->history, ctx->base, ectx);
    Value vr = roundValue(v.value, ctx->cfg);
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

   } catch (ControlRequest &e) {
    switch (e.kind) {
     case ControlRequest::Kind::Exit: return MMCAL_S_EXIT;
     case ControlRequest::Kind::Clear:
      ctx->history.clear();
      ctx->last_error.clear();
      ctx->last_error_pos = -1;
      return MMCAL_S_CLEARED;
     case ControlRequest::Kind::Explain:
     case ControlRequest::Kind::FileWrite:
     case ControlRequest::Kind::ClipboardCopy:
     default: return MMCAL_E_UNKNOWN;
    }
   } catch (const std::exception &e) {
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
    std::unique_ptr<ASTNode> ast = p.parseCompare();

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