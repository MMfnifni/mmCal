#include "dll_api.hpp"
#include "evaluate.hpp"
#include "formatter.hpp"
#include "functions.hpp"
#include "lexer_parser.hpp"

#include <cstring>
#include <string>
#include <unordered_map>
#include <vector>

namespace mm::cal {

 enum mmcal_status {
  MMCAL_OK = 0,
  MMCAL_E_NULL_CTX = -1,
  MMCAL_E_NULL_EXPR = -2,
  MMCAL_E_BUFFER_TOO_SMALL = -3,
  MMCAL_E_NULL_NAME = -4,
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
   std::unordered_map<std::string, Value> globals;
   std::unordered_map<std::string, UserFunction> userFunctions;
 };

 static int write_out(char *out, int cap, const std::string &s) {
  const int need = static_cast<int>(s.size()) + 1; // '\0' 含む

  if (!out || cap <= 0) return need;
  if (cap < need) return need;

  std::memcpy(out, s.data(), s.size());
  out[s.size()] = '\0';
  return need;
 }

 static void clear_error_state(mmcal_context *ctx) {
  ctx->last_error.clear();
  ctx->last_error_pos = -1;
 }

 static void set_error_state(mmcal_context *ctx, const std::string &msg, int pos) {
  ctx->last_error = msg;
  ctx->last_error_pos = pos;
 }

 static void reset_runtime_state(mmcal_context *ctx) {
  ctx->history.clear();
  ctx->globals.clear();
  ctx->userFunctions.clear();
  clear_error_state(ctx);
 }

 static int emit_error(mmcal_context *ctx, int status, const std::string &msg, int pos, char *err, int err_cap, int *err_need, int *err_pos) {
  set_error_state(ctx, msg, pos);

  const int need = static_cast<int>(ctx->last_error.size()) + 1;
  if (err_need) *err_need = need;
  if (err_pos) *err_pos = pos;

  if (!err || err_cap <= 0) return status;
  if (err_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;

  write_out(err, err_cap, ctx->last_error);
  return status;
 }

 static std::string build_success_output(const EvalResult &res, const SystemConfig &cfg) {
  std::string s;

  if (!res.explain.empty()) { s += res.explain; }

  if (res.suppressDisplay) {
   s += "evaluate success. display suppressed in silent mode";
   return s;
  }

  if (res.hasDisplayOverride()) {
   s += res.displayOverride;
   return s;
  }

  Value vr = roundValue(res.value, cfg);
  s += formatResult(vr, cfg);
  return s;
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

   const int need = static_cast<int>(ctx->last_error.size()) + 1;
   if (out_need) *out_need = need;

   if (out && out_cap > 0) {
    if (out_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;
    write_out(out, out_cap, ctx->last_error);
   }

   return MMCAL_OK;
  }

  MMCAL_API int mmcal_eval_ex(mmcal_context *ctx, const char *expr, char *out, int out_cap, int *out_need, char *err, int err_cap, int *err_need, int *err_pos) {
   if (!ctx) return MMCAL_E_NULL_CTX;

   if (out_need) *out_need = 0;
   if (err_need) *err_need = 0;
   if (err_pos) *err_pos = -1;

   clear_error_state(ctx);

   if (!expr) { return emit_error(ctx, MMCAL_E_NULL_EXPR, "expr is null", -1, err, err_cap, err_need, err_pos); }

   try {
    SessionState session{ctx->cfg, ctx->history, ctx->base, ctx->globals, ctx->userFunctions};
    EvaluationContext ectx{session};

    EvalResult res = evalLine(expr, ctx->cfg, ctx->history, ectx);
    applySideEffects(ectx, res);

    const std::string s = build_success_output(res, ctx->cfg);
    const int need = static_cast<int>(s.size()) + 1;

    if (out_need) *out_need = need;

    // In/Out 履歴は「表示文字列」ではなく値を持つ
    ctx->history.push_back({expr, res.value});

    if (!out || out_cap <= 0) return MMCAL_OK;
    if (out_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;

    write_out(out, out_cap, s);
    return MMCAL_OK;

   } catch (const ControlRequest &e) {
    switch (e.kind) {
     case ControlRequest::Kind::Exit: return MMCAL_S_EXIT;
     case ControlRequest::Kind::Clear: reset_runtime_state(ctx); return MMCAL_S_CLEARED;
     default: return emit_error(ctx, MMCAL_E_UNKNOWN, "unsupported control request in DLL mode", -1, err, err_cap, err_need, err_pos);
    }

   } catch (const CalcError &e) { return emit_error(ctx, MMCAL_E_CALC_ERROR, e.what(), static_cast<int>(e.pos), err, err_cap, err_need, err_pos); } catch (const std::exception &e) {
    return emit_error(ctx, MMCAL_E_STD_EXCEPTION, e.what(), -1, err, err_cap, err_need, err_pos);

   } catch (...) { return emit_error(ctx, MMCAL_E_UNKNOWN, "unknown error", -1, err, err_cap, err_need, err_pos); }
  }

  MMCAL_API int mmcal_syntax_check(mmcal_context *ctx, const char *expr, char *err, int err_cap, int *err_pos) {
   if (!ctx) return MMCAL_E_NULL_CTX;

   clear_error_state(ctx);
   if (err_pos) *err_pos = -1;

   if (!expr) {
    if (err_pos) *err_pos = -1;
    set_error_state(ctx, "expr is null", -1);

    const int need = static_cast<int>(ctx->last_error.size()) + 1;
    if (err && err_cap > 0) {
     if (err_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;
     write_out(err, err_cap, ctx->last_error);
    }
    return MMCAL_E_NULL_EXPR;
   }

   try {
    Parser p(ctx->cfg, expr);
    std::unique_ptr<ASTNode> ast = p.parse();

    (void)ast;

    if (p.cur.type != TokenType::End) { set_error_state(ctx, "Syntax Error: not closed", static_cast<int>(p.cur.pos)); }

   } catch (const CalcError &e) { set_error_state(ctx, e.what(), static_cast<int>(e.pos)); } catch (const std::exception &e) {
    set_error_state(ctx, e.what(), -1);

   } catch (...) { set_error_state(ctx, "unknown error", -1); }

   if (err_pos) *err_pos = ctx->last_error_pos;

   if (err && err_cap > 0) {
    const int need = static_cast<int>(ctx->last_error.size()) + 1;
    if (err_cap < need) return MMCAL_E_BUFFER_TOO_SMALL;
    write_out(err, err_cap, ctx->last_error);
   }

   return ctx->last_error_pos >= 0 ? MMCAL_E_CALC_ERROR : MMCAL_OK;
  }

  MMCAL_API int mmcal_unset_variable(mmcal_context *ctx, const char *name) {
   if (!ctx) return MMCAL_E_NULL_CTX;

   clear_error_state(ctx);

   if (!name) {
    set_error_state(ctx, "name is null", -1);
    return MMCAL_E_NULL_NAME;
   }

   ctx->globals.erase(name);
   return MMCAL_OK;
  }

  MMCAL_API int mmcal_undef_function(mmcal_context *ctx, const char *name) {
   if (!ctx) return MMCAL_E_NULL_CTX;

   clear_error_state(ctx);

   if (!name) {
    set_error_state(ctx, "name is null", -1);
    return MMCAL_E_NULL_NAME;
   }

   ctx->userFunctions.erase(name);
   return MMCAL_OK;
  }

  MMCAL_API int mmcal_clear_history(mmcal_context *ctx) {
   if (!ctx) return MMCAL_E_NULL_CTX;

   clear_error_state(ctx);
   ctx->history.clear();
   return MMCAL_OK;
  }

  MMCAL_API int mmcal_reset_session(mmcal_context *ctx) {
   if (!ctx) return MMCAL_E_NULL_CTX;

   reset_runtime_state(ctx);
   return MMCAL_OK;
  }

 } // extern "C"

} // namespace mm::cal