#pragma once

#include "core.hpp"
#include "syntax.hpp"
#include <iomanip>
#include <sstream>
#include <string>

namespace mm::cal {

 /* ============================
    評価マン
    ============================ */

 Value evaluate(const std::string &src, SystemConfig &cfg, const std::vector<InputEntry> &hist, int base); // sytaxにあるよ
 CalcResult calcEval(const std::string &expr, SystemConfig &cfg, std::vector<InputEntry> &history, int base);
 EvalResult evalLine(const std::string &line, SystemConfig &cfg, RuntimeState &rt, std::vector<InputEntry> &history);

 /* ============================
   フォーマッタ
   ============================ */

 std::string formatReal(double x, const SystemConfig &cfg);
 std::string formatComplex(const std::complex<double> &c, const SystemConfig &cfg);
 std::string formatResult(const Value &v, const SystemConfig &cfg);

 /* ============================
  表示
  ============================ */
 inline void setupStream(std::ostringstream &oss, const SystemConfig &cfg) {
  if (cfg.precision >= 0) {
   oss << std::fixed << std::setprecision(cfg.precision);
  } else {
   oss << std::setprecision(15);
  }
 }

 /* ============================
    単位変換のお姉さん
    ============================ */
 inline Value deg2rad_v(const Value &v, size_t pos) {
  constexpr double k = PI / 180.0;
  return mulReal(v, k, pos);
 }

 inline Value rad2deg_v(const Value &v, size_t pos) {
  constexpr double k = 180.0 / PI;
  return mulReal(v, k, pos);
 }

 inline Value deg2grad_v(const Value &v, size_t pos) { return mulReal(v, 400.0 / 360.0, pos); }

 inline Value grad2deg_v(const Value &v, size_t pos) { return mulReal(v, 360.0 / 400.0, pos); }

 inline Value rad2grad_v(const Value &v, size_t pos) {
  constexpr double k = 200.0 / PI;
  return mulReal(v, k, pos);
 }

 inline Value grad2rad_v(const Value &v, size_t pos) {
  constexpr double k = PI / 200.0;
  return mulReal(v, k, pos);
 }

 /* ============================
    ユーティリティの兄貴
    ============================ */
 Value div(const Value &a, const Value &b, size_t pos);

 // Value evalCompare(const Value &lhs, const Value &rhs, Parser::CmpOp op, FunctionContext &ctx);//syntaxにあるよ

 /* ============================
   DLL
   ============================ */
#ifdef __cplusplus
 extern "C" {
#endif

#ifdef _WIN32
# define MMCAL_API __declspec(dllexport)
#else
# define MMCAL_API
#endif

  typedef struct mmcal_context mmcal_context;

  MMCAL_API mmcal_context *mmcal_create(void);
  MMCAL_API void mmcal_destroy(mmcal_context *ctx);

  MMCAL_API void mmcal_set_base(mmcal_context *ctx, int base);
  MMCAL_API void mmcal_set_precision(mmcal_context *ctx, int precision);

  MMCAL_API int mmcal_eval_ex(mmcal_context *ctx, const char *expr, char *out, int out_cap, int *out_need, char *err, int err_cap, int *err_need, int *err_pos);

  MMCAL_API int mmcal_last_error(mmcal_context *ctx, char *out, int out_cap, int *out_need);

#ifdef __cplusplus
 }
#endif

} // namespace mm::cal
