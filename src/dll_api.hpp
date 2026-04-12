#pragma once

#include "core.hpp"

namespace mm::cal {

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

  MMCAL_API int mmcal_last_error_pos(mmcal_context *ctx);

  MMCAL_API int mmcal_eval_ex(mmcal_context *ctx, const char *expr, char *out, int out_cap, int *out_need, char *err, int err_cap, int *err_need, int *err_pos);

  MMCAL_API int mmcal_syntax_check(mmcal_context *ctx, const char *expr, char *err, int err_cap, int *err_pos);

  MMCAL_API int mmcal_last_error(mmcal_context *ctx, char *out, int out_cap, int *out_need);

  MMCAL_API int mmcal_unset_variable(mmcal_context *ctx, const char *name);
  MMCAL_API int mmcal_undef_function(mmcal_context *ctx, const char *name);
  MMCAL_API int mmcal_clear_history(mmcal_context *ctx);
  MMCAL_API int mmcal_reset_session(mmcal_context *ctx);

#ifdef __cplusplus
 }
#endif
} // namespace mm::cal