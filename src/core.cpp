#include "core.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

namespace mm::cal {

 Value roundValue(const Value &v, const SystemConfig &cfg) {
  auto roundDouble = [&](double d) { return std::round(d / cnst_precision_inv) * cnst_precision_inv; };

  return std::visit(Overloaded{// double
                               [&](double d) -> Value { return roundDouble(d); },

                               // complex
                               [&](std::complex<double> c) -> Value { return std::complex<double>(roundDouble(c.real()), roundDouble(c.imag())); },

                               // multi (vectorized)
                               [&](std::shared_ptr<MultiValue> m) -> Value {
                                auto result = std::make_shared<MultiValue>();
                                result->elems_.reserve(m->elems_.size());

                                for (const auto &e : m->elems_)
                                 result->elems_.push_back(roundValue(e, cfg));

                                return result;
                               },

                               // string / invalid
                               [&](auto &&x) -> Value {
                                if constexpr (std::is_same_v<std::decay_t<decltype(x)>, std::monostate>) {
                                 return Value();
                                } else if constexpr (std::is_same_v<std::decay_t<decltype(x)>, std::string>) {
                                 return Value(x);
                                } else {
                                 return Value(x); // double / complex / MultiPtr だけ
                                }
                               }},
                    v.storage());
 }

 bool Value::hasNestedMulti(size_t pos) const {
  const auto &mv = asMultiRef(pos);
  for (const auto &v : mv.elems()) {
   if (v.isMulti()) return true;
  }
  return false;
 }

 Value::Value(const EvalResult &r) : data_(std::monostate{}) { *this = r.value; }

 // multi util
 inline const MultiValue &Value::asMultiRef(size_t pos) const { return *asMulti(pos); }
 std::size_t Value::multiSize(size_t pos) const { return asMultiRef(pos).size(); }
 bool Value::multiEmpty(size_t pos) const { return asMultiRef(pos).size() == 0; }

 const Value &Value::multiAt(std::size_t i, size_t pos) const {
  const auto &mv = asMultiRef(pos);
  if (i >= mv.size()) throw CalcError(CalcErrorType::OutOfRange, "OutOfRange: multivalue index out of range", pos);
  return mv[i];
 }

 bool Value::isEmptyMulti() const noexcept {
  if (!isMulti()) return false;
  return std::get<MultiPtr>(data_)->size() == 0;
 }

 void Value::checkFiniteImpl(const MultiPtr &mv, size_t pos) {
  if (!mv) return;
  for (const auto &e : mv->elems_)
   Value::checkFinite(e, pos);
 }

 template <class F> Value MultiValue::map(F &&f) const {
  std::vector<Value> res;
  res.reserve(elems_.size());

  for (const auto &v : elems_)
   res.push_back(f(v));

  return Value(std::make_shared<MultiValue>(std::move(res)));
 }

 bool MultiValue::hasNested() const noexcept {
  for (const auto &v : elems_)
   if (v.isMulti()) return true;
  return false;
 }
 std::vector<std::vector<double>> MultiValue::toMatrix() const {
  std::vector<std::vector<double>> matrix;
  for (const auto &row : elems_) {
   if (!row.isMulti()) { throw CalcError(CalcErrorType::TypeError, "TypeError: Invalid matrix row", 0); }
   std::vector<double> row_vals;
   const auto &multi_row = row.asMultiRef(0);
   for (size_t i = 0; i < multi_row.size(); ++i) {
    row_vals.push_back(multi_row[i].asScalar(0));
   }
   matrix.push_back(std::move(row_vals));
  }
  return matrix;
 }

 std::vector<std::vector<double>> toMatrix(const Value &val, const FunctionContext &ctx) {
  const auto &M = val.asMultiRef(ctx.pos);

  size_t rows = M.size();
  if (rows == 0) throw CalcError(CalcErrorType::DomainError, "DomainError: empty matrix", ctx.pos);

  size_t cols = M[0].asMultiRef(ctx.pos).size();

  std::vector<std::vector<double>> A(rows, std::vector<double>(cols));

  for (size_t i = 0; i < rows; ++i) {
   const auto &row = M[i].asMultiRef(ctx.pos);

   if (row.size() != cols) throw CalcError(CalcErrorType::DomainError, "DomainError: irregular matrix", ctx.pos);

   for (size_t j = 0; j < cols; ++j)
    A[i][j] = row[j].asScalar(ctx.pos);
  }

  return A;
 }

 std::string whatByteUnit(size_t filesize) {
  const char *suffix[] = {"B", "KB", "MB", "GB", "TB"};
  size_t unit_index = 0;
  while (filesize > 1024 && unit_index < 4) {
   filesize /= 1024;
   ++unit_index;
  }
  return std::to_string(static_cast<int>(filesize)).append(suffix[unit_index]);
 }

 void serializeValueImpl(const Value &v, std::ostream &os, const SystemConfig &cfg, int depth = 0) {
  auto formatReal = [&](double x) -> std::string {
   if (!std::isfinite(x)) return "NaN";

   if (std::abs(x) < cnst_precision_inv) x = 0.0;

   std::ostringstream ss;
   ss << std::fixed << std::setprecision(cfg.precision) << x;

   std::string s = ss.str();

   if (cfg.trimTrailingZeros) {
    auto pos = s.find('.');
    if (pos != std::string::npos) {
     size_t end = s.size() - 1;

     while (end > pos && s[end] == '0')
      --end;

     if (s[end] == '.') --end;

     s.resize(end + 1);
    }
   }

   return s;
  };

  std::visit(
      [&](auto &&x) {
       using X = std::decay_t<decltype(x)>;

       if constexpr (std::is_same_v<X, std::monostate>) {
        os << "Invalid";
       }

       else if constexpr (std::is_same_v<X, double>) {
        os << formatReal(x);
       }

       else if constexpr (std::is_same_v<X, std::complex<double>>) {
        double re = x.real();
        double im = x.imag();

        bool hasRe = std::abs(re) > cnst_precision_inv;
        bool hasIm = std::abs(im) > cnst_precision_inv;

        if (!hasRe && !hasIm) {
         os << "0";
         return;
        }

        if (hasRe) os << formatReal(re);

        if (hasIm) {
         if (im > 0 && hasRe) os << "+";

         if (std::abs(im - 1.0) < cnst_precision_inv) os << "I";
         else if (std::abs(im + 1.0) < cnst_precision_inv) os << "-I";
         else os << formatReal(im) << "I";
        }
       }

       else if constexpr (std::is_same_v<X, std::string>) {
        os << "\"" << x << "\"";
       }

       else if constexpr (std::is_same_v<X, Value::MultiPtr>) {
        if (!x) {
         os << "{}";
         return;
        }

        os << "{";

        const auto &elems = x->elems();

        for (size_t i = 0; i < elems.size(); ++i) {
         if (i > 0) os << ",";

         serializeValueImpl(elems[i], os, cfg, depth + 1);
        }

        os << "}";
       }

       else {
        os << "Unknown";
       }
      },
      v.storage());
 }

 void applySideEffects(EvaluationContext &ectx, EvalResult &result) {
  for (const auto &e : ectx.sideEffects) {
   switch (e.kind) {
    case SideEffect::Kind::Explain:
     {
      result.explain += e.a;
      break;
     }

    case SideEffect::Kind::SuppressDisplay:
     {
      result.suppressDisplay = true;
      break;
     }

    case SideEffect::Kind::Message:
     {
      if (!result.explain.empty() && result.explain.back() != '\n') { result.explain += '\n'; }
      result.explain += e.a;
      if (!result.explain.empty() && result.explain.back() != '\n') { result.explain += '\n'; }
      break;
     }

    case SideEffect::Kind::FileWrite:
     {
      namespace fs = std::filesystem;

      const std::string &pathUtf8 = e.a;
      const std::string &content = e.b;

#ifdef _WIN32
      fs::path path(std::u8string_view(reinterpret_cast<const char8_t *>(pathUtf8.data()), pathUtf8.size()));
#else
      fs::path path(pathUtf8);
#endif

      fs::path tmp = path;
      tmp += ".tmp";

      try {
       std::ofstream ofs(tmp, std::ios::binary | std::ios::trunc);
       if (!ofs) { throw CalcError(CalcErrorType::IOError, "IOError: failed to open temp file", 0); }

       ofs.write(content.data(), static_cast<std::streamsize>(content.size()));
       if (!ofs) { throw CalcError(CalcErrorType::IOError, "IOError: failed to write file", 0); }

       ofs.close();
       if (!ofs) { throw CalcError(CalcErrorType::IOError, "IOError: failed to close file", 0); }

       std::error_code ec;
       fs::rename(tmp, path, ec);
       if (ec) {
        fs::remove(path, ec);
        ec.clear();
        fs::rename(tmp, path, ec);
        if (ec) {
         fs::remove(tmp, ec);
         throw CalcError(CalcErrorType::IOError, "IOError: failed to replace target file", 0);
        }
       }

       if (!result.explain.empty() && result.explain.back() != '\n') { result.explain += '\n'; }
       result.explain += "file write completed (" + whatByteUnit(content.size()) + ")\n";

      } catch (...) {
       std::error_code ec;
       fs::remove(tmp, ec);
       throw;
      }

      break;
     }

    case SideEffect::Kind::ClipboardCopy:
     {
#ifdef _WIN32
      const std::string &s = e.a;

      if (!OpenClipboard(nullptr)) { throw CalcError(CalcErrorType::IOError, "OpenClipboard failed", 0); }

      try {
       if (!EmptyClipboard()) { throw CalcError(CalcErrorType::InternalError, "EmptyClipboard failed", 0); }

       int size = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, nullptr, 0);
       if (size <= 0) { throw CalcError(CalcErrorType::InternalError, "UTF8->UTF16 convert failed", 0); }

       HGLOBAL hMem = GlobalAlloc(GMEM_MOVEABLE, size * sizeof(wchar_t));
       if (!hMem) { throw CalcError(CalcErrorType::InternalError, "GlobalAlloc failed", 0); }

       wchar_t *buf = static_cast<wchar_t *>(GlobalLock(hMem));
       if (!buf) {
        GlobalFree(hMem);
        throw CalcError(CalcErrorType::InternalError, "GlobalLock failed", 0);
       }

       MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, buf, size);
       GlobalUnlock(hMem);

       if (!SetClipboardData(CF_UNICODETEXT, hMem)) {
        GlobalFree(hMem);
        throw CalcError(CalcErrorType::InternalError, "SetClipboardData failed", 0);
       }

      } catch (...) {
       CloseClipboard();
       throw;
      }

      CloseClipboard();

      if (!result.explain.empty() && result.explain.back() != '\n') { result.explain += '\n'; }
      result.explain += "clipboard pasted (" + whatByteUnit(s.size()) + ")\n";
#else
      throw CalcError(CalcErrorType::InternalError, "clipboard is not supported on this platform", 0);
#endif
      break;
     }
   }
  }

  ectx.sideEffects.clear();
 }

} // namespace mm::cal