#include "formatter.hpp"
#include <algorithm>
#include <array>
#include <charconv>
#include <string>

namespace mm::cal {

 inline static void trimZeros(std::string &s) {
  if (auto pos = s.find('.'); pos != std::string::npos) {
   while (!s.empty() && s.back() == '0')
    s.pop_back();
   if (!s.empty() && s.back() == '.') s.pop_back();
  }
 }

 std::string formatReal(double x, const SystemConfig &cfg) {
  std::ostringstream oss;

  if (std::isnan(x)) return "nan";
  if (std::isinf(x)) return (x > 0 ? "inf" : "-inf");

  x = normalizeZero(x); // ゼロ補正

  oss << std::fixed << std::setprecision(cnst_precision) << x; // 指数表記禁止(暫定)

  std::string s = oss.str();

  if (s.find('.') != std::string::npos) { // 末尾のゼロ消し
   while (!s.empty() && s.back() == '0')
    s.pop_back();
   if (!s.empty() && s.back() == '.') s.pop_back();
  }

  if (s.empty()) s = "0";

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

 inline void appendReal(double x, const SystemConfig &cfg, std::string &out) { out += formatReal(x, cfg); }

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
        out += "[internal-error]";
       }
      },
      v.storage());
 }

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
        throw CalcError(CalcErrorType::InternalError, "InternalError: unknown Value type", 0);
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

} // namespace mm::cal