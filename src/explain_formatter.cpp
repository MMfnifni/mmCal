#include "explain_formatter.hpp"
#include "formatter.hpp"
#include <bit>
#include <cstdint>
#include <functional>
#include <iomanip>
#include <sstream>

namespace mm::cal {

 inline static std::string dataExplainDoubleBits(double x) {
  uint64_t bits = std::bit_cast<uint64_t>(x);
  std::ostringstream oss;
  oss << "0x" << std::hex << std::nouppercase << std::setw(16) << std::setfill('0') << bits;
  return oss.str();
 }

 struct DataExplainMatrixInfo {
   bool matrixLike = false;  // top-level が rows で、各要素が MultiValue
   bool rectangular = false; // すべての row 長が等しい
   bool square = false;      // rectangular かつ rows == cols
   size_t rows = 0;
   size_t cols = 0;
 };

 inline static DataExplainMatrixInfo dataExplainMatrixShape(const Value &value) {
  DataExplainMatrixInfo info;

  if (!value.isMulti()) return info;

  auto top = value.asMulti(0);
  info.rows = top->size();

  // 空は matrix-like とまでは言わないですよ
  if (info.rows == 0) return info;

  bool allRowsAreMulti = true;
  size_t cols = 0;
  bool first = true;
  bool rectangular = true;

  for (const auto &elem : top->elems_) {
   if (!elem.isMulti()) {
    allRowsAreMulti = false;
    break;
   }

   auto row = elem.asMulti(0);
   if (first) {
    cols = row->size();
    first = false;
   } else if (row->size() != cols) {
    rectangular = false;
   }
  }

  if (!allRowsAreMulti) return info;

  info.matrixLike = true;
  info.rectangular = rectangular;
  info.cols = cols;
  info.square = rectangular && (info.rows == info.cols);
  return info;
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
   return "N";
  };

  return dump(value);
 }

 size_t dataExplainMaxDepth(const Value &value) {
  std::function<size_t(const Value &, size_t)> walk = [&](const Value &v, size_t depth) -> size_t {
   if (!v.isMulti()) return depth;

   auto mv = v.asMulti(0);
   size_t best = depth;
   for (const auto &elem : mv->elems_) {
    best = std::max(best, walk(elem, depth + 1));
   }
   return best;
  };

  return walk(value, 0);
 }

 struct DataExplainSummary {
   bool hasReal = false;
   bool hasComplex = false;
   bool hasString = false;
   bool hasNone = false;
   bool homogeneous = true;
   bool firstSeen = false;
   char firstKind = '\0';
 };

 void dataExplainCollectSummary(const Value &value, DataExplainSummary &summary) {
  std::function<void(const Value &)> walk = [&](const Value &v) {
   char kind = 'N';

   if (v.isMulti()) {
    for (const auto &elem : v.asMulti(0)->elems_) {
     walk(elem);
    }
    return;
   }

   if (v.isScalar()) {
    kind = 'R';
    summary.hasReal = true;
   } else if (v.isComplex()) {
    kind = 'C';
    summary.hasComplex = true;
   } else if (v.isString()) {
    kind = 'S';
    summary.hasString = true;
   } else {
    kind = 'N';
    summary.hasNone = true;
   }

   if (!summary.firstSeen) {
    summary.firstSeen = true;
    summary.firstKind = kind;
   } else if (summary.firstKind != kind) {
    summary.homogeneous = false;
   }
  };

  walk(value);
 }

 std::string dataExplainContainsLine(const DataExplainSummary &summary) {
  std::vector<std::string> parts;

  if (summary.hasReal) parts.push_back("Real");
  if (summary.hasComplex) parts.push_back("Complex");
  if (summary.hasString) parts.push_back("String");
  if (summary.hasNone) parts.push_back("None");

  if (parts.empty()) return "Contains = None";

  std::ostringstream oss;
  oss << "Contains = ";
  for (size_t i = 0; i < parts.size(); ++i) {
   if (i > 0) oss << ", ";
   oss << parts[i];
  }
  return oss.str();
 }

 std::string dataExplainQuotedString(const std::string &s) {
  std::ostringstream oss;
  oss << '"';

  for (char ch : s) {
   switch (ch) {
    case '\\': oss << "\\\\"; break;
    case '"': oss << "\\\""; break;
    case '\n': oss << "\\n"; break;
    case '\r': oss << "\\r"; break;
    case '\t': oss << "\\t"; break;
    default: oss << ch; break;
   }
  }

  oss << '"';
  return oss.str();
 }

 std::string dataExplain(Value value, const SystemConfig &cfg) {
  std::ostringstream oss;

  oss << "\n[Data Explain]\n\n";

  const std::string displayed = formatResult(value, cfg);
  const std::string shape = dataExplainStructure(value);

  oss << "Type = ";
  if (value.isScalar()) {
   oss << "Real Number\n";
  } else if (value.isComplex()) {
   oss << "Complex Number\n";
  } else if (value.isString()) {
   oss << "String\n";
  } else if (value.isMulti()) {
   oss << "MultiValue\n";
  } else {
   oss << "None\n";
  }

  oss << "Displayed = " << displayed << "\n";

  if (value.isScalar()) {
   double raw = value.asScalar(0);
   double x = normalizeZero(raw);

   oss << "Canonical form = " << formatReal(x, cfg) << "\n";

   std::ostringstream ro, rabs;
   ro << std::setprecision(17) << std::scientific << raw;
   rabs << std::setprecision(17) << std::scientific << std::abs(raw);

   oss << "Raw value = " << ro.str() << "\n";
   oss << "Raw value (bits) = " << dataExplainDoubleBits(raw) << "\n";
   oss << "Raw absolute value = " << rabs.str() << "\n";
   oss << "Zero-normalized absolute value = " << formatReal(std::abs(x), cfg) << "\n";
   oss << "Zero-normalization threshold = " << formatReal(cnst_precision_inv, cfg) << "\n";

   if (raw != 0.0 && x == 0.0) {
    oss << "Interpretation = Nonzero value hidden by zero-normalization\n";
   } else if (raw == 0.0) {
    oss << "Interpretation = Exactly zero\n";
   } else if (std::abs(raw) < cnst_precision_inv) {
    oss << "Interpretation = Very small real value\n";
   }

  } else if (value.isComplex()) {
   std::complex<double> z = value.asComplex(0);
   double rawRe = z.real();
   double rawIm = z.imag();
   double re = normalizeZero(rawRe);
   double im = normalizeZero(rawIm);

   std::ostringstream rre, rim, rabs;
   rre << std::setprecision(17) << std::scientific << rawRe;
   rim << std::setprecision(17) << std::scientific << rawIm;
   rabs << std::setprecision(17) << std::scientific << std::abs(z);

   std::complex<double> normalized(re, im);

   oss << "Canonical form = " << formatComplex(normalized, cfg) << "\n";
   oss << "Real part = " << formatReal(re, cfg) << "\n";
   oss << "Imaginary part = " << formatReal(im, cfg) << "\n";
   oss << "Raw real part = " << rre.str() << "\n";
   oss << "Raw real part (bits) = " << dataExplainDoubleBits(rawRe) << "\n";
   oss << "Raw imaginary part = " << rim.str() << "\n";
   oss << "Raw imaginary part (bits) = " << dataExplainDoubleBits(rawIm) << "\n";
   oss << "Raw absolute value = " << rabs.str() << "\n";
   oss << "Zero-normalized absolute value = " << formatReal(std::abs(normalized), cfg) << "\n";
   oss << "Zero-normalization threshold = " << formatReal(cnst_precision_inv, cfg) << "\n";

   bool hiddenRe = (rawRe != 0.0 && re == 0.0);
   bool hiddenIm = (rawIm != 0.0 && im == 0.0);

   if (hiddenRe || hiddenIm) {
    if (hiddenRe && hiddenIm) {
     oss << "Interpretation = Nonzero real and imaginary components hidden by zero-normalization\n";
    } else if (hiddenRe) {
     oss << "Interpretation = Nonzero real component hidden by zero-normalization\n";
    } else {
     oss << "Interpretation = Nonzero imaginary component hidden by zero-normalization\n";
    }
   } else if (rawRe == 0.0 && rawIm == 0.0) {
    oss << "Interpretation = Exactly zero complex value\n";
   } else if (im == 0.0) {
    oss << "Interpretation = Purely real after zero-normalization\n";
   } else if (re == 0.0) {
    oss << "Interpretation = Purely imaginary after zero-normalization\n";
   }

  } else if (value.isString()) {
   const std::string &s = value.asString(0);

   oss << "Length = " << s.size() << "\n";
   oss << "Quoted = " << dataExplainQuotedString(s) << "\n";

   if (s.empty()) { oss << "Interpretation = Empty string\n"; }

  } else if (value.isMulti()) {
   auto mv = value.asMulti(0);

   DataExplainSummary summary;
   dataExplainCollectSummary(value, summary);

   DataExplainMatrixInfo matrixInfo = dataExplainMatrixShape(value);

   oss << "Top-level size = " << mv->size() << "\n";
   oss << "Maximum depth = " << dataExplainMaxDepth(value) << "\n";
   oss << "Homogeneous = " << (summary.homogeneous ? "true" : "false") << "\n";
   oss << dataExplainContainsLine(summary) << "\n";

   oss << "Matrix-like = " << (matrixInfo.matrixLike ? "true" : "false") << "\n";
   if (matrixInfo.matrixLike) {
    oss << "Rectangular = " << (matrixInfo.rectangular ? "true" : "false") << "\n";
    oss << "Rows = " << matrixInfo.rows << "\n";
    oss << "Cols = " << matrixInfo.cols << "\n";
    oss << "Square = " << (matrixInfo.square ? "true" : "false") << "\n";
   }
  } else {
   oss << "Interpretation = No value\n";
  }

  oss << "Shape = " << shape << "\n";
  oss << "Legend = R:Real, C:Complex, S:String, N:None\n";

  oss << "\n[Data Explain END]\n\n";
  return oss.str();
 }
} // namespace mm::cal