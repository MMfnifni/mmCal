#include "explain_formatter.hpp"
#include "formatter.hpp"
#include <functional>
#include <sstream>

namespace mm::cal {
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
   double x = normalizeZero(value.asScalar(0));

   oss << "Canonical form = " << formatReal(x, cfg) << "\n";
   oss << "Absolute value = " << formatReal(std::abs(x), cfg) << "\n";
   oss << "Zero-normalization threshold = " << formatReal(cnst_precision_inv, cfg) << "\n";

   if (x == 0.0) {
    oss << "Interpretation = Treated as zero\n";
   } else if (std::abs(x) < cnst_precision_inv) {
    oss << "Interpretation = Very small real value\n";
   }

  } else if (value.isComplex()) {
   std::complex<double> z = value.asComplex(0);
   double re = normalizeZero(z.real());
   double im = normalizeZero(z.imag());
   std::complex<double> normalized(re, im);

   oss << "Canonical form = " << formatComplex(normalized, cfg) << "\n";
   oss << "Real part = " << formatReal(re, cfg) << "\n";
   oss << "Imaginary part = " << formatReal(im, cfg) << "\n";
   oss << "Absolute value = " << formatReal(std::abs(normalized), cfg) << "\n";
   oss << "Zero-normalization threshold = " << formatReal(cnst_precision_inv, cfg) << "\n";

   if (im == 0.0) {
    oss << "Interpretation = Purely real complex value\n";
   } else if (re == 0.0) {
    oss << "Interpretation = Purely imaginary value\n";
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

   oss << "Top-level size = " << mv->size() << "\n";
   oss << "Maximum depth = " << dataExplainMaxDepth(value) << "\n";
   oss << "Homogeneous = " << (summary.homogeneous ? "true" : "false") << "\n";
   oss << dataExplainContainsLine(summary) << "\n";

  } else {
   oss << "Interpretation = No value\n";
  }

  oss << "Shape = " << shape << "\n";
  oss << "Legend = R:Real, C:Complex, S:String, N:None\n";

  oss << "\n[Data Explain END]\n\n";
  return oss.str();
 }
} // namespace mm::cal