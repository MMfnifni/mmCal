#include "ast.hpp"

namespace mm::cal {
 /* ============================
    AST ヘルパくん
    ============================ */

 CmpOp tokenToCmpOp(TokenType t) {
  switch (t) {
   case TokenType::Less: return CmpOp::Less;
   case TokenType::LessEq: return CmpOp::LessEq;
   case TokenType::Greater: return CmpOp::Greater;
   case TokenType::GreaterEq: return CmpOp::GreaterEq;
   case TokenType::Equal: return CmpOp::Equal;
   case TokenType::NotEqual: return CmpOp::NotEqual;
   default: throw std::logic_error("not a compare token");
  }
 }
 BinOp tokenToBinOp(TokenType t) {
  switch (t) {
   case TokenType::Plus: return BinOp::Add;
   case TokenType::Minus: return BinOp::Sub;
   case TokenType::Mul: return BinOp::Mul;
   case TokenType::Div: return BinOp::Div;
   case TokenType::Pow: return BinOp::Pow;
   default: throw CalcError(CalcErrorType::InvalidOperation, errorMessage(CalcErrorType::InvalidOperation), 0);
  }
 }

 /* ============================
    AST
    ============================ */

 // AST ノード系
 NumberNode::NumberNode(Value v, size_t p) : value(std::move(v)) { pos = p; }
 Value NumberNode::evalImpl(EvaluationContext &ectx) const { return value; }
 AssignNode::AssignNode(std::string n, std::unique_ptr<ASTNode> e) : name(std::move(n)), expr(std::move(e)) {}
 Value AssignNode::evalImpl(EvaluationContext &ctx) const {
  Value v = expr->eval(ctx);
  ctx.variables[name] = v;
  return v;
 }
 BinaryNode::BinaryNode(BinOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }
 Value BinaryNode::evalImpl(EvaluationContext &ectx) const {
  auto a = lhs->eval(ectx);
  auto b = rhs->eval(ectx);
  if (a.isMulti() || b.isMulti()) throw CalcError(CalcErrorType::TypeError, "operator not defined for multivalue", pos);
  switch (op) {
   case BinOp::Add: return add(a, b, pos);
   case BinOp::Sub: return sub(a, b, pos);
   case BinOp::Mul: return mul(a, b, pos);
   case BinOp::Div: return div(a, b, pos);
   case BinOp::Pow: return power(a, b, pos);
  }

  throw CalcError(CalcErrorType::InvalidOperation, "invalid op", pos);
 }

 UnaryNode::UnaryNode(UnaryOp o, std::unique_ptr<ASTNode> r, size_t p) : op(o), rhs(std::move(r)) { pos = p; }
 Value UnaryNode::evalImpl(EvaluationContext &ectx) const {
  auto v = rhs->eval(ectx);
  if (v.isMulti()) throw CalcError(CalcErrorType::TypeError, "unary operator not defined for multivalue", pos);
  if (op == UnaryOp::Minus) return negate(v, pos);
  return v;
 }

 PostfixUnaryNode::PostfixUnaryNode(char op, std::unique_ptr<ASTNode> e, size_t p) : op(op), expr(std::move(e)) { pos = p; }

 Value PostfixUnaryNode::evalImpl(EvaluationContext &ectx) const {
  Value v = expr->eval(ectx);

  if (op == '!') return mm::cal::factorial(v, pos);

  throw CalcError(CalcErrorType::InvalidOperation, "invalid postfix", pos);
 }

 CompareNode::CompareNode(CmpOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }
 Value CompareNode::evalImpl(EvaluationContext &ectx) const {
  Value a = lhs->eval(ectx);
  Value b = rhs->eval(ectx);

  if (a.isMulti() || b.isMulti()) throw CalcError(CalcErrorType::TypeError, "comparison not defined for multivalue", pos);
  mm::cal::CompareOp cop;

  switch (op) {
   case CmpOp::Equal: cop = CompareOp::Eq; break;
   case CmpOp::NotEqual: cop = CompareOp::Ne; break;
   case CmpOp::Less: cop = CompareOp::Lt; break;
   case CmpOp::LessEq: cop = CompareOp::Le; break;
   case CmpOp::Greater: cop = CompareOp::Gt; break;
   case CmpOp::GreaterEq: cop = CompareOp::Ge; break;
   default: throw CalcError(CalcErrorType::InvalidOperation, "invalid operator", pos);
  }

  return mm::cal::compare(a, b, cop, pos);
 }

 Value FunctionCallNode::evalImpl(EvaluationContext &ectx) const {

  auto it = ectx.cfg.functions.find(name);
  if (it == ectx.cfg.functions.end()) throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), pos);

  auto &f = it->second;
  int argc = static_cast<int>(args.size());
  if (!f.validArgc(argc)) throw CalcError(CalcErrorType::FunctionMissing, errorMessage(CalcErrorType::FunctionMissing), pos);
  std::vector<Value> v;
  v.reserve(args.size());
  for (auto &a : args)
   v.push_back(a->eval(ectx));

  FunctionContext ctx{ectx.cfg, ectx.hist, ectx.base, pos, ectx.sideEffects};
  Value r = f.f(v, ctx);
  // double のときだけ有限チェック -> 一括捕捉に変更
  /* if (std::holds_alternative<double>(r)) checkFinite(std::get<double>(r), pos);*/
  return r;
 }

 SymbolNode::SymbolNode(std::string n, size_t p) : name(std::move(n)) { pos = p; }

 Value SymbolNode::evalImpl(EvaluationContext &ectx) const {
  if (ectx.variables.contains(name)) return ectx.variables.at(name); // ユーザ変数(パーサは置換するものじゃない)
  if (constants.count(name)) return constants.at(name);              // 定数はここで解決できる（parse段階で解決しなくてよくなる）
  if (symbols.count(name)) return name;                              // symbols(deg, rad, mm...) はオプション指定子としては文字列扱いで返す
  // return name; // それ以外も「オプション指定子」として許可するならここで返す(要検討!!!)

  throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), pos);
 }

 UnitApplyNode::UnitApplyNode(std::unique_ptr<ASTNode> e, std::string u, size_t p) : expr(std::move(e)), unit(std::move(u)) { pos = p; }

 Value UnitApplyNode::evalImpl(EvaluationContext &ectx) const {
  Value v = expr->eval(ectx);
  double x = v.asScalar(pos);

  // angle
  if (unit == "deg") return x;
  if (unit == "rad") { return x * 180.0 / PI; }
  if (unit == "grad") { return x * 0.9; }

  // length (将来)
  if (unit == "mm" || unit == "cm" || unit == "m" || unit == "inch") { throw CalcError(CalcErrorType::NotImplemented, "length units are not implemented yet", pos); }

  throw CalcError(CalcErrorType::UnknownIdentifier, "unknown unit: " + unit, pos);
 }

 MultiLiteralNode::MultiLiteralNode(std::vector<std::unique_ptr<ASTNode>> e, size_t p) {
  pos = p;
  elems = std::move(e);
 }

 Value MultiLiteralNode::evalImpl(EvaluationContext &ectx) const {
  std::vector<Value> values;
  values.reserve(elems.size());

  for (const auto &n : elems) {
   Value v = n->eval(ectx);
   values.push_back(std::move(v));
  }

  return Value(std::make_shared<MultiValue>(std::move(values)));
 }

 OutNode::OutNode(int idx, size_t p) : index(idx) { pos = p; }

 Value OutNode::evalImpl(EvaluationContext &ectx) const {
  if (index == 0) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos);

  int i;
  if (index > 0) i = index;
  else i = static_cast<int>(ectx.hist.size()) + index + 1;

  if (i <= 0 || i > static_cast<int>(ectx.hist.size())) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos);

  return ectx.hist[i - 1].value;
 }
 Value evalCompare(const Value &lhs, const Value &rhs, CmpOp op, FunctionContext &ctx) {
  // 複素数は大小比較不可
  if (lhs.isComplex() || rhs.isComplex()) throw CalcError(CalcErrorType::NotImplemented, "complex comparison is not implemented", ctx.pos);
  double a = lhs.asScalar(ctx.pos), b = rhs.asScalar(ctx.pos);
  bool r = false;
  int prec = ctx.cfg.precision;

  switch (op) {
   case CmpOp::Equal: r = tolerantEqual(a, b, prec); break;
   case CmpOp::NotEqual: r = !tolerantEqual(a, b, prec); break;
   case CmpOp::Less: r = a < b && !tolerantEqual(a, b, prec); break;
   case CmpOp::LessEq: r = a < b || tolerantEqual(a, b, prec); break;
   case CmpOp::Greater: r = a > b && !tolerantEqual(a, b, prec); break;
   case CmpOp::GreaterEq: r = a > b || tolerantEqual(a, b, prec); break;
  }
  return r ? 1.0 : 0.0;
 }
} // namespace mm::cal