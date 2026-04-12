#include "ast.hpp"

namespace mm::cal {
 UserFunction::UserFunction() = default;
 UserFunction::~UserFunction() = default;
 UserFunction::UserFunction(UserFunction &&) noexcept = default;
 UserFunction &UserFunction::operator=(UserFunction &&) noexcept = default;

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
 AssignNode::AssignNode(std::string n, std::unique_ptr<ASTNode> e, size_t p) : name(std::move(n)), expr(std::move(e)) { pos = p; }
 Value AssignNode::evalImpl(EvaluationContext &ctx) const {
  if (ctx.session.userFunctions.contains(name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: name already used as function", pos);
  Value v = expr->eval(ctx);
  if (ctx.callStack.empty()) {
   ctx.session.globals[name] = v;
  } else {
   ctx.locals[name] = v;
  }
  return v;
 }

 FunctionDefNode::FunctionDefNode(std::string n, std::vector<std::string> s, std::unique_ptr<ASTNode> b, size_t p) : name(std::move(n)), params(std::move(s)), body(std::move(b)) { pos = p; }
 Value FunctionDefNode::evalImpl(EvaluationContext &ctx) const {
  if (ctx.locals.contains(name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: name already used as variable", pos);           // 変数との衝突
  if (ctx.session.globals.contains(name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: name already used as variable", pos);  // 変数との衝突
  if (ctx.session.cfg.functions.contains(name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: cannot override builtin", pos);  // builtinとの衝突
  if (ctx.session.userFunctions.contains(name)) throw CalcError(CalcErrorType::DefinitionError, "DefinitionError: function already defined", pos); // 再定義禁止

  { // 仮引数重複チェック
   std::unordered_set<std::string> seen;
   for (const auto &p : params)
    if (!seen.insert(p).second) throw CalcError(CalcErrorType::RuntimeError, "RuntimeError: duplicate parameter name", pos);
  }

  UserFunction fn;
  fn.params = params;

  // bodyは所有権を持たせる(現在の設計ではASTは一回評価用なのでmoveできない。)
  // clone関数がないならdeep copyが必要。 ここでは簡易的に const_cast + clone前提で記述：

  fn.body = body->clone(); // clone() をASTNodeに用意すること

  ctx.session.userFunctions[name] = std::move(fn);

  return Value(); // none
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
  if (std::find(ectx.callStack.begin(), ectx.callStack.end(), name) != ectx.callStack.end()) { throw CalcError(CalcErrorType::DefinitionError, "recursion disabled", pos); } // 自己再帰検出

  // user function 優先
  if (auto uit = ectx.session.userFunctions.find(name); uit != ectx.session.userFunctions.end()) {
   const UserFunction &fn = uit->second;
   if (args.size() != fn.params.size()) throw CalcError(CalcErrorType::FunctionMissing, errorMessage(CalcErrorType::FunctionMissing), pos);
   std::vector<Value> values;
   values.reserve(args.size());
   for (auto &a : args)
    values.push_back(a->eval(ectx));

   // ローカルコンテキスト生成
   EvaluationContext local(ectx.session);
   local.callStack = ectx.callStack;
   local.callStack.push_back(name);

   for (size_t i = 0; i < fn.params.size(); ++i) {
    local.locals[fn.params[i]] = values[i];
   }

   return fn.body->eval(local);
  }

  // 内蔵関数
  auto it = ectx.session.cfg.functions.find(name);
  if (it == ectx.session.cfg.functions.end()) throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), pos);

  auto &f = it->second;
  int argc = static_cast<int>(args.size());

  if (!f.validArgc(argc)) throw CalcError(CalcErrorType::FunctionMissing, errorMessage(CalcErrorType::FunctionMissing), pos);

  std::vector<Value> v;
  v.reserve(args.size());
  for (auto &a : args)
   v.push_back(a->eval(ectx));

  FunctionContext ctx{ectx.session, pos, ectx.sideEffects};

  return f.f(v, ctx);
 }

 SymbolNode::SymbolNode(std::string n, size_t p) : name(std::move(n)) { pos = p; }

 Value SymbolNode::evalImpl(EvaluationContext &ectx) const {
  if (ectx.locals.contains(name)) return ectx.locals.at(name);
  if (ectx.session.globals.contains(name)) return ectx.session.globals.at(name);
  if (constants.count(name)) return constants.at(name);
  if (symbols.count(name)) return name;
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
  else i = static_cast<int>(ectx.session.history.size()) + index + 1;

  if (i <= 0 || i > static_cast<int>(ectx.session.history.size())) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos);

  return ectx.session.history[i - 1].value;
 }
 Value evalCompare(const Value &lhs, const Value &rhs, CmpOp op, FunctionContext &ctx) {
  // 複素数は大小比較不可
  if (lhs.isComplex() || rhs.isComplex()) throw CalcError(CalcErrorType::NotImplemented, "complex comparison is not implemented", ctx.pos);
  double a = lhs.asScalar(ctx.pos), b = rhs.asScalar(ctx.pos);
  bool r = false;
  int prec = ctx.session.cfg.precision;

  switch (op) {
   case CmpOp::Equal: r = nearly_equal(a, b); break;
   case CmpOp::NotEqual: r = !nearly_equal(a, b); break;
   case CmpOp::Less: r = a < b && !nearly_equal(a, b); break;
   case CmpOp::LessEq: r = a < b || nearly_equal(a, b); break;
   case CmpOp::Greater: r = a > b && !nearly_equal(a, b); break;
   case CmpOp::GreaterEq: r = a > b || nearly_equal(a, b); break;
  }
  return r ? 1.0 : 0.0;
 }
} // namespace mm::cal