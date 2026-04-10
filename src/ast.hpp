#pragma once

#include "core.hpp"
#include "math_util.hpp"
#include "value_arithmetic.hpp"
#include "value_compare.hpp"
#include <cctype>
#include <charconv>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

namespace mm::cal {

 inline const std::unordered_set<std::string> symbols = {"deg", "rad", "grad", "mm", "cm", "m", "inch"};

 enum class TokenType {
  End,
  Number,
  Identifier,
  Plus,  // +
  Minus, // -
  Mul,   // *
  Div,   // /
  Pow,   // ^
  // Star,
  // Slash,
  // Caret,
  LParen,    // (
  RParen,    // )
  LBracket,  // [
  RBracket,  // ]
  LBrace,    // {
  RBrace,    // }
  Comma,     // ,
  Bang,      // !
  Percent,   // %
  Less,      // <
  LessEq,    // <=
  Greater,   // >
  GreaterEq, // >=
  Equal,     // ==
  NotEqual,  // !=
  Assign,    // =
  String,    // ""
 };

 /* ============================
    AST
    ============================ */

 enum class UnaryOp { Plus, Minus };
 enum class BinOp { Add, Sub, Mul, Div, Pow };
 enum class CmpOp { Less, LessEq, Greater, GreaterEq, Equal, NotEqual };

 CmpOp tokenToCmpOp(TokenType t);
 BinOp tokenToBinOp(TokenType t);

 struct ASTNode {
   size_t pos = 0;
   virtual ~ASTNode() = default;
   // virtual Value eval(EvaluationContext &) const = 0;

   // 外から呼ばれる唯一の eval
   virtual Value eval(EvaluationContext &ectx) const final {
    Value v = evalImpl(ectx);
    Value::checkFinite(v, pos);
    return v;
   }
   // virtual std::unique_ptr<ASTNode> clone() const = 0;

  protected:
   // 派生クラスはこれだけ実装する
   virtual Value evalImpl(EvaluationContext &) const = 0;

  public:
   std::unique_ptr<ASTNode> clone() const { return cloneImpl(); }

  private:
   virtual std::unique_ptr<ASTNode> cloneImpl() const = 0;
 };

 template <class Derived> struct ASTNodeCRTP : ASTNode {
  private:
   std::unique_ptr<ASTNode> cloneImpl() const override { return std::make_unique<Derived>(static_cast<const Derived &>(*this)); }
 };

 /* ---------- Number ---------- */
 struct NumberNode : ASTNodeCRTP<NumberNode> {
   Value value;

   NumberNode(Value v, size_t p);
   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ------ User variables and func ------- */

 struct AssignNode : ASTNodeCRTP<AssignNode> {
   std::string name;
   std::unique_ptr<ASTNode> expr;

   AssignNode(std::string n, std::unique_ptr<ASTNode> e, size_t p);
   AssignNode(const AssignNode &other) : name(other.name), expr(other.expr ? other.expr->clone() : nullptr) { pos = other.pos; }
   Value evalImpl(EvaluationContext &ectx) const override;
 };

 struct FunctionDefNode : ASTNodeCRTP<FunctionDefNode> {
   std::string name;
   std::vector<std::string> params;
   std::unique_ptr<ASTNode> body;

   FunctionDefNode(std::string n, std::vector<std::string> s, std::unique_ptr<ASTNode> b, size_t p);
   FunctionDefNode(const FunctionDefNode &other) : name(other.name), params(other.params), body(other.body ? other.body->clone() : nullptr) { pos = other.pos; }

   Value evalImpl(EvaluationContext &ctx) const override;
 };

 /* ---------- Binary ---------- */

 struct BinaryNode : ASTNodeCRTP<BinaryNode> {
   BinOp op;
   std::unique_ptr<ASTNode> lhs, rhs;

   BinaryNode(BinOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p);
   BinaryNode(const BinaryNode &other) : op(other.op), lhs(other.lhs ? other.lhs->clone() : nullptr), rhs(other.rhs ? other.rhs->clone() : nullptr) { pos = other.pos; }

   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ---------- Unary ---------- */
 struct UnaryNode : ASTNodeCRTP<UnaryNode> {
   UnaryOp op;
   std::unique_ptr<ASTNode> rhs;
   UnaryNode(UnaryOp o, std::unique_ptr<ASTNode> r, size_t p);
   UnaryNode(const UnaryNode &other) : op(other.op), rhs(other.rhs ? other.rhs->clone() : nullptr) { pos = other.pos; }

   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ---------- Postfix ---------- */
 struct PostfixUnaryNode : ASTNodeCRTP<PostfixUnaryNode> {
   char op;
   std::unique_ptr<ASTNode> expr;

   PostfixUnaryNode(char op, std::unique_ptr<ASTNode> e, size_t p);
   PostfixUnaryNode(const PostfixUnaryNode &other) : op(other.op), expr(other.expr ? other.expr->clone() : nullptr) { pos = other.pos; }

   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ---------- Compare ---------- */
 struct CompareNode : ASTNodeCRTP<CompareNode> {
   CmpOp op;
   std::unique_ptr<ASTNode> lhs, rhs;

   CompareNode(CmpOp o, std::unique_ptr<ASTNode> l, std::unique_ptr<ASTNode> r, size_t p);
   CompareNode(const CompareNode &other) : op(other.op), lhs(other.lhs ? other.lhs->clone() : nullptr), rhs(other.rhs ? other.rhs->clone() : nullptr) { pos = other.pos; }

   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ---------- Function ---------- */

 struct FunctionCallNode : ASTNodeCRTP<FunctionCallNode> {
   std::string name;
   std::vector<std::unique_ptr<ASTNode>> args;
   std::unordered_map<std::string, std::unique_ptr<ASTNode>> options;

   FunctionCallNode() = default;

   FunctionCallNode(const FunctionCallNode &other) : name(other.name) {
    pos = other.pos;
    args.reserve(other.args.size());
    for (const auto &a : other.args)
     args.push_back(a ? a->clone() : nullptr);

    for (const auto &[k, v] : other.options)
     options.emplace(k, v ? v->clone() : nullptr);
   }

   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ---------- SymbolNode ---------- */

 struct SymbolNode : ASTNodeCRTP<SymbolNode> {
   std::string name;
   SymbolNode(std::string n, size_t p);

   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ---------- Unit ---------- */

 struct UnitApplyNode : ASTNodeCRTP<UnitApplyNode> {
   std::unique_ptr<ASTNode> expr;
   std::string unit;

   UnitApplyNode(std::unique_ptr<ASTNode> e, std::string u, size_t p);
   UnitApplyNode(const UnitApplyNode &other) : expr(other.expr ? other.expr->clone() : nullptr), unit(other.unit) { pos = other.pos; }
   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ---------- MultiLiteral ---------- */
 struct MultiLiteralNode : ASTNodeCRTP<MultiLiteralNode> {
   std::vector<std::unique_ptr<ASTNode>> elems;

   MultiLiteralNode(std::vector<std::unique_ptr<ASTNode>> e, size_t p);
   MultiLiteralNode(const MultiLiteralNode &other) {
    pos = other.pos;
    elems.reserve(other.elems.size());
    for (const auto &e : other.elems)
     elems.push_back(e ? e->clone() : nullptr);
   }

  protected:
   Value evalImpl(EvaluationContext &ectx) const override;
 };

 /* ---------- History Out[n] ---------- */
 struct OutNode : ASTNodeCRTP<OutNode> {
   int index; // 正: Out[n], 負: % / %%
   OutNode(int idx, size_t p);

   Value evalImpl(EvaluationContext &ectx) const override;
 };
 Value evalCompare(const Value &lhs, const Value &rhs, CmpOp op, FunctionContext &ctx);

 // ASTヘルパ・その他
 inline std::unique_ptr<ASTNode> cloneAST(const ASTNode *n) {
  if (!n) return nullptr;
  return n->clone();
 }
 struct UserFunction;

} // namespace mm::cal