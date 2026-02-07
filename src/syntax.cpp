#include "syntax.hpp"
#include "core.hpp"

namespace mm::cal {
 /* ============================
    Token
============================ */

 /* ============================
    Lexer
============================ */

 const Token &Lexer::peek() {
  if (!hasPeek) {
   lookahead = nextToken();
   hasPeek = true;
  }
  return lookahead;
 }

 Token Lexer::get() {
  if (hasPeek) {
   hasPeek = false;
   return lookahead;
  }
  return nextToken();
 }

 bool Lexer::eof() { return peek().type == TokenType::End; }

 void Lexer::skipSpaces() {
  while (p < src.size() && std::isspace((unsigned char)src[p]))
   ++p;
 }

 Token Lexer::nextToken() {
  skipSpaces();
  size_t start = p;

  if (p >= src.size()) return {TokenType::End, "", start};

  char c = src[p++];

  // number
  if (std::isdigit(c) || c == '.') {
   --p;

   const char *begin = src.data() + p;
   const char *end = src.data() + src.size();

   double value;
   auto [ptr, ec] = std::from_chars(begin, end, value);

   if (ec == std::errc::invalid_argument) throw CalcError(CalcErrorType::InvalidNumber, "invalid number", p);

   if (ec == std::errc::result_out_of_range) throw CalcError(CalcErrorType::Overflow, "number out of range", p);

   // 1..1, 1.2.3潰し
   if (ptr < end && *ptr == '.') { throw CalcError(CalcErrorType::InvalidNumber, "multiple '.' in number", static_cast<size_t>(ptr - src.data())); }

   size_t len = static_cast<size_t>(ptr - begin);
   std::string s(src.substr(p, len));

   p += len;
   return {TokenType::Number, s, value, start};
  }

  // identifier
  if (std::isalpha(c) || c == '_') {
   std::string s(1, c);
   while (p < src.size() && (std::isalnum(src[p]) || src[p] == '_'))
    s += src[p++];
   return {TokenType::Identifier, s, start};
  }

  if (c == '<') {
   if (p < src.size() && src[p] == '=') {
    ++p;
    return {TokenType::LessEq, "<=", start};
   }
   return {TokenType::Less, "<", start};
  }

  if (c == '>') {
   if (p < src.size() && src[p] == '=') {
    ++p;
    return {TokenType::GreaterEq, ">=", start};
   }
   return {TokenType::Greater, ">", start};
  }

  if (c == '=') {
   if (p < src.size() && src[p] == '=') {
    ++p;
    return {TokenType::Equal, "==", start};
   }
   throw CalcError(CalcErrorType::InvalidCharacter, "single '=' not allowed", start);
  }

  switch (c) {
   case '+': return {TokenType::Plus, "+", start};
   case '-': return {TokenType::Minus, "-", start};
   case '*': return {TokenType::Mul, "*", start};
   case '/': return {TokenType::Div, "/", start};
   case '^': return {TokenType::Pow, "^", start};
   case '(': return {TokenType::LParen, "(", start};
   case ')': return {TokenType::RParen, ")", start};
   case '[': return {TokenType::LBracket, "[", start};
   case ']': return {TokenType::RBracket, "]", start};
   case ',': return {TokenType::Comma, ",", start};
   case '!': return {TokenType::Bang, "!", start};
   case '%': return {TokenType::Percent, "%", start};
  }

  throw CalcError(CalcErrorType::InvalidCharacter, std::string("unexpected character: ") + c, start);
 }

 /* ============================
   Parser（暗黙乗算）
   ============================ */

 Parser::Parser(SystemConfig &cfg, const std::string &s) : cfg(cfg), lex(s) { cur = lex.get(); }

 void Parser::advance() {
  prev = cur;
  cur = lex.get();
 }

 bool Parser::accept(TokenType t) {
  if (cur.type == t) {
   advance();
   return true;
  }
  return false;
 }

 void Parser::expect(TokenType t) {
  if (!accept(t)) throw CalcError(CalcErrorType::Syntax, errorMessage(CalcErrorType::Syntax), cur.pos);
 }

 bool Parser::startsPrimary() const { return cur.type == TokenType::Number || cur.type == TokenType::Identifier || cur.type == TokenType::LParen || cur.type == TokenType::Percent; }
 bool Parser::isFunctionName(const std::string &name) const { return cfg.functions.find(name) != cfg.functions.end(); }
 bool Parser::isConstantName(const std::string &name) const { return constants.find(name) != constants.end(); }

 bool isValueEnd(TokenType t) {
  switch (t) {
   case TokenType::Number:
   case TokenType::Identifier:
   case TokenType::RParen:
   case TokenType::RBracket:
   case TokenType::Bang: return true;
   default: return false;
  }
 }

 bool isValueStart(TokenType t) {
  switch (t) {
   case TokenType::Number:
   case TokenType::Identifier:
   case TokenType::LParen:
   case TokenType::LBracket:
   case TokenType::Percent: return true;
   default: return false;
  }
 }

 bool Parser::isImplicitMul(const Token &prev, const Token &cur) const {
  // 右が「値の開始」でなければダメ
  auto isValueStart = [&](TokenType t) { return t == TokenType::Number || t == TokenType::Identifier || t == TokenType::LParen; };
  // 左が「値の終了」でなければダメ
  auto isValueEnd = [&](TokenType t) { return t == TokenType::Number || t == TokenType::Identifier || t == TokenType::RParen || t == TokenType::Bang; };
  if (!isValueEnd(prev.type) || !isValueStart(cur.type)) return false;
  // Identifier + Number → 定数のみ許可（sin30 みたいなのはルカルカ★ナイトフィーバー）
  if (prev.type == TokenType::Identifier && cur.type == TokenType::Number) { return isConstantName(prev.text); }
  // Identifier + '('
  if (prev.type == TokenType::Identifier && cur.type == TokenType::LParen) {
   // 関数なら implicit mul しない（= 関数呼び出しに任せる）
   return !isFunctionName(prev.text);
  }
  return true;
 }

 std::unique_ptr<Parser::ASTNode> Parser::parseExpression() {
  auto n = parseTerm();
  while (cur.type == TokenType::Plus || cur.type == TokenType::Minus) {
   BinOp op = tokenToBinOp(cur.type);
   size_t p = cur.pos;
   advance();
   auto r = parseTerm();
   n = std::make_unique<BinaryNode>(op, std::move(n), std::move(r), p);
  }
  return n;
 }

 std::unique_ptr<Parser::ASTNode> Parser::parseTerm() {
  auto n = parseUnary();

  while (true) {
   if (cur.type == TokenType::Mul || cur.type == TokenType::Div) {
    BinOp op = tokenToBinOp(cur.type);
    size_t p = cur.pos;
    advance();
    auto r = parseUnary();
    n = std::make_unique<BinaryNode>(op, std::move(n), std::move(r), p);
   }
   // 暗黙乗算
   else if (isImplicitMul(prev, cur)) {
    size_t p = cur.pos;
    auto r = parseUnary();
    n = std::make_unique<BinaryNode>(BinOp::Mul, std::move(n), std::move(r), p);
   } else {
    break;
   }
  }
  return n;
 }

 std::unique_ptr<Parser::ASTNode> Parser::parsePower() {
  auto n = parsePostfix();
  if (cur.type == TokenType::Pow) {
   size_t p = cur.pos;
   advance();
   auto r = parseUnary();
   // auto r = parsePower(); // 厳格な右辺結合ならpsersePowerがいいが[2^-2]の-が読めなくなる
   n = std::make_unique<BinaryNode>(BinOp::Pow, std::move(n), std::move(r), p);
  }
  return n;
 }
 std::unique_ptr<Parser::ASTNode> Parser::parseUnary() {
  if (cur.type == TokenType::Plus || cur.type == TokenType::Minus) {
   UnaryOp op = (cur.type == TokenType::Minus) ? UnaryOp::Minus : UnaryOp::Plus;
   size_t p = cur.pos;
   advance();
   return std::make_unique<UnaryNode>(op, parseUnary(), p);
  }
  return parsePower();
 }
 std::unique_ptr<Parser::ASTNode> Parser::parsePostfix() {
  auto node = parsePrimary();

  while (cur.type == TokenType::Bang) {
   size_t p = cur.pos;
   advance();
   node = std::make_unique<PostfixUnaryNode>('!', std::move(node), p);
  }

  return node;
 }

 std::unique_ptr<Parser::ASTNode> Parser::parsePrimary() {
  size_t p = cur.pos;

  // ---- number ----
  if (cur.type == TokenType::Number) {
   double v = cur.number;
   advance();
   return std::make_unique<NumberNode>(v, p);
  }

  // ---- history % ----
  if (cur.type == TokenType::Percent) {
   int c = 0;
   while (cur.type == TokenType::Percent) {
    advance();
    ++c;
   }
   return std::make_unique<OutNode>(-c, p);
  }

  // ---- identifier ----
  if (cur.type == TokenType::Identifier) {
   std::string name = cur.text;
   advance();

   // ----- followed by ( or [ -----
   if (cur.type == TokenType::LParen || cur.type == TokenType::LBracket) {
    TokenType open = cur.type;
    size_t callPos = cur.pos;
    advance();

    // 関数なら普通に call
    if (cfg.functions.count(name)) {
     auto f = std::make_unique<FunctionCallNode>();
     f->name = name;
     f->pos = p;

     // no args
     if ((open == TokenType::LParen && accept(TokenType::RParen)) || (open == TokenType::LBracket && accept(TokenType::RBracket))) { return f; }

     while (true) {
      f->args.push_back(parseExpression());
      if (accept(TokenType::Comma)) continue;
      if (open == TokenType::LParen) expect(TokenType::RParen);
      else expect(TokenType::RBracket);
      break;
     }
     return f;
    }

    // ★ 定数なら「定数 * (expr)」
    if (constants.count(name)) {
     auto lhs = std::make_unique<NumberNode>(constants.at(name), p);
     auto rhs = parseExpression();
     if (open == TokenType::LParen) expect(TokenType::RParen);
     else expect(TokenType::RBracket);
     return std::make_unique<BinaryNode>(BinOp::Mul, std::move(lhs), std::move(rhs), callPos);
    }

    throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), p);
   }

   // ----- identifier only -----
   if (constants.count(name)) return std::make_unique<NumberNode>(constants.at(name), p);

   throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), p);
  }

  // ---- grouped expression ----
  if (cur.type == TokenType::LParen || cur.type == TokenType::LBracket) {
   TokenType open = cur.type;
   advance();

   auto n = parseExpression();

   if (open == TokenType::LParen) expect(TokenType::RParen);
   else expect(TokenType::RBracket);

   return n;
  }

  throw CalcError(CalcErrorType::Syntax, errorMessage(CalcErrorType::Syntax), cur.pos);
 }

 // <,>,<=,>=,==,!=
 std::unique_ptr<Parser::ASTNode> Parser::parseCompare() {
  auto n = parseExpression();
  while (cur.type == TokenType::Less || cur.type == TokenType::LessEq || cur.type == TokenType::Greater || cur.type == TokenType::GreaterEq || cur.type == TokenType::Equal || cur.type == TokenType::NotEqual) {
   CmpOp op = tokenToCmpOp(cur.type);
   size_t p = cur.pos;
   advance();
   auto r = parseExpression();
   n = std::make_unique<CompareNode>(op, std::move(n), std::move(r), p);
  }
  return n;
 }

 /* ============================
    AST ヘルパくん
    ============================ */

 Parser::CmpOp Parser::tokenToCmpOp(TokenType t) {
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
 Parser::BinOp Parser::tokenToBinOp(TokenType t) {
  switch (t) {
   case TokenType::Plus: return BinOp::Add;
   case TokenType::Minus: return BinOp::Sub;
   case TokenType::Mul: return BinOp::Mul;
   case TokenType::Div: return BinOp::Div;
   case TokenType::Pow: return BinOp::Pow;
   default: throw CalcError(CalcErrorType::InvalidOperation, errorMessage(CalcErrorType::InvalidOperation), 0);
  }
 }

 Value evalCompare(const Value &lhs, const Value &rhs, Parser::CmpOp op, FunctionContext &ctx) {
  // 複素数は大小比較不可
  if (isComplex(lhs) || isComplex(rhs)) throw CalcError(CalcErrorType::NotImplemented, "complex comparison is not implemented", ctx.pos);
  double a = asDouble(lhs, ctx.pos), b = asDouble(rhs, ctx.pos);
  bool r = false;
  int prec = ctx.cfg.precision;

  switch (op) {
   case Parser::CmpOp::Equal: r = tolerantEqual(a, b, prec); break;
   case Parser::CmpOp::NotEqual: r = !tolerantEqual(a, b, prec); break;
   case Parser::CmpOp::Less: r = a < b && !tolerantEqual(a, b, prec); break;
   case Parser::CmpOp::LessEq: r = a < b || tolerantEqual(a, b, prec); break;
   case Parser::CmpOp::Greater: r = a > b && !tolerantEqual(a, b, prec); break;
   case Parser::CmpOp::GreaterEq: r = a > b || tolerantEqual(a, b, prec); break;
  }

  return r ? 1.0 : 0.0;
 }

 /* ============================
    AST
    ============================ */

 // AST ノード系
 Parser::NumberNode::NumberNode(Value v, size_t p) : value(std::move(v)) { pos = p; }
 Value Parser::NumberNode::evalImpl(SystemConfig &, const std::vector<InputEntry> &, int) const { return value; }
 Parser::BinaryNode::BinaryNode(BinOp o, std::unique_ptr<Parser::ASTNode> l, std::unique_ptr<Parser::ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }

 Value Parser::BinaryNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {
  auto a = lhs->eval(cfg, hist, base);
  auto b = rhs->eval(cfg, hist, base);

  // ===== complex が含まれる場合 =====
  if (isComplex(a) || isComplex(b)) {
   Complex x = toComplex(a);
   Complex y = toComplex(b);

   // nan注意
   if (!std::isfinite(x.real()) || !std::isfinite(x.imag()) || !std::isfinite(y.real()) || !std::isfinite(y.imag())) { throw CalcError(CalcErrorType::DomainError, "nan or inf", pos); }

   switch (op) {
    case BinOp::Add: return x + y;
    case BinOp::Sub: return x - y;
    case BinOp::Mul: return x * y;
    case BinOp::Div:
     if (std::abs(y) == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
     return x / y;
    case BinOp::Pow: return std::pow(x, y);
   }
  }

  // ===== real 同士 =====
  double x = asDouble(a, pos);
  double y = asDouble(b, pos);

  // nan注意（real）
  if (!std::isfinite(x) || !std::isfinite(y)) throw CalcError(CalcErrorType::DomainError, "nan or inf", pos);

  switch (op) {
   case BinOp::Add: return x + y;
   case BinOp::Sub: return x - y;
   case BinOp::Mul: return x * y;
   case BinOp::Div:
    if (y == 0.0) throw CalcError(CalcErrorType::DivisionByZero, errorMessage(CalcErrorType::DivisionByZero), pos);
    return x / y;
   case BinOp::Pow:
    // real→complex 昇格判定はここだけ
    if (x < 0.0 && std::floor(y) != y) {
     // 警告（stderr）
     calcWarn(cfg, pos, "(-a)^(p/q) : principal value only; other branches may exist");
     return std::pow(Complex(x, 0.0), Complex(y, 0.0));
    }
    return std::pow(x, y);
  }

  throw CalcError(CalcErrorType::InvalidOperation, "invalid op", pos);
 }

 Parser::UnaryNode::UnaryNode(UnaryOp o, std::unique_ptr<Parser::ASTNode> r, size_t p) : op(o), rhs(std::move(r)) { pos = p; }
 Value Parser::UnaryNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {
  auto v = rhs->eval(cfg, hist, base);
  if (op == UnaryOp::Minus) {
   if (isComplex(v)) return -toComplex(v);
   return -asDouble(v, pos);
  }
  return v;
 }

 Parser::PostfixUnaryNode::PostfixUnaryNode(char op, std::unique_ptr<Parser::ASTNode> e, size_t p) : op(op), expr(std::move(e)) { pos = p; }

 Value Parser::PostfixUnaryNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {
  Value v = expr->eval(cfg, hist, base);
  // 後置演算子は実数のみ対応
  if (!std::holds_alternative<double>(v)) { throw CalcError(CalcErrorType::TypeError, errorMessage(CalcErrorType::TypeError), pos); }
  double d = std::get<double>(v), r;
  switch (op) {
   case '!': r = factorial(d, pos); break;
   default: throw CalcError(CalcErrorType::Syntax, "unknown postfix operator", pos);
  }
  return r;
 }

 Parser::CompareNode::CompareNode(CmpOp o, std::unique_ptr<Parser::ASTNode> l, std::unique_ptr<Parser::ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }
 Value Parser::CompareNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {
  Value a = lhs->eval(cfg, hist, base);
  Value b = rhs->eval(cfg, hist, base);

  FunctionContext ctx{cfg, hist, base, pos};
  return evalCompare(a, b, op, ctx);
 }

 Value Parser::FunctionCallNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {

  auto it = cfg.functions.find(name);
  if (it == cfg.functions.end()) throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), pos);

  auto &f = it->second;
  int argc = static_cast<int>(args.size());
  if (!f.validArgc(argc)) throw CalcError(CalcErrorType::FunctionMissing, errorMessage(CalcErrorType::FunctionMissing), pos);
  std::vector<Value> v;
  v.reserve(args.size());
  for (auto &a : args)
   v.push_back(a->eval(cfg, hist, base));

  FunctionContext ctx{cfg, hist, base, pos};
  Value r = f.f(v, ctx);
  // double のときだけ有限チェック -> 一括捕捉に変更
  /* if (std::holds_alternative<double>(r)) checkFinite(std::get<double>(r), pos);*/
  return r;
 }

 Parser::OutNode::OutNode(int idx, size_t p) : index(idx) { pos = p; }

 Value Parser::OutNode::evalImpl(SystemConfig &, const std::vector<InputEntry> &hist, int) const {
  if (index == 0) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos);

  int i;
  if (index > 0) i = index;
  else i = static_cast<int>(hist.size()) + index + 1;

  if (i <= 0 || i > static_cast<int>(hist.size())) throw CalcError(CalcErrorType::OutOfRange, errorMessage(CalcErrorType::OutOfRange), pos);

  return hist[i - 1].value;
 }

} // namespace mm::cal