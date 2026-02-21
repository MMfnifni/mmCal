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
   while (p < src.size() && (std::isalnum((unsigned char)src[p]) || src[p] == '_'))
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

  if (c == '"') {
   std::string s;
   while (p < src.size() && src[p] != '"') {
    // エスケープ無しでOK
    s.push_back(src[p]);
    ++p;
   }

   if (p >= src.size()) { throw CalcError(CalcErrorType::SyntaxError, "unterminated string", start); }

   ++p; // closing "
   return {TokenType::String, s, 0.0, start};
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
   case '%':
    return {TokenType::Percent, "%", start};
    // case '"': return {TokenType::String, '"', start};
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
  if (!accept(t)) throw CalcError(CalcErrorType::SyntaxError, errorMessage(CalcErrorType::SyntaxError), cur.pos);
 }

 bool Parser::startsPrimary() const { return cur.type == TokenType::Number || cur.type == TokenType::String || cur.type == TokenType::Identifier || cur.type == TokenType::LParen || cur.type == TokenType::Percent; }
 bool Parser::isFunctionName(const std::string &name) const { return cfg.functions.find(name) != cfg.functions.end(); }
 bool Parser::isConstantName(const std::string &name) const { return constants.find(name) != constants.end(); }
 bool Parser::isUnitName(const std::string &s) const { return symbols.count(s) != 0; }

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
   case TokenType::Percent: return true; // ←ここは方針次第なので要検討
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

    // ★ unit の場合は UnitApplyNode
    if (cur.type == TokenType::Identifier && symbols.count(cur.text)) {

     // unitを付けていい左辺を制限する
     // OK: 30deg, 30 rad, PI rad
     // NG: (30)rad, sin(30)rad, xrad
     bool okLhs = false;

     // 直前トークンが Number ならOK
     if (prev.type == TokenType::Number) okLhs = true;

     // 直前トークンが Identifier の場合は「定数だけOK」
     if (prev.type == TokenType::Identifier) { okLhs = isConstantName(prev.text); }

     if (!okLhs) { throw CalcError(CalcErrorType::SyntaxError, "unit can follow only a number or constant (e.g. 30deg, PI rad)", p); }

     std::string u = cur.text;
     advance(); // unit消費
     n = std::make_unique<UnitApplyNode>(std::move(n), std::move(u), p);
     continue;
    }

    // 通常の暗黙乗算
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

  // ---- string ----
  if (cur.type == TokenType::String) {
   std::string s = cur.text;
   advance();
   return std::make_unique<NumberNode>(Value(std::move(s)), p);
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
      f->args.push_back(parseCompare());

      if (accept(TokenType::Comma)) continue;

      if (open == TokenType::LParen) expect(TokenType::RParen);
      else expect(TokenType::RBracket);
      break;
     }
     return f;
    }

    // 定数なら「定数 * (expr)」
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
   if (constants.count(name)) { return std::make_unique<NumberNode>(constants.at(name), p); }
   if (symbols.count(name)) { throw CalcError(CalcErrorType::SyntaxError, "unit must follow a value (e.g. 30deg)", p); } // 単位は単体では許可しない（30deg の形だけ）
   return std::make_unique<SymbolNode>(std::move(name), p);                                                              // それ以外は「オプション指定子」としてASTに残す
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

  throw CalcError(CalcErrorType::SyntaxError, errorMessage(CalcErrorType::SyntaxError), cur.pos);
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
  if (lhs.isComplex() || rhs.isComplex()) throw CalcError(CalcErrorType::NotImplemented, "complex comparison is not implemented", ctx.pos);
  double a = lhs.asScalar(ctx.pos), b = rhs.asScalar(ctx.pos);
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
  switch (op) {
   case BinOp::Add: return add(a, b, pos);
   case BinOp::Sub: return sub(a, b, pos);
   case BinOp::Mul: return mul(a, b, pos);
   case BinOp::Div: return div(a, b, pos);
   case BinOp::Pow: return power(a, b, pos);
  }

  throw CalcError(CalcErrorType::InvalidOperation, "invalid op", pos);
 }

 Parser::UnaryNode::UnaryNode(UnaryOp o, std::unique_ptr<Parser::ASTNode> r, size_t p) : op(o), rhs(std::move(r)) { pos = p; }
 Value Parser::UnaryNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {
  auto v = rhs->eval(cfg, hist, base);
  if (op == UnaryOp::Minus) return negate(v, pos);
  return v;
 }

 Parser::PostfixUnaryNode::PostfixUnaryNode(char op, std::unique_ptr<Parser::ASTNode> e, size_t p) : op(op), expr(std::move(e)) { pos = p; }

 Value Parser::PostfixUnaryNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {
  Value v = expr->eval(cfg, hist, base);

  if (op == '!') return mm::cal::factorial(v, pos);

  throw CalcError(CalcErrorType::InvalidOperation, "invalid postfix", pos);
 }

 Parser::CompareNode::CompareNode(CmpOp o, std::unique_ptr<Parser::ASTNode> l, std::unique_ptr<Parser::ASTNode> r, size_t p) : op(o), lhs(std::move(l)), rhs(std::move(r)) { pos = p; }
 Value Parser::CompareNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {

  Value a = lhs->eval(cfg, hist, base);
  Value b = rhs->eval(cfg, hist, base);

  mm::cal::CompareOp cop;

  switch (op) {
   case CmpOp::Equal: cop = CompareOp::Eq; break;
   case CmpOp::NotEqual: cop = CompareOp::Ne; break;
   case CmpOp::Less: cop = CompareOp::Lt; break;
   case CmpOp::LessEq: cop = CompareOp::Le; break;
   case CmpOp::Greater: cop = CompareOp::Gt; break;
   case CmpOp::GreaterEq: cop = CompareOp::Ge; break;
  }

  return mm::cal::compare(a, b, cop, pos);
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

 Parser::SymbolNode::SymbolNode(std::string n, size_t p) : name(std::move(n)) { pos = p; }

 Value Parser::SymbolNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &, int) const {
  if (constants.count(name)) return constants.at(name); // 定数はここで解決できる（parse段階で解決しなくてよくなる）
  if (symbols.count(name)) return name;                 // symbols(deg, rad, mm...) はオプション指定子としては文字列扱いで返す

  // 変数（将来用）
  // if (cfg.variables.count(name)) return cfg.variables.at(name);

  // それ以外も「オプション指定子」として許可するならここで返す(要検討!!!)
  // return name;

  throw CalcError(CalcErrorType::UnknownIdentifier, errorMessage(CalcErrorType::UnknownIdentifier), pos);
 }

 Parser::UnitApplyNode::UnitApplyNode(std::unique_ptr<ASTNode> e, std::string u, size_t p) : expr(std::move(e)), unit(std::move(u)) { pos = p; }

 Value Parser::UnitApplyNode::evalImpl(SystemConfig &cfg, const std::vector<InputEntry> &hist, int base) const {
  Value v = expr->eval(cfg, hist, base);
  double x = v.asScalar(pos);

  // angle
  if (unit == "deg") return x;
  if (unit == "rad") { return x * 180.0 / PI; }
  if (unit == "grad") { return x * 0.9; }

  // length (将来)
  if (unit == "mm" || unit == "cm" || unit == "m" || unit == "inch") { throw CalcError(CalcErrorType::NotImplemented, "length units are not implemented yet", pos); }

  throw CalcError(CalcErrorType::UnknownIdentifier, "unknown unit: " + unit, pos);
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