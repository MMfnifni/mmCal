#include "lexer_parser.hpp"
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
   throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: '=' is not allowed; use ':=' for assignment/definition", start);
  }

  if (c == ':') {
   if (p < src.size() && src[p] == '=') {
    ++p;
    return {TokenType::Assign, ":=", start};
   }
   throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: unexpected ':'", start);
  }

  if (c == '"') {
   std::string s;
   while (p < src.size() && src[p] != '"') {
    // エスケープ無しでOK
    s.push_back(src[p]);
    ++p;
   }

   if (p >= src.size()) { throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: unterminated string", start); }

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
   case '{': return {TokenType::LBrace, "{", start};
   case '}': return {TokenType::RBrace, "}", start};
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

 bool Parser::startsPrimary() const {
  return cur.type == TokenType::Number || cur.type == TokenType::String || cur.type == TokenType::Identifier || cur.type == TokenType::LParen || cur.type == TokenType::LBracket || cur.type == TokenType::LBrace || cur.type == TokenType::Percent;
 }
 bool Parser::isFunctionName(const std::string &name) const { return cfg.functions.find(name) != cfg.functions.end(); }
 bool Parser::isConstantName(const std::string &name) const { return constants.find(name) != constants.end(); }
 bool Parser::isUnitName(const std::string &s) const { return symbols.count(s) != 0; }

 bool Parser::isValueEnd(TokenType t) {
  switch (t) {
   case TokenType::Number:
   case TokenType::Identifier:
   case TokenType::RParen:
   case TokenType::RBracket:
   case TokenType::RBrace:
   case TokenType::Bang: return true;
   default: return false;
  }
 }

 bool Parser::isValueStart(TokenType t) {
  switch (t) {
   case TokenType::Number:
   case TokenType::Identifier:
   case TokenType::LParen:
   case TokenType::LBracket:
   case TokenType::LBrace:
   case TokenType::Percent: return true; // ←ここは方針次第なので要検討
   default: return false;
  }
 }

 bool Parser::isImplicitMul(const Token &prev, const Token &cur) const { // ちなみにMultiValueは暗黙乗算の対象外
  // 右が「値の開始」でなければダメ
  auto isValueStart = [&](TokenType t) { return t == TokenType::Number || t == TokenType::Identifier || t == TokenType::LParen; };
  // 左が「値の終了」でなければダメダメよ
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

 std::unique_ptr<ASTNode> Parser::parse() {
  auto node = parseAssignment(); // パースの入口が変わったここを変えること(どうせ私は忘れてる)
  if (cur.type != TokenType::End) throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: unexpected token", cur.pos);
  return node;
 }

 std::unique_ptr<ASTNode> Parser::parseAssignment() {
  auto left = parseCompare();

  if (cur.type == TokenType::Assign) {
   size_t assignPos = cur.pos;
   advance();
   auto rhs = parseAssignment();

   // 変数代入: x = expr
   if (auto ident = dynamic_cast<SymbolNode *>(left.get())) { return std::make_unique<AssignNode>(ident->name, std::move(rhs), ident->pos); }

   // 関数定義: f(x, y) = expr
   if (auto call = dynamic_cast<FunctionCallNode *>(left.get())) {
    std::vector<std::string> params;
    params.reserve(call->args.size());

    for (const auto &arg : call->args) {
     auto sym = dynamic_cast<SymbolNode *>(arg.get());
     if (!sym) { throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: function parameters must be identifiers", assignPos); }
     params.push_back(sym->name);
    }

    return std::make_unique<FunctionDefNode>(call->name, std::move(params), std::move(rhs), call->pos);
   }

   throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: left side of assignment must be identifier or function signature", assignPos);
  }

  return left;
 }

 std::unique_ptr<ASTNode> Parser::parseExpression() {
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

 std::unique_ptr<ASTNode> Parser::parseTerm() {
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

    // unit の場合は UnitApplyNode
    if (cur.type == TokenType::Identifier && symbols.count(cur.text)) {

     // unitを付けていい左辺を制限する
     // OK: 30deg, 30 rad, PI rad
     // NG: (30)rad, sin(30)rad, xrad
     bool okLhs = false;

     // 直前トークンが Number ならOK
     if (prev.type == TokenType::Number) okLhs = true;

     // 直前トークンが Identifier の場合は「定数だけOK」
     if (prev.type == TokenType::Identifier) { okLhs = isConstantName(prev.text); }

     if (!okLhs) { throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: nit can follow only a number or constant (e.g. 30deg, PI rad)", p); }

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

 std::unique_ptr<ASTNode> Parser::parsePower() {
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
 std::unique_ptr<ASTNode> Parser::parseUnary() {
  if (cur.type == TokenType::Plus || cur.type == TokenType::Minus) {
   UnaryOp op = (cur.type == TokenType::Minus) ? UnaryOp::Minus : UnaryOp::Plus;
   size_t p = cur.pos;
   advance();
   return std::make_unique<UnaryNode>(op, parseUnary(), p);
  }
  return parsePower();
 }
 std::unique_ptr<ASTNode> Parser::parsePostfix() {
  auto node = parsePrimary();

  while (cur.type == TokenType::Bang) {
   size_t p = cur.pos;
   advance();
   node = std::make_unique<PostfixUnaryNode>('!', std::move(node), p);
  }

  return node;
 }

 std::unique_ptr<ASTNode> Parser::parsePrimary() {

  // ---- number ----
  if (cur.type == TokenType::Number) {
   size_t p = cur.pos;
   double v = cur.number;
   advance();
   return std::make_unique<NumberNode>(v, p);
  }

  // ---- string ----
  if (cur.type == TokenType::String) {
   size_t p = cur.pos;
   std::string s = cur.text;
   advance();
   return std::make_unique<NumberNode>(Value(std::move(s)), p);
  }

  // ---- history % ----
  if (cur.type == TokenType::Percent) {
   size_t p = cur.pos;
   int c = 0;
   while (cur.type == TokenType::Percent) {
    advance();
    ++c;
   }
   return std::make_unique<OutNode>(-c, p);
  }

  // ----------- {multivalue} -----------
  if (cur.type == TokenType::LBrace) {
   size_t startPos = cur.pos;
   advance();
   std::vector<std::unique_ptr<ASTNode>> elems;
   if (cur.type != TokenType::RBrace) {
    while (true) {
     elems.push_back(parseAssignment());
     if (accept(TokenType::Comma)) continue;
     expect(TokenType::RBrace);
     break;
    }
   } else {
    advance();
   }
   return std::make_unique<MultiLiteralNode>(std::move(elems), startPos);
  }

  // ---- identifier ----
  if (cur.type == TokenType::Identifier) {
   std::string name = cur.text;
   size_t p = cur.pos;
   advance();

   // ----- followed by ( or [ -----
   if (cur.type == TokenType::LParen || cur.type == TokenType::LBracket) {
    TokenType open = cur.type;
    size_t callPos = cur.pos;
    advance();

    // 定数なら「定数 * (expr)」
    if (constants.count(name)) {
     auto lhs = std::make_unique<NumberNode>(constants.at(name), p);
     auto rhs = parseExpression();
     if (open == TokenType::LParen) expect(TokenType::RParen);
     else expect(TokenType::RBracket);
     return std::make_unique<BinaryNode>(BinOp::Mul, std::move(lhs), std::move(rhs), callPos);
    }

    // それ以外は builtin / user function を区別せず関数呼び出しとして構文木化
    auto f = std::make_unique<FunctionCallNode>();
    f->name = name;
    f->pos = p;

    // no args
    if ((open == TokenType::LParen && accept(TokenType::RParen)) || (open == TokenType::LBracket && accept(TokenType::RBracket))) { return f; }

    while (true) {
     f->args.push_back(parseAssignment());

     if (accept(TokenType::Comma)) continue;

     if (open == TokenType::LParen) expect(TokenType::RParen);
     else expect(TokenType::RBracket);
     break;
    }
    return f;
   }

   // ----- identifier only -----
   if (constants.count(name)) { return std::make_unique<NumberNode>(constants.at(name), p); }
   if (symbols.count(name)) { throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: unit must follow a value (e.g. 30deg)", p); }
   return std::make_unique<SymbolNode>(std::move(name), p);
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

  throw CalcError(CalcErrorType::SyntaxError, "SyntaxError: invalid parentheses", cur.pos);
 }

 // <,>,<=,>=,==,!=
 std::unique_ptr<ASTNode> Parser::parseCompare() {
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

} // namespace mm::cal