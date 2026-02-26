#include "core.hpp"
#include "evaluate.hpp"
#include "functions.hpp"

#include <iostream>
#include <string>

using namespace mm::cal;

/* ============================
  main
  ============================ */

int main(int argc, char *argv[]) {
 SystemConfig syscfg;
 RuntimeState rtmstt;
 std::vector<InputEntry> history;

 initFunctions(syscfg);

 // --batchモード(tester用)
 if (argc > 1 && std::string(argv[1]) == "--batch") {
  std::cout << std::unitbuf;

  std::ios::sync_with_stdio(false);
  std::cin.tie(nullptr);

  std::string line;

  while (std::getline(std::cin, line)) {
   if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

   try {
    auto res = evalLine(line, syscfg, rtmstt, history);
    std::cout << formatResult(res.value, syscfg) << '\n' << std::flush;
   } catch (const std::exception &e) { std::cout << "Error: " << e.what() << '\n'; }
  }

  return 0;
 }

 // ================================
 // CLI モード
 // ================================
 if (argc >= 2) {
  std::string line = argv[1];

  try {
   EvalResult res = evalLine(line, syscfg, rtmstt, history);

   if (res.kind == EvalKind::Value) { std::cout << formatResult(res.value, syscfg) << "\n"; }
   // Clear / None / Exit は何も出さず終了
  } catch (const CalcError &e) {
   std::cout << "Error: " << e.what() << "\n" << line << "\n" << std::string(e.pos, ' ') << "^\n";
   return 1;
  }

  return 0;
 }

 // ================================
 // REPL モード
 // ================================
 std::cout << "================================\n"
              "  mm Calculator\n"
              "        mmKreutzef 2021-2026\n"
              "================================\n\n";

 std::string line;
 std::cout << std::unitbuf;
 while (true) {
  std::cout << "In [" << history.size() + 1 << "] := " << std::flush;
  if (!std::getline(std::cin, line)) break;

  if (line.find_first_not_of(" \t\r\n") == std::string::npos) continue;

  // ---- メタコマンド ----
  if (line.starts_with("SetFix")) {
   syscfg.precision = std::stoi(line.substr(6));
   continue;
  }

  try {
   EvalResult res = evalLine(line, syscfg, rtmstt, history);

   switch (res.kind) {
    case EvalKind::Clear: history.clear(); break;

    case EvalKind::Exit: return 0;

    case EvalKind::Value:
     history.push_back({line, res.value});
     std::cout << "Out[" << history.size() << "] := " << formatResultMulti(res.value, syscfg) << "\n\n";
     break;

    case EvalKind::None: break;
   }
  } catch (const CalcError &e) { std::cout << "\nError: " << e.what() << "\n" << line << "\n" << std::string(e.pos, ' ') << "^\n\n"; }
 }
 std::cout << "bye...nara\n";
}