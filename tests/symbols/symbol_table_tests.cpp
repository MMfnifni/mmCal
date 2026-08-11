// 文字列とSymbolIDのintern管理の回帰テスト
#include "symbol_table_tests.hpp"

#include "expression/symbol.hpp"
#include "symbols/symbol_table.hpp"
#include "test_framework.hpp"

#include <string>

namespace mmcal::tests {

void runSymbolTableTests(TestRunner& tests) {
    using expression::Symbol;
    using symbols::SymbolTable;

    SymbolTable table;
    const Symbol x1 = table.intern("x");
    const Symbol x2 = table.intern("x");
    const Symbol y = table.intern("y");

    tests.expect(x1.valid() && x1.id().valid(),
        "SymbolTable: intern returns a valid symbol");
    tests.expectEqual(x1.name(), std::string{"x"},
        "SymbolTable: preserves the source name");
    tests.expect(x1.sameIdentity(x2) && x1.id() == x2.id(),
        "SymbolTable: repeated intern reuses the same identity");
    tests.expect(!x1.sameIdentity(y) && x1.id() != y.id(),
        "SymbolTable: different names receive different identities");
    tests.expectEqual(table.size(), std::size_t{2},
        "SymbolTable: counts unique names only");
    tests.expect(table.find("x").sameIdentity(x1),
        "SymbolTable: find returns the interned symbol");
    tests.expect(!table.find("missing").valid(),
        "SymbolTable: missing lookup returns an invalid symbol");

    SymbolTable otherTable;
    const Symbol otherX = otherTable.intern("x");
    tests.expect(x1 == otherX,
        "Symbol: equal names are semantically equal across tables");
    tests.expect(!x1.sameIdentity(otherX) && x1.id() != otherX.id(),
        "Symbol: intern identity remains table-specific");

    const Symbol defaultX1{"standalone_x"};
    const Symbol defaultX2{"standalone_x"};
    tests.expect(defaultX1.sameIdentity(defaultX2),
        "Symbol: standalone construction uses the default intern table");
}

} // namespace mmcal::tests
