import subprocess
import sys
import re
import math
import glob

EPS = 1e-12
EXE = r"..\build\x64\Release\mmCal_x64.exe"
TESTS = sys.argv[1:] or glob.glob("test*.txt")

# ---------------- Floating Point Helpers ----------------
def normalize_float(x, eps=EPS):
    if math.isinf(x):
        return x
    return 0.0 if abs(x) < eps else x

def almost_equal(a, b, eps=EPS):
    a, b = normalize_float(a, eps), normalize_float(b, eps)
    return a == b if math.isinf(a) or math.isinf(b) else abs(a - b) < eps

def parse_float(s):
    s = s.strip().lower()
    if s in ("inf", "+inf"): return float("inf")
    if s == "-inf": return float("-inf")
    return float(s)

# ---------------- Complex Helpers ----------------
_complex_re = re.compile(r"""
^\s*
(?:
  (?P<real>[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?
  (?:
    (?P<imag_sign>[+-])?
    (?P<imag>(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?
    [iI]
  )?
|
  (?P<pure_sign>[+-])?
  (?P<pure_imag>(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)?
  [iI]
)
\s*$
""", re.VERBOSE)

def parse_complex(s: str):
    s = re.sub(r"\s+", "", s)
    m = _complex_re.match(s)
    if not m:
        return float(s), 0.0  # fallback: real only

    # pure imaginary like "I", "-I", "2I"
    if m.group("pure_imag") or m.group("pure_sign"):
        sign = -1.0 if m.group("pure_sign") == "-" else 1.0
        mag = m.group("pure_imag")
        return 0.0, sign * (float(mag) if mag else 1.0)

    real = float(m.group("real")) if m.group("real") else 0.0
    if not m.group("imag") and not m.group("imag_sign"):
        return real, 0.0

    sign = -1.0 if m.group("imag_sign") == "-" else 1.0
    imag = float(m.group("imag")) if m.group("imag") else 1.0
    return real, sign * imag

def complex_equal(a, b):
    ar, ai = map(normalize_float, parse_complex(a))
    br, bi = map(normalize_float, parse_complex(b))
    return almost_equal(ar, br) and almost_equal(ai, bi)

def looks_complex(s):
    return "i" in s.lower()

# ---------------- Test Runner ----------------
def run_test(expr, expect):
    try:
        p = subprocess.run([EXE, expr], capture_output=True, text=True, timeout=2)
    except Exception:
        return "runner_error", None

    out, ret = p.stdout.strip(), p.returncode

    if expect == "Error":
        return ("pass", out) if ret != 0 else ("fail", out or "(empty)")

    if ret != 0:
        return "fail", "Error"

    try:
        if looks_complex(out) or looks_complex(expect):
            return ("pass", out) if complex_equal(out, expect) else ("fail", out)
        else:
            return ("pass", out) if almost_equal(parse_float(out), parse_float(expect)) else ("fail", out)
    except Exception:
        return "fail", out or "(empty)"

# ---------------- Main ----------------
PASS = FAIL = TOTAL = 0
print("Running tests...\n")

for testfile in TESTS:
    with open(testfile, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            TOTAL += 1

            if ":=" not in line:
                print(f"[SKIP] {line}")
                continue

            expr, expect = map(str.strip, line.split(":=", 1))
            if expect.lower() == "skip":
                print(f"[SKIP] {expr}")
                continue

            result, got = run_test(expr, expect)

            if result == "pass":
                print(f"[PASS] {expr} = {got}" if got else f"[PASS] {expr}")
                PASS += 1
            else:
                print(f"[FAIL] {expr}")
                print(f"  expected: {expect}")
                print(f"  got     : {got}")
                FAIL += 1

print("\n=====================")
print(f"TOTAL: {TOTAL}")
print(f"PASS : {PASS}")
print(f"FAIL : {FAIL}")
print("=====================")

sys.exit(FAIL != 0)
