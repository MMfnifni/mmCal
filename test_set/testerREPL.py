import subprocess
import sys
import re
import math
import glob
import time
import datetime

RED = "\033[31m"
ORANGE = "\033[38;5;208m"   # 256色モード
RESET = "\033[0m"
EPS = 1e-10
EXE = r"..\build\x64\Release\mmCal_x64.exe"
TESTS = sys.argv[1:] or glob.glob("test*.txt")

# ================= Floating Point Helpers =================

def is_quoted_string(s: str) -> bool:
    s = s.strip()
    return len(s) >= 2 and s[0] == '"' and s[-1] == '"'

def unquote_string(s: str) -> str:
    s = s.strip()
    return s[1:-1]

def normalize_float(x, eps=EPS):
    if math.isinf(x):
        return x
    return 0.0 if abs(x) < eps else x

def almost_equal(a, b, abs_eps=EPS, rel_eps=1e-10):
    a, b = normalize_float(a, abs_eps), normalize_float(b, abs_eps)

    if math.isinf(a) or math.isinf(b):
        return a == b

    diff = abs(a - b)
    scale = max(1.0, abs(a), abs(b))
    return diff <= max(abs_eps, rel_eps * scale)

def parse_atom(s: str):
    s = s.strip()

    try:
        return ("real", parse_float(s))
    except Exception:
        pass

    try:
        re_part, im_part = parse_complex(s)
        # complex表記らしいものだけ complex として採用
        if looks_complex(s):
            return ("complex", (re_part, im_part))
    except Exception:
        pass

    return ("text", s)

def parse_float(s):
    s = s.strip().lower()
    if s in ("inf", "+inf"): return float("inf")
    if s == "-inf": return float("-inf")
    return float(s)

# ================= Complex Helpers =================

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
        return float(s), 0.0

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

# ================= Structure Parser =================

def parse_structure(s: str):
    s = s.strip()

    if not s.startswith("{"):
        kind, value = parse_atom(s)
        return (kind, value)

    if s == "{}":
        return []

    inner = s[1:-1].strip()
    if not inner:
        return []

    elems = []
    depth = 0
    token = ""

    for ch in inner:
        if ch == "," and depth == 0:
            elems.append(parse_structure(token))
            token = ""
        else:
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
            token += ch

    if token:
        elems.append(parse_structure(token))

    return elems

def structure_equal(a, b):
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return False
        return all(structure_equal(x, y) for x, y in zip(a, b))

    if isinstance(a, tuple) and isinstance(b, tuple):
        ka, va = a
        kb, vb = b

        # real-real
        if ka == "real" and kb == "real":
            return almost_equal(float(va), float(vb))

        # complex-complex
        if ka == "complex" and kb == "complex":
            ar, ai = va
            br, bi = vb
            return almost_equal(ar, br) and almost_equal(ai, bi)

        # real-complex / complex-real は虚部ゼロなら同値扱い
        if ka == "real" and kb == "complex":
            br, bi = vb
            return almost_equal(float(va), br) and almost_equal(0.0, bi)

        if ka == "complex" and kb == "real":
            ar, ai = va
            return almost_equal(ar, float(vb)) and almost_equal(ai, 0.0)

        # text-text
        if ka == "text" and kb == "text":
            return va == vb

        return False

    return False

def normalize_structure(s):
    s = s.strip()
    s = re.sub(r'\s*([{},])\s*', r'\1', s)
    return s

# ================= Batch Control =================

def start_batch():
    return subprocess.Popen(
        [EXE, "--batch"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1
    )

def evaluate_expr(p, expr):
    p.stdin.write(expr + "\n")
    p.stdin.flush()
    line = p.stdout.readline()
    if not line:
        return "", 1
    line = line.rstrip()
    if line.startswith("Error:"):
        return line, 1
    return line, 0

# ================= Test Runner =================

def run_test_batch(p, expr, expect):
    out, ret = evaluate_expr(p, expr)

    if expect == "Error":
        return ("pass", out) if ret != 0 else ("fail", out or "(empty)")

    if ret != 0:
        return "fail", out

    # 明示的な文字列期待値: "..."
    if is_quoted_string(expect):
        expected_text = unquote_string(expect).strip()
        got_text = out.strip()
        return ("pass", out) if got_text == expected_text else ("fail", out or "(empty)")

    try:
        norm_out = normalize_structure(out)
        norm_expect = normalize_structure(expect)

        if norm_out.startswith("{") or norm_expect.startswith("{"):
            parsed_out = parse_structure(norm_out)
            parsed_expect = parse_structure(norm_expect)
            return ("pass", out) if structure_equal(parsed_out, parsed_expect) else ("fail", out)

        if looks_complex(norm_out) or looks_complex(norm_expect):
            return ("pass", out) if complex_equal(norm_out, norm_expect) else ("fail", out)

        return ("pass", out) if almost_equal(parse_float(norm_out), parse_float(norm_expect)) else ("fail", out)

    except Exception:
        # 最後の保険: 生文字列比較
        return ("pass", out) if out.strip() == expect.strip() else ("fail", out or "(empty)")

# ================= Main =================

PASS = FAIL = TOTAL = 0
now = datetime.datetime.now() 
today = now.strftime("'%y/%m/%d/%H:%M:%S") 
print(today+" start\nRunning tests...\n")
start = time.perf_counter()


p = start_batch()

for testfile in TESTS:
    with open(testfile, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            TOTAL += 1

            if "=>" not in line:
                print(f"{ORANGE}[SKIP]{RESET} {line}")
                continue

            expr, expect = map(str.strip, line.split("=>", 1))

            if expect.lower() == "skip":
                print(f"{ORANGE}[SKIP]{RESET} {expr}")
                continue

            result, got = run_test_batch(p, expr, expect)

            if result == "pass":
                print(f"[PASS] {expr} => {got}" if got else f"[PASS] {expr}")
                PASS += 1
            else:
                print(f"{RED}[FAIL]{RESET} {expr}")
                print(f"  expected: {expect}")
                print(f"  got     : {got}")
                FAIL += 1

end = time.perf_counter()
elapsed = (end - start)* 1000
print("\n=====================")
print(f"TOTAL: {TOTAL}")
print(f"PASS : {PASS}")
print(f"FAIL : {FAIL}")
print(f"\nElapsed time: {elapsed:.3f} [ms]")
print("=====================")

p.stdin.close()
p.terminate()
sys.exit(FAIL != 0)