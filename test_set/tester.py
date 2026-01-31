import subprocess
import sys
import re

EPS = 1e-9

def almost_equal(a, b, eps=EPS):
    return abs(a - b) < eps

def parse_complex(s):
    s = s.strip()

    if s.endswith("I"):
        body = s[:-1]

        if body == "" or body == "+":
            return 0.0, 1.0
        if body == "-":
            return 0.0, -1.0

        if "+" in body[1:]:
            r, i = body.rsplit("+", 1)
            return float(r), float(i)
        elif "-" in body[1:]:
            r, i = body.rsplit("-", 1)
            return float(r), -float(i)
        else:
            return 0.0, float(body)
    else:
        return float(s), 0.0


def complex_equal(a, b):
    ar, ai = parse_complex(a)
    br, bi = parse_complex(b)
    return almost_equal(ar, br) and almost_equal(ai, bi)

EXE = r"..\x64\Release\cal.exe"
TESTS = "test.txt"

PASS = FAIL = TOTAL = 0

print("Running tests...\n")

with open(TESTS, encoding="utf-8") as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#"):
            continue

        TOTAL += 1

        if "=" not in line:
            print(f"[SKIP] {line}")
            continue

        expr, expect = line.split(":=", 1)
        expect = expect.strip()

        try:
            p = subprocess.run(
                [EXE, expr],
                capture_output=True,
                text=True,
                timeout=2
            )
        except Exception:
            print(f"[FAIL] {expr} (runner error)")
            FAIL += 1
            continue

        out = p.stdout.strip()
        ret = p.returncode

        # ---- Error expected ----
        if expect == "Error":
            if ret != 0:
                print(f"[PASS] {expr} (Error)")
                PASS += 1
            else:
                print(f"[FAIL] {expr}")
                print(f"  expected: Error")
                print(f"  got     : {out if out else '(empty)'}")
                FAIL += 1

        # ---- Normal result ----
        else:
            if ret == 0:
                try:
                    if "I" in out or "I" in expect:
                        if complex_equal(out, expect):
                            print(f"[PASS] {expr} = {out}")
                            PASS += 1
                        else:
                            raise ValueError
                    else:
                        if almost_equal(float(out), float(expect)):
                            print(f"[PASS] {expr} = {out}")
                            PASS += 1
                        else:
                            raise ValueError
                except Exception:
                    print(f"[FAIL] {expr}")
                    print(f"  expected: {expect}")
                    print(f"  got     : {out if out else '(empty)'}")
                    FAIL += 1
            else:
                print(f"[FAIL] {expr}")
                print(f"  expected: {expect}")
                print(f"  got     : Error")
                FAIL += 1

print("\n=====================")
print(f"TOTAL: {TOTAL}")
print(f"PASS : {PASS}")
print(f"FAIL : {FAIL}")
print("=====================")

sys.exit(FAIL != 0)
