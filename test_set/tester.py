# mmCal V1.5 black-box test runner.
# Python 3.7+ / Windows and POSIX compatible.

from __future__ import print_function

import argparse
import datetime
from fractions import Fraction
import glob
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time

ABS_EPS = 1e-10
REL_EPS_STRICT = 1e-10
REL_EPS_LOOSE = 1e-6

ERROR_TYPES = (
    "SyntaxError", "ResourceLimitError", "DomainError", "TypeError", "NameError",
    "OverflowError", "EvaluationError", "InternalError"
)
PUBLIC_ERROR_TYPES = tuple(name for name in ERROR_TYPES if name != "InternalError")


class RunnerProtocolError(RuntimeError):
    pass

# V1.5 final is ``In[n]>`` / ``Out[n]>``.  Older black-box targets used
# ``In [n]>`` in some builds, so accept optional whitespace around the index.
# This keeps the runner useful across the historical test targets as requested.
PROMPT_RE = re.compile(r"(?:^|\n)In\s*\[\s*\d+\s*\]>[ \t]*")
OUT_RE = re.compile(r"(?:^|\n)Out\s*\[\s*\d+\s*\]>[ \t]*([^\n]*)")


def script_dir():
    return Path(__file__).resolve().parent


def find_executable(explicit=None):
    if explicit:
        path = Path(explicit).expanduser().resolve()
        if path.is_file():
            return str(path)
        raise RuntimeError("mmCal executable not found: {}".format(path))

    env = os.environ.get("MMCAL_EXE")
    if env:
        path = Path(env).expanduser().resolve()
        if path.is_file():
            return str(path)
        raise RuntimeError("MMCAL_EXE does not exist: {}".format(path))

    base = script_dir()

    # test/_set や test_set 配下からでもプロジェクト本体の最新buildを先に探す。
    # cwdのmmCal.exeは古いコピーである可能性があるため最後に回す。
    search_roots = []
    current = base
    for _ in range(5):
        if current not in search_roots:
            search_roots.append(current)
        if current.parent == current:
            break
        current = current.parent

    candidates = []
    for root in search_roots:
        candidates.extend([
            root / "x64" / "Release" / "mmCal.exe",
            root / "build" / "x64" / "Release" / "mmCal.exe",
            root / "build" / "Release" / "mmCal.exe",
            root / "build" / "mmCal.exe",
            root / "x64" / "Release" / "mmCal",
            root / "build" / "x64" / "Release" / "mmCal",
            root / "build" / "Release" / "mmCal",
            root / "build" / "mmCal",
        ])
    candidates.extend([Path.cwd() / "mmCal.exe", Path.cwd() / "mmCal"])

    seen = set()
    for path in candidates:
        key = str(path)
        if key in seen:
            continue
        seen.add(key)
        if path.is_file():
            return str(path.resolve())

    raise RuntimeError(
        "mmCal executable was not found. Use --exe PATH or set MMCAL_EXE."
    )


def parse_test_sessions(path, text):
    """Parse one physical test file into isolated mmCal sessions.

    ``# @session [ARGS...]`` starts a fresh mmCal process.  This allows related
    regression groups to live in one file without leaking definitions, angle
    mode, history, or diagnostics across the former file boundaries.

    Legacy ``# @args ...`` remains supported.  It sets the startup arguments
    for the current session and must appear before that session's first case.
    """
    sessions = []
    startup_args = []
    cases = []

    def flush_session():
        nonlocal startup_args, cases
        if cases:
            sessions.append((startup_args, cases))
        startup_args = []
        cases = []

    for line_number, raw in enumerate(text.splitlines(), 1):
        line = raw.strip()
        if not line:
            continue

        if line == "# @session" or line.startswith("# @session "):
            flush_session()
            arg_text = line[len("# @session"):].strip()
            startup_args = shlex.split(arg_text, posix=(os.name != "nt")) if arg_text else []
            continue

        if line.startswith("# @args "):
            if cases:
                raise RuntimeError(
                    "{}:{}: # @args must precede cases in a session; use # @session for a new isolated session"
                    .format(path, line_number)
                )
            startup_args = shlex.split(line[len("# @args "):], posix=(os.name != "nt"))
            continue

        if line.startswith("#"):
            continue

        if "==>" in line:
            expr, expected = map(str.strip, line.split("==>", 1))
            mode = "exact"
        elif "=>>" in line:
            expr, expected = map(str.strip, line.split("=>>", 1))
            mode = "loose"
        elif "=>" in line:
            expr, expected = map(str.strip, line.split("=>", 1))
            mode = "strict"
        else:
            raise RuntimeError(
                "{}:{}: expected '==>', '=>' or '=>>'".format(path, line_number)
            )

        cases.append((line_number, expr, expected, mode))

    flush_session()
    return sessions


def parse_test_text(path, text):
    """Backward-compatible single-session parser helper."""
    sessions = parse_test_sessions(path, text)
    if not sessions:
        return [], []
    if len(sessions) != 1:
        raise RuntimeError("{}: contains multiple # @session groups".format(path))
    return sessions[0]


def preload_test_files(test_files):
    loaded = []
    for path in test_files:
        text = path.read_text(encoding="utf-8")
        sessions = parse_test_sessions(path, text)
        loaded.append((path, sessions))
    return loaded


def run_session(executable, startup_args, cases, timeout):
    payload = "\n".join(case[1] for case in cases)
    if payload:
        payload += "\n"

    process = subprocess.run(
        [executable] + startup_args,
        input=payload,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=timeout,
    )

    stdout = process.stdout.replace("\r\n", "\n").replace("\r", "\n")
    stderr = process.stderr.replace("\r\n", "\n").replace("\r", "\n")
    parts = PROMPT_RE.split(stdout)

    # banner is parts[0]. For N inputs the REPL normally emits N+1 prompts,
    # because it prints the next prompt before observing EOF.
    if len(parts) < len(cases) + 1:
        preview = stdout[:1200]
        if len(stdout) > len(preview):
            preview += "\n... (stdout truncated)"
        err_preview = stderr[:600]
        if len(stderr) > len(err_preview):
            err_preview += "\n... (stderr truncated)"
        raise RunnerProtocolError(
            "REPL output could not be segmented: {} cases, {} prompts\n"
            "Check the selected executable and prompt format.\nstdout preview:\n{}\nstderr preview:\n{}"
            .format(len(cases), max(0, len(parts) - 1), preview, err_preview)
        )

    responses = []
    for index in range(len(cases)):
        chunk = parts[index + 1].strip()
        match = OUT_RE.search(chunk)
        if match:
            responses.append((match.group(1).strip(), chunk))
        else:
            responses.append((chunk, chunk))

    return responses, process.returncode, stderr


def unquote(text):
    text = text.strip()
    if len(text) >= 2 and text[0] == '"' and text[-1] == '"':
        body = text[1:-1]
        return bytes(body, "utf-8").decode("unicode_escape")
    return None


def first_error_line(text):
    return text.strip().split("\n", 1)[0].strip()


def first_error_type(text):
    first = first_error_line(text)
    for name in ERROR_TYPES:
        if first.startswith(name + ":"):
            return name
    return None


def parse_real(text):
    text = text.strip()
    low = text.lower()
    if low in ("inf", "+inf"):
        return float("inf")
    if low == "-inf":
        return float("-inf")
    if re.match(r"^[+-]?\d+/\d+$", text):
        return float(Fraction(text))
    return float(text)


_COMPLEX_RE = re.compile(
    r"^([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|[+-]?\d+/\d+)?"
    r"([+-])?"
    r"((?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?|\d+/\d+)?I$"
)


def parse_complex(text):
    s = re.sub(r"\s+", "", text)
    if s in ("I", "+I"):
        return 0.0, 1.0
    if s == "-I":
        return 0.0, -1.0

    # Pure imaginary with explicit magnitude.
    m = re.match(r"^([+-]?)(\d+(?:\.\d*)?|\.\d+|\d+/\d+)(?:[eE]([+-]?\d+))?I$", s)
    if m:
        sign = -1.0 if m.group(1) == "-" else 1.0
        mag_text = m.group(2)
        if m.group(3) is not None:
            mag_text += "e" + m.group(3)
        return 0.0, sign * parse_real(mag_text)

    # Split real and imaginary at the last sign not belonging to an exponent.
    split_at = None
    for index in range(1, len(s) - 1):
        if s[index] in "+-" and s[index - 1] not in "eE":
            split_at = index
    if split_at is not None and s.endswith("I"):
        real_text = s[:split_at]
        imag_text = s[split_at:-1]
        if imag_text in ("+", "-"):
            imag_text += "1"
        return parse_real(real_text), parse_real(imag_text)

    raise ValueError("not a complex literal")


def numeric_equal(a, b, rel_eps):
    if math.isinf(a) or math.isinf(b):
        return a == b
    diff = abs(a - b)
    scale = max(1.0, abs(a), abs(b))
    return diff <= ABS_EPS + rel_eps * scale


def split_top_level(inner):
    result = []
    depth = 0
    token = []
    in_string = False
    escape = False
    for ch in inner:
        if in_string:
            token.append(ch)
            if escape:
                escape = False
            elif ch == "\\":
                escape = True
            elif ch == '"':
                in_string = False
            continue
        if ch == '"':
            in_string = True
            token.append(ch)
        elif ch == "{":
            depth += 1
            token.append(ch)
        elif ch == "}":
            depth -= 1
            token.append(ch)
        elif ch == "," and depth == 0:
            result.append("".join(token).strip())
            token = []
        else:
            token.append(ch)
    result.append("".join(token).strip())
    return result


def compare_structure(got, expected, rel_eps):
    got = got.strip()
    expected = expected.strip()
    if not (got.startswith("{") and got.endswith("}") and
            expected.startswith("{") and expected.endswith("}")):
        return False
    gi = got[1:-1].strip()
    ei = expected[1:-1].strip()
    if not gi or not ei:
        return gi == ei
    gs = split_top_level(gi)
    es = split_top_level(ei)
    if len(gs) != len(es):
        return False
    return all(compare_value(a, b, rel_eps) for a, b in zip(gs, es))


def normalize_symbolic(text):
    return re.sub(r"\s+", "", text.strip())


def compare_value(got, expected, rel_eps):
    got = got.strip()
    expected = expected.strip()

    quoted = unquote(expected)
    if quoted is not None:
        return got == quoted

    if got.startswith("{") or expected.startswith("{"):
        return compare_structure(got, expected, rel_eps)

    try:
        return numeric_equal(parse_real(got), parse_real(expected), rel_eps)
    except Exception:
        pass

    if "I" in got or "I" in expected:
        try:
            try:
                ar, ai = parse_complex(got)
            except Exception:
                ar, ai = parse_real(got), 0.0
            try:
                br, bi = parse_complex(expected)
            except Exception:
                br, bi = parse_real(expected), 0.0
            return numeric_equal(ar, br, rel_eps) and numeric_equal(ai, bi, rel_eps)
        except Exception:
            pass

    return normalize_symbolic(got) == normalize_symbolic(expected)


def compare_response(got, full_chunk, expected, mode):
    expected = expected.strip()
    error_line = first_error_line(full_chunk)
    error_type = first_error_type(full_chunk)

    if expected == "Error":
        ok = error_type in PUBLIC_ERROR_TYPES
        return ok, error_line if error_type is not None else "no error"
    if expected in ERROR_TYPES:
        return error_type == expected, error_line if error_type is not None else "no error"

    # ErrorType: exact message 形式では、分類だけでなく先頭のエラー文も固定する。
    for name in ERROR_TYPES:
        if expected.startswith(name + ":"):
            return error_line == expected, error_line if error_type is not None else "no error"

    if error_type is not None:
        return False, error_line

    if mode == "exact":
        quoted = unquote(expected)
        if quoted is not None:
            return got == quoted, None
        return normalize_symbolic(got) == normalize_symbolic(expected), None

    rel_eps = REL_EPS_LOOSE if mode == "loose" else REL_EPS_STRICT
    return compare_value(got, expected, rel_eps), None


def collect_test_files(arguments):
    if arguments:
        files = []
        for item in arguments:
            matched = glob.glob(item)
            files.extend(matched if matched else [item])
        return [Path(path) for path in files]
    return sorted(script_dir().glob("test*.txt"))


RUNNER_VERSION = "V1.5 exact black-box r4"

def main(argv=None):
    parser = argparse.ArgumentParser(description="mmCal V1.5 black-box tests")
    parser.add_argument("tests", nargs="*", help="test files or glob patterns")
    parser.add_argument("--exe", help="path to mmCal executable")
    parser.add_argument("--timeout", type=float, default=120.0, help="timeout per test file")
    parser.add_argument(
        "--timings", nargs="?", const=10, type=int, default=0, metavar="N",
        help="show the N slowest test files by mmCal process wall time (default: 10)")
    args = parser.parse_args(argv)

    executable = find_executable(args.exe)
    test_files = collect_test_files(args.tests)
    if not test_files:
        print("No test files found.", file=sys.stderr)
        return 2

    # ファイルI/Oとtest記述のparseは計測前に完了させる。
    # Elapsed timeはmmCalの起動・評価・結果比較だけを概ね反映する。
    try:
        loaded_files = preload_test_files(test_files)
    except Exception as exc:
        print("Test preload failed: {}".format(exc), file=sys.stderr)
        return 2

    total = passed = failed = skipped = 0
    active_total = 0
    session_total = 0
    file_timings = []
    for _path, sessions in loaded_files:
        session_total += len(sessions)
        for _startup_args, cases in sessions:
            active_total += sum(1 for case in cases if case[2].lower() != "skip")
            skipped += sum(1 for case in cases if case[2].lower() == "skip")

    print("mmCal {}".format(RUNNER_VERSION))
    print("Executable: {}".format(executable))
    print("Preloaded {} files / {} isolated sessions / {} tests.".format(
        len(loaded_files), session_total, active_total))
    print("Test-file I/O and parsing are excluded from Elapsed time.")
    print(datetime.datetime.now().strftime("'%y/%m/%d/%H:%M:%S") + " start")
    print("Running tests...\n")
    started = time.perf_counter()

    for test_file, sessions in loaded_files:
        file_elapsed_ms = 0.0
        file_active_count = 0

        for session_index, (startup_args, cases) in enumerate(sessions, 1):
            active = [case for case in cases if case[2].lower() != "skip"]
            total += len(active)
            file_active_count += len(active)
            if not active:
                continue

            try:
                session_started = time.perf_counter()
                responses, returncode, stderr = run_session(
                    executable, startup_args, active, args.timeout
                )
                file_elapsed_ms += (time.perf_counter() - session_started) * 1000.0
            except RunnerProtocolError as exc:
                print("[RUNNER ERROR] {} session {}: {}".format(
                    test_file.name, session_index, exc))
                print("Aborted: this is a runner/REPL protocol error, not {} test failures.".format(len(active)))
                return 2
            except Exception as exc:
                print("[RUNNER ERROR] {} session {}: {}".format(
                    test_file.name, session_index, exc))
                print("Aborted before assigning PASS/FAIL to this session.")
                return 2

            for case, response in zip(active, responses):
                line_number, expr, expected, mode = case
                got, chunk = response
                ok, info = compare_response(got, chunk, expected, mode)
                if ok:
                    passed += 1
                    shown = first_error_line(chunk) if first_error_type(chunk) else got
                    print("[PASS] {} => {}".format(expr, shown) if shown else "[PASS] {}".format(expr))
                else:
                    failed += 1
                    print("[FAIL] {}".format(expr))
                    print("  file    : {}:{} (session {})".format(
                        test_file.name, line_number, session_index))
                    print("  expected: {}".format(expected))
                    if first_error_type(chunk):
                        print("  got     : {}".format(first_error_line(chunk)))
                    else:
                        print("  got     : {}".format(got or "(empty)"))
                    if info and info != first_error_line(chunk):
                        print("  error   : {}".format(info))

            if returncode != 0:
                print("[WARN] mmCal exited with code {} ({} session {})".format(
                    returncode, test_file.name, session_index))
            if stderr.strip():
                print("[WARN] stderr ({} session {}):\n{}".format(
                    test_file.name, session_index, stderr.strip()))

        if file_active_count:
            file_timings.append((file_elapsed_ms, test_file.name, file_active_count))

    elapsed_ms = (time.perf_counter() - started) * 1000.0
    print("\n=====================")
    print("TOTAL: {}".format(total))
    print("PASS : {}".format(passed))
    print("FAIL : {}".format(failed))
    if skipped:
        print("SKIP : {}".format(skipped))
    print("\nElapsed time: {:.3f} [ms]".format(elapsed_ms))
    if args.timings > 0 and file_timings:
        print("\nSlowest test files (mmCal process wall time):")
        for elapsed, name, count in sorted(file_timings, reverse=True)[:args.timings]:
            per_test = elapsed / count if count else 0.0
            print("  {:10.3f} ms  {:8.3f} ms/test  {:4d}  {}".format(
                elapsed, per_test, count, name
            ))
    print("=====================")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
