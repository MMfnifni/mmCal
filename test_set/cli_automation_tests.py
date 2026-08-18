#!/usr/bin/env python3
"""Black-box contract tests for mmCal --eval/--batch automation modes."""

from __future__ import print_function

import argparse
from pathlib import Path
import subprocess
import sys


def run(executable, arguments, input_text="", timeout=20):
    completed = subprocess.run(
        [executable] + list(arguments),
        input=input_text,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        timeout=timeout,
    )
    return (
        completed.returncode,
        completed.stdout.replace("\r\n", "\n").replace("\r", "\n"),
        completed.stderr.replace("\r\n", "\n").replace("\r", "\n"),
    )


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def clean_stdout(stdout):
    require("mm Calculator" not in stdout, "automation stdout contains the banner")
    require("In [" not in stdout, "automation stdout contains a prompt")
    require("Out[" not in stdout, "automation stdout contains an interactive label")
    require("bye..nara" not in stdout, "automation stdout contains the farewell")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", required=True, help="path to mmCal executable")
    args = parser.parse_args()
    executable = str(Path(args.exe).expanduser().resolve())
    require(Path(executable).is_file(), "mmCal executable was not found")

    checked = 0

    code, stdout, stderr = run(executable, ["--eval", "1+2*3"])
    clean_stdout(stdout)
    require((code, stdout, stderr) == (0, "7\n", ""), "--eval success contract changed")
    checked += 1

    code, stdout, stderr = run(executable, ["--angle", "deg", "--eval", "sin[30]"])
    clean_stdout(stdout)
    require(code == 0 and stdout == "1/2\n" and stderr == "", "--eval lost startup options")
    checked += 1

    code, stdout, stderr = run(executable, ["--eval", "1+"])
    clean_stdout(stdout)
    require(code == 3 and stdout == "", "syntax failure must use exit code 3 and empty stdout")
    require(stderr.startswith("SyntaxError:"), "syntax failure was not sent to stderr")
    checked += 1

    code, stdout, stderr = run(executable, ["--eval", ""])
    clean_stdout(stdout)
    require(code == 3 and stdout == "", "an empty --eval expression must be a syntax failure")
    require(stderr.startswith("SyntaxError:"), "empty --eval diagnostics were not sent to stderr")
    checked += 1

    code, stdout, stderr = run(executable, ["--eval", "D[f[x],x]"])
    clean_stdout(stdout)
    require(code == 0 and stdout == "D[f[x], x]\n", "warning-bearing result changed")
    require(stderr.startswith("WARN:"), "warning was not isolated on stderr")
    checked += 1

    code, stdout, stderr = run(executable, ["--eval", "1/0"])
    clean_stdout(stdout)
    require(code == 4 and stdout == "", "evaluation failure must use exit code 4")
    require("Error:" in stderr, "evaluation failure was not sent to stderr")
    checked += 1

    code, stdout, stderr = run(executable, ["--eval"])
    clean_stdout(stdout)
    require(code == 2 and stdout == "", "argument failure must use exit code 2")
    require(stderr.startswith("Argument error:"), "argument failure was not sent to stderr")
    checked += 1

    payload = "x:=2\nx^3\n1+\n4\n"
    code, stdout, stderr = run(executable, ["--batch"], payload)
    clean_stdout(stdout)
    require(code == 3, "batch must retain the highest syntax exit code")
    require(stdout == "2\n8\n4\n", "batch did not preserve state or continue after an error")
    require(stderr.startswith("SyntaxError:"), "batch syntax failure was not sent to stderr")
    checked += 1

    code, stdout, stderr = run(executable, ["--bach"], "2+5\n")
    clean_stdout(stdout)
    require((code, stdout, stderr) == (0, "7\n", ""), "--bach compatibility alias changed")
    checked += 1

    hostile = "-" * 50000 + "1\n"
    code, stdout, stderr = run(executable, ["--batch"], hostile)
    clean_stdout(stdout)
    require(code == 3 and stdout == "", "hostile parser input did not fail safely")
    require(stderr.startswith("ResourceLimitError:"), "resource limit is not externally identifiable")
    checked += 1

    code, stdout, stderr = run(executable, ["--batch"], "1+\n1/0\n5\n")
    clean_stdout(stdout)
    require(code == 4 and stdout == "5\n", "batch did not aggregate exit severity or continue")
    require("SyntaxError:" in stderr and "Error:" in stderr, "batch lost one of multiple diagnostics")
    checked += 1

    print("CLI automation: PASS ({} contracts)".format(checked))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, subprocess.TimeoutExpired) as error:
        print("CLI automation: FAIL: {}".format(error), file=sys.stderr)
        sys.exit(1)
