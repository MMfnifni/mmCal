#!/usr/bin/env python3
"""Black-box contract tests for mmCal --eval/--batch automation modes."""

from __future__ import print_function

import argparse
from pathlib import Path
import subprocess
import sys
import signal
import time


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



def run_interrupted(executable, expression, delay=0.05, timeout=10):
    creationflags = 0
    if sys.platform == "win32":
        creationflags = subprocess.CREATE_NEW_PROCESS_GROUP
    process = subprocess.Popen(
        [executable, "--eval", expression],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
        creationflags=creationflags,
    )
    time.sleep(delay)
    require(process.poll() is None, "interrupt probe finished before Ctrl-C could be sent")
    if sys.platform == "win32":
        process.send_signal(signal.CTRL_BREAK_EVENT)
    else:
        process.send_signal(signal.SIGINT)
    stdout, stderr = process.communicate(timeout=timeout)
    return (
        process.returncode,
        stdout.replace("\r\n", "\n").replace("\r", "\n"),
        stderr.replace("\r\n", "\n").replace("\r", "\n"),
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

    code, stdout, stderr = run(executable, ["--help"])
    clean_stdout(stdout)
    require(code == 0 and stderr == "", "--help must succeed without diagnostics")
    require(
        stdout
        == "Usage: mmCal [--fix <0..1000>] [--angle <deg|rad|grad>] [--layout <auto|single|multi>]\n"
        "       mmCal [options] --eval <expression>\n"
        "       mmCal [options] --batch\n"
        "Interactive help: :help [function]\n",
        "--help is no longer the concise startup summary",
    )
    checked += 1

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

    code, stdout, stderr = run(executable, ["--eval", "log[0]"])
    clean_stdout(stdout)
    require(code == 4 and stdout == "", "evaluation failure must use exit code 4")
    require("Error:" in stderr, "evaluation failure was not sent to stderr")
    checked += 1

    code, stdout, stderr = run(executable, ["--eval"])
    clean_stdout(stdout)
    require(code == 2 and stdout == "", "argument failure must use exit code 2")
    require(stderr.startswith("Argument error:"), "argument failure was not sent to stderr")
    checked += 1

    code, stdout, stderr = run(executable, ["--bach"])
    clean_stdout(stdout)
    require(code == 2 and stdout == "", "removed --bach alias must be rejected")
    require(stderr.startswith("Argument error:"), "removed --bach alias must report an argument error")
    checked += 1

    code, stdout, stderr = run(executable, ["--layout", "multi", "--batch"])
    clean_stdout(stdout)
    require(code == 2 and stdout == "", "--layout must not alter automation output modes")
    require(
        stderr.startswith("Argument error:") and "only available in interactive mode" in stderr,
        "--layout automation rejection must be explicit",
    )
    checked += 1

    payload = "x:=2\nx^3\n1+\n4\n"
    code, stdout, stderr = run(executable, ["--batch"], payload)
    clean_stdout(stdout)
    require(code == 3, "batch must retain the highest syntax exit code")
    require(stdout == "2\n8\n4\n", "batch did not preserve state or continue after an error")
    require(stderr.startswith("SyntaxError:"), "batch syntax failure was not sent to stderr")
    checked += 1

    hostile = "-" * 50000 + "1\n"
    code, stdout, stderr = run(executable, ["--batch"], hostile)
    clean_stdout(stdout)
    require(code == 3 and stdout == "", "hostile parser input did not fail safely")
    require(stderr.startswith("ResourceLimitError:"), "resource limit is not externally identifiable")
    checked += 1

    code, stdout, stderr = run(executable, ["--batch"], "1+\nlog[0]\n5\n")
    clean_stdout(stdout)
    require(code == 4 and stdout == "5\n", "batch did not aggregate exit severity or continue")
    require("SyntaxError:" in stderr and "Error:" in stderr, "batch lost one of multiple diagnostics")
    checked += 1

    code, stdout, stderr = run(
        executable, ["--batch"], ":help sin\n2+3\nIn[1]\n"
    )
    clean_stdout(stdout)
    require(code == 0 and stderr == "", ":help sin must be a successful REPL command")
    require(
        stdout
        == "sin\n"
        "Computes sine with exact special-angle simplification where possible.\n"
        "Usage:\n"
        "  sin[x]\n"
        "Arguments: 1\n"
        "Inputs:\n"
        "  x: A real or complex angle. A bare real uses angleMode[]; Deg, Rad, or Grad overrides it.\n"
        "Notes:\n"
        "  The default session angle mode is Rad.\n"
        "Examples:\n"
        "  sin[Pi/6]  ->  1/2\n"
        "  sin[30 Deg]  ->  1/2\n"
        "  sin[100 Grad]  ->  1\n"
        "5\n"
        "5\n",
        ":help sin output changed or it consumed In[1]",
    )
    checked += 1

    code, stdout, stderr = run(executable, ["--batch"], ":quit\n2+3\n")
    clean_stdout(stdout)
    require((code, stdout, stderr) == (0, "", ""), ":quit must stop a batch session cleanly")
    checked += 1

    code, stdout, stderr = run(executable, ["--batch"], ":exit\n2+3\n")
    clean_stdout(stdout)
    require((code, stdout, stderr) == (0, "", ""), ":exit must alias :quit")
    checked += 1

    code, stdout, stderr = run(executable, ["--batch"], ":help ln\n:help noSuchFunction\n")
    clean_stdout(stdout)
    require(code == 0 and stderr == "", "help lookup must not be an evaluation failure")
    require("log (alias: ln)\n" in stdout, ":help did not resolve an alias")
    require(
        "No help for 'noSuchFunction'.\n"
        "Use :help functions or :help constants to list available topics.\n"
        in stdout,
        "unknown help lookup lost its guidance",
    )
    checked += 1

    code, stdout, stderr = run(
        executable, ["--batch"], ":help Pi\n:help sdv\n:help qr\n2+3\nIn[1]\n"
    )
    clean_stdout(stdout)
    require(code == 0 and stderr == "", "constant and suggested help must succeed")
    require("Pi\nThe exact circle constant" in stdout, ":help Pi is missing")
    require("Did you mean 'svd'?" in stdout, "sdv typo suggestion is missing")
    require(
        "Did you mean 'qrDecomposition'?" in stdout,
        "qr decomposition suggestion is missing",
    )
    require(stdout.endswith("5\n5\n"), "help lookup consumed an In[n] history slot")
    checked += 1

    code, stdout, stderr = run_interrupted(executable, "inverse[identity[500]]")
    clean_stdout(stdout)
    require(code == 3 and stdout == "", "Ctrl-C cancellation must use resource exit code 3")
    require(
        "ResourceLimitError: Evaluation cancelled by frontend" in stderr,
        "Ctrl-C was not routed through EvaluationCancellationToken",
    )
    checked += 1

    print("CLI automation: PASS ({} contracts)".format(checked))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (AssertionError, subprocess.TimeoutExpired) as error:
        print("CLI automation: FAIL: {}".format(error), file=sys.stderr)
        sys.exit(1)
