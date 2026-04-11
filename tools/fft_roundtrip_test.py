from __future__ import annotations

import argparse
import math
import random
import statistics
import subprocess
import sys
import tempfile
import time
from pathlib import Path

from matrix_validator import compare_arrays


def write_real_random_matrix(filename: Path, rows: int, cols: int, seed: int = 1234) -> None:
    rng = random.Random(seed)
    with filename.open("w", encoding="utf-8") as f:
        f.write("{")
        for i in range(rows):
            f.write("{")
            for j in range(cols):
                val = rng.uniform(-1.0, 1.0)
                f.write(f"{val:.10f}")
                if j != cols - 1:
                    f.write(",")
            f.write("}")
            if i != rows - 1:
                f.write(",")
        f.write("}")


def write_sine_matrix(filename: Path, rows: int, cols: int, freq_x: int = 3, freq_y: int = 5) -> None:
    with filename.open("w", encoding="utf-8") as f:
        f.write("{")
        for i in range(rows):
            f.write("{")
            for j in range(cols):
                x = math.sin(2.0 * math.pi * freq_x * i / rows)
                y = math.sin(2.0 * math.pi * freq_y * j / cols)
                val = x + y
                f.write(f"{val:.10f}")
                if j != cols - 1:
                    f.write(",")
            f.write("}")
            if i != rows - 1:
                f.write(",")
        f.write("}")


def write_impulse_matrix(
    filename: Path,
    rows: int,
    cols: int,
    impulse_row: int | None = None,
    impulse_col: int | None = None,
    amplitude: float = 1.0,
) -> None:
    if impulse_row is None:
        impulse_row = rows // 2
    if impulse_col is None:
        impulse_col = cols // 2

    with filename.open("w", encoding="utf-8") as f:
        f.write("{")
        for i in range(rows):
            f.write("{")
            for j in range(cols):
                val = amplitude if (i == impulse_row and j == impulse_col) else 0.0
                f.write(f"{val:.10f}")
                if j != cols - 1:
                    f.write(",")
            f.write("}")
            if i != rows - 1:
                f.write(",")
        f.write("}")

def write_complex_random_matrix(filename: Path, rows: int, cols: int, seed: int = 1234) -> None:
    rng = random.Random(seed)

    with filename.open("w", encoding="utf-8") as f:
        f.write("{")

        for i in range(rows):
            f.write("{")
            for j in range(cols):
                re = rng.uniform(-1.0, 1.0)
                im = rng.uniform(-1.0, 1.0)

                if abs(im) < 1e-15:
                    s = f"{re:.10f}"
                elif im >= 0:
                    s = f"{re:.10f}+{im:.10f}I"
                else:
                    s = f"{re:.10f}{im:.10f}I"

                f.write(s)

                if j != cols - 1:
                    f.write(",")
            f.write("}")
            if i != rows - 1:
                f.write(",")

        f.write("}")

def make_mmcalc_batch_script(input_path: Path, output_path: Path) -> str:
    in_s = str(input_path).replace("\\", "\\\\")
    out_s = str(output_path).replace("\\", "\\\\")
    return "\n".join([
        f'fileout(ifft(fft(filein("{in_s}"))), "{out_s}")',
        "Exit[]",
        "",
    ])


def run_mmcalc(exe: Path, script_text: str, timeout_sec: int = 3600) -> tuple[subprocess.CompletedProcess[str], float]:
    start = time.perf_counter()
    proc = subprocess.run(
        [str(exe), "--batch"],
        input=script_text,
        text=True,
        capture_output=True,
        timeout=timeout_sec,
        check=False,
    )
    elapsed = time.perf_counter() - start
    return proc, elapsed


def judge(
    result: dict[str, float],
    max_abs: float,
    mean_abs: float,
    max_rel: float,
    energy_rel: float,
) -> tuple[bool, list[str]]:
    reasons: list[str] = []

    if result["max_err"] > max_abs:
        reasons.append(f"max_err {result['max_err']} > {max_abs}")
    if result["mean_err"] > mean_abs:
        reasons.append(f"mean_err {result['mean_err']} > {mean_abs}")
    if result["max_rel_err"] > max_rel:
        reasons.append(f"max_rel_err {result['max_rel_err']} > {max_rel}")
    if result["energy_rel_err"] > energy_rel:
        reasons.append(f"energy_rel_err {result['energy_rel_err']} > {energy_rel}")

    return (len(reasons) == 0, reasons)


def print_timing_summary(label: str, values: list[float]) -> None:
    if not values:
        return
    print(f"{label}:")
    print(f"  min : {min(values):.6f} s")
    print(f"  mean: {statistics.mean(values):.6f} s")
    print(f"  max : {max(values):.6f} s")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--exe", default="../build/x64/Release/mmCal_x64.exe",
                        help="Path to mmCalculator executable")
    parser.add_argument("--size", type=int, default=512)
    parser.add_argument("--kind", choices=["random", "sin", "impulse", "complex_random"], default="random")
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument("--max-abs", type=float, default=1e-11)
    parser.add_argument("--mean-abs", type=float, default=1e-13)
    parser.add_argument("--max-rel", type=float, default=1e-7)
    parser.add_argument("--energy-rel", type=float, default=1e-15)
    parser.add_argument("--repeat", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--keep-files", action="store_true")
    args = parser.parse_args()

    exe = Path(args.exe).resolve()
    if not exe.exists():
        print(f"ERROR: executable not found: {exe}", file=sys.stderr)
        return 2

    total_start = time.perf_counter()

    if args.keep_files:
        tmp = Path.cwd() / f"fft_rt_{args.kind}_{args.size}"
        tmp.mkdir(parents=True, exist_ok=True)
        temp_ctx = None
    else:
        temp_ctx = tempfile.TemporaryDirectory(prefix="mmcalc_fft_")
        tmp = Path(temp_ctx.name)

    try:
        input_path = tmp / f"fft_{args.kind}_{args.size}.txt"
        output_path = tmp / f"out_{args.kind}_{args.size}.txt"

        gen_start = time.perf_counter()
        if args.kind == "random":
            write_real_random_matrix(input_path, args.size, args.size, seed=args.seed)
        elif args.kind == "sin":
            write_sine_matrix(input_path, args.size, args.size)
        elif args.kind == "impulse":
            write_impulse_matrix(input_path, args.size, args.size)
        else:
            write_complex_random_matrix(input_path, args.size, args.size)
        gen_elapsed = time.perf_counter() - gen_start

        print("===================================")
        print("FFT Round-trip Test - mmKreutzef")
        print("===================================")
        print(f"Executable : {exe}")
        print(f"Input kind : {args.kind}")
        print(f"Matrix size: {args.size} x {args.size}")
        print(f"Input file : {input_path}")
        print(f"Output file: {output_path}")
        print(f"Generate time: {gen_elapsed:.6f} s")
        print("-----------------------------------")

        script = make_mmcalc_batch_script(input_path, output_path)

        mmcalc_times: list[float] = []
        for i in range(args.repeat):
            if output_path.exists():
                output_path.unlink()

            proc, elapsed = run_mmcalc(exe, script, timeout_sec=args.timeout)
            mmcalc_times.append(elapsed)

            print(f"Run {i + 1}/{args.repeat}: mmCalc fft->ifft elapsed = {elapsed:.6f} s")

            if proc.returncode != 0:
                print("mmCalculator process failed", file=sys.stderr)
                print("--- stdout ---")
                print(proc.stdout)
                print("--- stderr ---", file=sys.stderr)
                print(proc.stderr, file=sys.stderr)
                return 3

            if not output_path.exists():
                print("ERROR: output file was not created", file=sys.stderr)
                print("--- stdout ---")
                print(proc.stdout)
                print("--- stderr ---", file=sys.stderr)
                print(proc.stderr, file=sys.stderr)
                return 4

        print("-----------------------------------")
        print_timing_summary("mmCalc execution time", mmcalc_times)
        print("-----------------------------------")

        compare_start = time.perf_counter()
        result = compare_arrays(str(input_path), str(output_path))
        compare_elapsed = time.perf_counter() - compare_start

        ok, reasons = judge(
            result,
            max_abs=args.max_abs,
            mean_abs=args.mean_abs,
            max_rel=args.max_rel,
            energy_rel=args.energy_rel,
        )

        print("-----------------------------------")
        print(f"Validator elapsed: {compare_elapsed:.6f} s")
        total_elapsed = time.perf_counter() - total_start
        print(f"Total harness time: {total_elapsed:.6f} s")
        print("===================================")

        if ok:
            print(f"[PASS] FFT round-trip ({args.kind}, {args.size}x{args.size})")
            return 0

        print(f"[FAIL] FFT round-trip ({args.kind}, {args.size}x{args.size})")
        for reason in reasons:
            print(" -", reason)
        return 1

    finally:
        if temp_ctx is not None:
            temp_ctx.cleanup()


if __name__ == "__main__":
    raise SystemExit(main())