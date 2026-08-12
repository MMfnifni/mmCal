# -*- coding: utf-8 -*-
"""
選択したフォルダ以下を再帰走査し、対象ファイルの
  - 行数
  - 文字数（改行コードを除く）
  - コメント行数
を集計して CSV に保存する。

Python 3.7.9 対応。
対象拡張子: .txt, .md, .py, .cpp, .hpp
"""

import csv
import io
import os
import sys
import tokenize
import tkinter as tk
from tkinter import filedialog, messagebox


TARGET_EXTENSIONS = {".txt", ".md", ".py", ".cpp", ".hpp"}
OUTPUT_CSV_NAME = "source_stats.csv"


def read_text_file(path, ext):
    """ファイルを文字列として読む。戻り値は (text, encoding)。"""
    with open(path, "rb") as f:
        data = f.read()

    # Python は PEP 263 の encoding 宣言を優先する。
    if ext == ".py":
        try:
            encoding, _ = tokenize.detect_encoding(io.BytesIO(data).readline)
            return data.decode(encoding), encoding
        except (SyntaxError, UnicodeDecodeError, LookupError):
            pass

    # 日本語環境でありがちな文字コードを順に試す。
    for encoding in ("utf-8-sig", "utf-8", "cp932", "shift_jis"):
        try:
            return data.decode(encoding), encoding
        except UnicodeDecodeError:
            continue

    # 最後の保険。処理自体は止めず、置換文字を使う。
    return data.decode("utf-8", errors="replace"), "utf-8(replace)"


def count_python_comment_lines(text):
    """Python の COMMENT トークンが存在する物理行数を数える。"""
    comment_lines = set()
    try:
        tokens = tokenize.generate_tokens(io.StringIO(text).readline)
        for tok in tokens:
            if tok.type == tokenize.COMMENT:
                comment_lines.add(tok.start[0])
    except (tokenize.TokenError, IndentationError, SyntaxError):
        # 壊れた/編集中の Python ファイルでも、それまで読めたコメントは残す。
        pass
    return len(comment_lines)


def count_cpp_comment_lines(text):
    """
    C/C++ の // と /* ... */ を走査し、コメントを含む物理行数を数える。
    文字列・文字リテラル内の // や /* は無視する。
    """
    lines = text.splitlines()
    count = 0
    in_block_comment = False

    for line in lines:
        i = 0
        n = len(line)
        in_string = False
        in_char = False
        escaped = False
        has_comment = in_block_comment

        while i < n:
            if in_block_comment:
                end = line.find("*/", i)
                if end == -1:
                    # この行は最後までブロックコメント。
                    i = n
                    break
                in_block_comment = False
                has_comment = True
                i = end + 2
                continue

            ch = line[i]
            nxt = line[i + 1] if i + 1 < n else ""

            if in_string:
                if escaped:
                    escaped = False
                elif ch == "\\":
                    escaped = True
                elif ch == '"':
                    in_string = False
                i += 1
                continue

            if in_char:
                if escaped:
                    escaped = False
                elif ch == "\\":
                    escaped = True
                elif ch == "'":
                    in_char = False
                i += 1
                continue

            if ch == '"':
                in_string = True
                i += 1
                continue

            if ch == "'":
                in_char = True
                i += 1
                continue

            if ch == "/" and nxt == "/":
                has_comment = True
                break

            if ch == "/" and nxt == "*":
                has_comment = True
                in_block_comment = True
                i += 2
                continue

            i += 1

        if has_comment:
            count += 1

    return count


def count_stats(path):
    ext = os.path.splitext(path)[1].lower()
    text, encoding = read_text_file(path, ext)

    line_count = len(text.splitlines())
    char_count = len(text.replace("\r", "").replace("\n", ""))

    if ext == ".py":
        comment_count = count_python_comment_lines(text)
    elif ext in (".cpp", ".hpp"):
        comment_count = count_cpp_comment_lines(text)
    else:
        # .txt / .md にはプログラム言語としてのコメント構文を定義しない。
        comment_count = 0

    return line_count, char_count, comment_count, encoding


def scan_folder(root_folder):
    rows = []

    for current_dir, dirnames, filenames in os.walk(root_folder):
        # よくある不要ディレクトリは対象外にする。
        dirnames[:] = [
            d for d in dirnames
            if d not in {".git", ".svn", "__pycache__", ".venv", "venv"}
        ]

        for filename in filenames:
            ext = os.path.splitext(filename)[1].lower()
            if ext not in TARGET_EXTENSIONS:
                continue

            path = os.path.join(current_dir, filename)
            rel_path = os.path.relpath(path, root_folder)

            try:
                line_count, char_count, comment_count, encoding = count_stats(path)
                rows.append({
                    "relative_path": rel_path,
                    "extension": ext,
                    "lines": line_count,
                    "characters": char_count,
                    "comment_lines": comment_count,
                    "encoding": encoding,
                    "error": "",
                })
            except Exception as e:
                rows.append({
                    "relative_path": rel_path,
                    "extension": ext,
                    "lines": "",
                    "characters": "",
                    "comment_lines": "",
                    "encoding": "",
                    "error": "{}: {}".format(type(e).__name__, e),
                })

    rows.sort(key=lambda r: r["relative_path"].lower())
    return rows


def save_csv(root_folder, rows):
    output_path = os.path.join(root_folder, OUTPUT_CSV_NAME)

    with open(output_path, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "relative_path",
                "extension",
                "lines",
                "characters",
                "comment_lines",
                "encoding",
                "error",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)

    return output_path


def main():
    root = tk.Tk()
    root.withdraw()
    root.update()

    folder = filedialog.askdirectory(title="集計するフォルダを選択してください")
    if not folder:
        root.destroy()
        return

    rows = scan_folder(folder)
    output_path = save_csv(folder, rows)

    ok_rows = [r for r in rows if not r["error"]]
    error_rows = [r for r in rows if r["error"]]

    total_lines = sum(r["lines"] for r in ok_rows)
    total_chars = sum(r["characters"] for r in ok_rows)
    total_comments = sum(r["comment_lines"] for r in ok_rows)

    msg = (
        "集計が完了しました。\n\n"
        "対象ファイル数: {files}\n"
        "総行数: {lines:,}\n"
        "総文字数: {chars:,}\n"
        "総コメント行数: {comments:,}\n"
        "エラー: {errors}\n\n"
        "CSV:\n{csv_path}"
    ).format(
        files=len(ok_rows),
        lines=total_lines,
        chars=total_chars,
        comments=total_comments,
        errors=len(error_rows),
        csv_path=output_path,
    )

    messagebox.showinfo("ソース集計", msg)
    root.destroy()


if __name__ == "__main__":
    main()
