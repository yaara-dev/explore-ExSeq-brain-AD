#!/usr/bin/env python3
"""Convert README.md to HTML suitable for Zenodo's rich-text description field."""

from __future__ import annotations

import html as H
import re
import subprocess
import sys
from pathlib import Path


def inline(text: str) -> str:
    codes: list[str] = []

    def save_code(m: re.Match[str]) -> str:
        codes.append(m.group(1))
        return f"\x00CODE{len(codes) - 1}\x00"

    text = re.sub(r"`([^`]+)`", save_code, text)
    text = H.escape(text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2">\1</a>', text)
    for i, c in enumerate(codes):
        text = text.replace(f"\x00CODE{i}\x00", f"<code>{H.escape(c)}</code>")
    return text


def md_to_html(md: str) -> str:
    lines = md.splitlines()
    out: list[str] = []
    i = 0
    in_code = False
    code_buf: list[str] = []
    ul_open = False
    ol_open = False
    # Allow 1-column tables like |------|
    sep_re = re.compile(r"^\s*\|?\s*:?-+:?\s*(\|\s*:?-+:?\s*)*\|?\s*$")

    def close_lists() -> None:
        nonlocal ul_open, ol_open
        if ul_open:
            out.append("</ul>")
            ul_open = False
        if ol_open:
            out.append("</ol>")
            ol_open = False

    while i < len(lines):
        line = lines[i]
        if line.startswith("```"):
            close_lists()
            if not in_code:
                in_code = True
                code_buf = []
            else:
                out.append("<pre><code>" + H.escape("\n".join(code_buf)) + "</code></pre>")
                in_code = False
            i += 1
            continue
        if in_code:
            code_buf.append(line)
            i += 1
            continue

        if (
            line.strip().startswith("|")
            and i + 1 < len(lines)
            and sep_re.match(lines[i + 1])
        ):
            close_lists()
            rows: list[str] = []
            while i < len(lines) and "|" in lines[i]:
                rows.append(lines[i])
                i += 1
            out.append("<table>")
            for ridx, row in enumerate(rows):
                if ridx == 1:
                    continue
                cells = [c.strip() for c in row.strip().strip("|").split("|")]
                tag = "th" if ridx == 0 else "td"
                cells_html = "".join(f"<{tag}>{inline(c)}</{tag}>" for c in cells)
                out.append(f"<tr>{cells_html}</tr>")
            out.append("</table>")
            continue

        if line.strip() == "---":
            close_lists()
            out.append("<hr>")
            i += 1
            continue

        m = re.match(r"^(#{1,6})\s+(.*)$", line)
        if m:
            close_lists()
            level = len(m.group(1))
            out.append(f"<h{level}>{inline(m.group(2))}</h{level}>")
            i += 1
            continue

        m = re.match(r"^[-*]\s+(.*)$", line)
        if m:
            if ol_open:
                out.append("</ol>")
                ol_open = False
            if not ul_open:
                out.append("<ul>")
                ul_open = True
            out.append(f"<li>{inline(m.group(1))}</li>")
            i += 1
            continue

        m = re.match(r"^\d+\.\s+(.*)$", line)
        if m:
            if ul_open:
                out.append("</ul>")
                ul_open = False
            if not ol_open:
                out.append("<ol>")
                ol_open = True
            out.append(f"<li>{inline(m.group(1))}</li>")
            i += 1
            continue

        if not line.strip():
            close_lists()
            i += 1
            continue

        close_lists()
        paras = [line]
        i += 1
        while i < len(lines):
            nxt = lines[i]
            if not nxt.strip():
                break
            if nxt.startswith("#") or nxt.startswith("```") or nxt.strip() == "---":
                break
            if re.match(r"^[-*]\s+", nxt) or re.match(r"^\d+\.\s+", nxt):
                break
            if nxt.strip().startswith("|") and i + 1 < len(lines) and sep_re.match(lines[i + 1]):
                break
            paras.append(nxt)
            i += 1
        out.append("<p>" + "<br>".join(inline(p) for p in paras) + "</p>")

    close_lists()
    if in_code:
        out.append("<pre><code>" + H.escape("\n".join(code_buf)) + "</code></pre>")
    return "\n".join(out)


GITHUB_BLOB = "https://github.com/yaara-dev/early-AD-data-code/blob/main/"
GITHUB_TREE = "https://github.com/yaara-dev/early-AD-data-code/tree/main/"


def absolutize_repo_links(html: str) -> str:
    """Point relative README links at the public GitHub repo."""

    def repl(m: re.Match[str]) -> str:
        href = m.group(1)
        if href.startswith(("http://", "https://", "mailto:", "#")):
            return m.group(0)
        base = GITHUB_TREE if href.endswith("/") else GITHUB_BLOB
        return f'href="{base}{href}"'

    return re.sub(r'href="([^"]+)"', repl, html)


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    src = root / "README.md"
    dst = root / "zenodo_description.html"
    html = absolutize_repo_links(md_to_html(src.read_text()))
    dst.write_text(html)
    try:
        subprocess.run(["pbcopy"], input=html.encode(), check=True)
        copied = True
    except Exception:
        copied = False
    print(f"Wrote {dst} ({len(html)} chars)")
    if copied:
        print("Copied HTML to clipboard — paste into Zenodo Description (rich text / HTML mode).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
