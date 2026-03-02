#!/usr/bin/env python3
"""
Convert manual.tex to GitHub-friendly Markdown files (one per chapter) with
internal links. Requires pandoc. Run from the manual/ directory:

    python3 tex2md.py

Produces:
  - md/01-general-remarks.md, md/02-global-optimization-clusters.md, ...
  - md/README.md (index with links to all chapters)

GitHub renders these natively; links work across files.
"""

import re
import subprocess
import sys
from pathlib import Path
from typing import Optional

MANUAL_DIR = Path(__file__).resolve().parent
MANUAL_TEX = MANUAL_DIR / "manual.tex"
MANUAL_BIB = MANUAL_DIR / "manual.bib"
MD_DIR = MANUAL_DIR / "md"
FULL_MD = MANUAL_DIR / "_manual_full.md"


def strip_html(text: str) -> str:
    """Remove HTML tags for plain-text slug/title."""
    return re.sub(r"<[^>]+>", "", text)


def title_to_slug(title: str) -> str:
    """e.g. 'General remarks' -> 'general-remarks'. Safe for filenames."""
    t = strip_html(title)
    s = re.sub(r"[^\w\s-]", "", t.lower())
    return re.sub(r"[-\s]+", "-", s).strip("-") or "chapter"


def slug_from_index_and_title(index: int, title: str) -> str:
    """e.g. 1, 'General remarks' -> '01-general-remarks'."""
    return f"{index:02d}-{title_to_slug(title)}"


def run_pandoc() -> bool:
    """Convert manual.tex to a single Markdown file. Returns True on success."""
    args = [
        "pandoc",
        str(MANUAL_TEX),
        "-f", "latex",
        "-t", "gfm",
        "--number-sections",
        "--wrap=auto",
        "-o", str(FULL_MD),
    ]
    if MANUAL_BIB.exists():
        args.extend(["--bibliography", str(MANUAL_BIB), "--citeproc"])
    try:
        subprocess.run(
            args,
            check=True,
            capture_output=True,
            cwd=str(MANUAL_DIR),
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        if isinstance(e, FileNotFoundError):
            print("error: pandoc not found. Install pandoc to convert LaTeX to Markdown.", file=sys.stderr)
        else:
            print("error: pandoc failed:", e.stderr.decode() if e.stderr else e, file=sys.stderr)
        return False
    return True


def split_and_collect_chapters() -> list[tuple[str, str, str]]:
    """
    Split _manual_full.md by top-level # (chapter) headings.
    Pandoc outputs "# Title" without section numbers. We number chapters by order.
    Returns list of (slug, chapter_num_str, markdown_content).
    """
    if not FULL_MD.exists():
        return []

    text = FULL_MD.read_text(encoding="utf-8", errors="replace")
    # Split on any level-1 heading: "# Title" (pandoc does not add numbers in GFM)
    pattern = re.compile(r"^# (.+)$", re.MULTILINE)
    parts = pattern.split(text)

    # parts[0] = content before first chapter (often empty or title/TOC)
    # parts[1], parts[2] = "General remarks", "\n\n## This manual\n..."
    # parts[3], parts[4] = "Global optimization...", "\n\n..."
    chapters = []
    i = 1
    chapter_num = 1
    while i + 1 < len(parts):
        title, content = parts[i], parts[i + 1]
        slug = slug_from_index_and_title(chapter_num, title)
        full_content = f"# {title}\n{content}"
        chapters.append((slug, str(chapter_num), full_content))
        chapter_num += 1
        i += 2

    return chapters


def normalize_paragraph_breaks(content: str) -> str:
    """
    Follow LaTeX convention: only double newline means a real line/paragraph break.
    Within a paragraph, collapse single newlines to spaces so prose flows as one block.
    Leave list items, code blocks, and tables unchanged.
    """
    blocks = re.split(r"\n\n+", content)
    result = []
    for block in blocks:
        block = block.rstrip("\n")
        if not block:
            result.append("")
            continue
        lines = block.split("\n")
        first_line = lines[0]
        # Leave lists, indented code (4+ spaces), fenced code, tables as-is
        if re.match(r"^[\s]*[-*+]\s", first_line):
            result.append(block)
            continue
        if re.match(r"^[\s]{4,}", first_line):
            result.append(block)
            continue
        if first_line.strip().startswith("```"):
            result.append(block)
            continue
        if first_line.strip().startswith("|"):
            result.append(block)
            continue
        if re.match(r"^[\s]*\d+\.\s", first_line):
            result.append(block)
            continue
        # Heading on first line: keep heading, collapse the rest as one paragraph
        if re.match(r"^#{1,6}\s", first_line):
            rest = " ".join(lines[1:]) if len(lines) > 1 else ""
            result.append(f"{first_line}\n\n{rest}" if rest else first_line)
            continue
        # Prose paragraph: join all lines with a single space
        result.append(" ".join(lines))
    return "\n\n".join(result)


REFERENCE_SLUG = "10-references"


def extract_references_block(content: str) -> tuple[str, Optional[str]]:
    """
    If content ends with a pandoc-citeproc references block (div id="refs" or
    class="references csl-bib-body"), extract it and return (content_without_refs,
    refs_as_markdown). Otherwise return (content, None).
    """
    refs_start = re.search(
        r'<div\s+id="refs"\s+class="references[^"]*"[^>]*>',
        content,
        re.IGNORECASE,
    )
    if not refs_start:
        refs_start = re.search(r'<div[^>]*class="references\s+csl-bib-body"[^>]*>', content)
    if not refs_start:
        return content, None

    start_pos = refs_start.start()
    # Find matching closing </div> by counting nested divs
    pos = refs_start.end()
    depth = 1
    while depth > 0 and pos < len(content):
        next_open = content.find("<div", pos)
        next_close = content.find("</div>", pos)
        if next_close == -1:
            break
        if next_open != -1 and next_open < next_close:
            depth += 1
            pos = next_open + 4
        else:
            depth -= 1
            pos = next_close + 6
            if depth == 0:
                refs_html = content[start_pos:next_close + 6]
                break
    else:
        return content, None

    content_without = content[:start_pos].rstrip()
    refs_md = references_html_to_markdown(refs_html)
    return content_without, refs_md


def references_html_to_markdown(refs_div_html: str) -> str:
    """Convert pandoc refs div to a Markdown section: # References and a list of entries."""
    entries = re.findall(
        r'<div[^>]*class="csl-entry"[^>]*>\s*(.*?)\s*</div>',
        refs_div_html,
        re.DOTALL | re.IGNORECASE,
    )
    if not entries:
        return ""
    lines = ["# References", ""]
    for entry in entries:
        # Normalize whitespace (single newlines to space, keep one line per entry)
        text = " ".join(entry.split())
        # Keep * for italics and <url> for links (valid in Markdown)
        lines.append(f"- {text}")
    return "\n".join(lines)


def rewrite_cross_chapter_links(
    content: str, this_anchor: str, slug_for_anchor: dict[str, str]
) -> str:
    """
    Rewrite [text](#anchor) that point to another chapter to [text](slug.md#anchor).
    slug_for_anchor: map of GFM anchor -> filename slug for each chapter.
    """
    def repl(m: re.Match) -> str:
        anchor = m.group(1).lstrip("#")
        if anchor == this_anchor:
            return m.group(0)
        slug = slug_for_anchor.get(anchor)
        if slug:
            return f"]({slug}.md#{anchor})"
        return m.group(0)

    content = re.sub(r"\]\((#[-\w]+)\)", repl, content)
    return content


def write_chapters(chapters: list[tuple[str, str, str]]) -> tuple[list[tuple[str, str]], bool]:
    """Write each chapter to md/<slug>.md. Returns (chapter_slugs, had_references)."""
    MD_DIR.mkdir(parents=True, exist_ok=True)
    slug_for_anchor = {}
    for slug, _num, content in chapters:
        first_line = content.split("\n")[0].strip().lstrip("#").strip()
        anchor = title_to_slug(strip_html(first_line))
        slug_for_anchor[anchor] = slug

    had_references = False
    for idx, (slug, _num, content) in enumerate(chapters):
        first_line = content.split("\n")[0].strip().lstrip("#").strip()
        this_anchor = title_to_slug(strip_html(first_line))
        content = rewrite_cross_chapter_links(content, this_anchor, slug_for_anchor)
        content = normalize_paragraph_breaks(content)

        # Last chapter: extract references into their own section if present
        if idx == len(chapters) - 1:
            content, refs_md = extract_references_block(content)
            if refs_md:
                (MD_DIR / f"{REFERENCE_SLUG}.md").write_text(refs_md, encoding="utf-8")
                had_references = True

        (MD_DIR / f"{slug}.md").write_text(content, encoding="utf-8")

    return [(slug, num) for slug, num, _ in chapters], had_references


def write_index(
    chapter_slugs: list[tuple[str, str]],
    chapters: list[tuple[str, str, str]],
    has_references: bool = False,
) -> None:
    """Write md/README.md as the manual index."""
    title = "OGOLEM Manual"
    lines = [
        f"# {title}",
        "",
        "This manual is generated from the LaTeX source (`manual.tex`). "
        "Below: links to each chapter (viewable on GitHub).",
        "",
        "## Chapters",
        "",
    ]
    for (slug, _num), (_, _, content) in zip(chapter_slugs, chapters):
        first_line = content.split("\n")[0].strip()
        title_text = strip_html(first_line.lstrip("#").strip())
        lines.append(f"- [{title_text}]({slug}.md)")
    if has_references:
        lines.append(f"- [References]({REFERENCE_SLUG}.md)")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*To rebuild from LaTeX, run from the `manual/` directory:* `python3 tex2md.py`")
    (MD_DIR / "README.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    if not run_pandoc():
        return 1

    chapters = split_and_collect_chapters()
    if not chapters:
        print("error: no chapters found in pandoc output.", file=sys.stderr)
        return 1

    chapter_slugs, has_references = write_chapters(chapters)
    write_index(chapter_slugs, chapters, has_references=has_references)

    try:
        FULL_MD.unlink()
    except OSError:
        pass

    n = len(chapters) + (1 if has_references else 0)
    print(f"Wrote {n} section(s) (including references) to {MD_DIR}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())
