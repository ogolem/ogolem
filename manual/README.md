# OGOLEM Manual

LaTeX source for the OGOLEM manual (input formats, keywords, backends, and usage).

## How to get the manual

- **Pre-built PDF**  
  The manual is built automatically on every push/PR. To get the PDF:
  1. Open the **Actions** tab of this repository (on GitHub).
  2. Open the latest successful workflow run.
  3. Download the **ogolem_artifacts** artifact; it contains `manual.pdf` (along with the jar and other build outputs).

- **Build locally**  
  From the repository root:
  ```bash
  cd manual && pdflatex manual.tex && bibtex manual && pdflatex manual.tex && pdflatex manual.tex
  ```
  Then open or copy `manual/manual.pdf`.

## Markdown version (GitHub-friendly)

A Markdown version of the manual lives in **[md/](md/)**. It is split into one file per chapter with internal links, so you can read and navigate it directly on GitHub (no build step).

- **Index:** [md/README.md](md/README.md) — links to all chapters.
- **Regenerate from LaTeX:** Install [pandoc](https://pandoc.org/), then from this directory run:
  ```bash
  python3 tex2md.py
  ```
  This overwrites the contents of `md/` with a fresh conversion from `manual.tex`. Commit the updated `md/` files to keep the on-GitHub manual in sync with the LaTeX source.

## Contents

The manual covers:

- Building and running OGOLEM, command-line modes, and JRE requirements  
- Global optimization of clusters: input structure, building blocks, constraints, keywords, local optimizers and backends  
- Global parametrization (potentials, pseudopotentials, etc.)  
- Molecular design (switchable molecules, general design)  
- Lionbench meta-benchmarking  
- Lid/threshold capabilities  
- JVM-based distributed parallelization (RMI)  
- FAQs and troubleshooting  

For a quick start, see the main [README](../README.md) in the repository root.
