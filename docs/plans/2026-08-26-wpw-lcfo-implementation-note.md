# WF+PW Divided-SCF/LCFO Technical Note Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Produce a self-contained Japanese LaTeX technical note that derives the mathematics and maps it precisely to the current SALMON implementation and validation state.

**Architecture:** Use one standalone source file under `docs/notes/`, organized mathematics-first with implementation and validation appendices. Treat source, tests, commits, and retained Si64 evidence as primary records; explicitly distinguish completed unit-level work from the still-running physical validation.

**Tech Stack:** LaTeX, `latexmk`/LuaLaTeX when available, Python textual contract checks, Git.

---

### Task 1: Freeze the note contract

**Files:**
- Create: `tests/dg/check_wpw_lcfo_implementation_note.py`
- Create: `docs/notes/wpw_lcfo_divided_scf_implementation.tex`

**Step 1: Write the failing textual contract**

Require the LaTeX source to contain the four design invariants, all principal mathematical objects (`H`, `S`, window functions, density, generalized eigenproblem), MPI physical-point ownership, the Si64 execution conditions, implementation status labels, and source/test path references.

**Step 2: Run the contract and verify RED**

Run: `python3 tests/dg/check_wpw_lcfo_implementation_note.py`

Expected: failure because the LaTeX note does not exist.

**Step 3: Add the minimal LaTeX skeleton**

Create a compilable document with title, abstract, table of contents, section headings, notation table, and explicit status notice.

**Step 4: Run the contract and verify GREEN**

Run: `python3 tests/dg/check_wpw_lcfo_implementation_note.py`

Expected: `WF+PW LCFO implementation note contract: PASS`.

### Task 2: Write the mathematical and algorithmic core

**Files:**
- Modify: `docs/notes/wpw_lcfo_divided_scf_implementation.tex`

**Step 1: Add continuous-to-discrete derivations**

Derive the Kohn--Sham weak form, quadrature-weighted overlap and Hamiltonian matrices, windowed plane waves, fragment bases, density reconstruction, mixing map, and the final generalized eigenproblem. Define every index and normalization.

**Step 2: Add the route state machine**

Prove from the algorithm ordering that divided SCF converges before exactly one LCFO solve, followed only by validation/publication, with no density recomputation or convergence gate.

**Step 3: Run the textual contract**

Run: `python3 tests/dg/check_wpw_lcfo_implementation_note.py`

Expected: PASS.

### Task 3: Add source mapping and diagnostic history

**Files:**
- Modify: `docs/notes/wpw_lcfo_divided_scf_implementation.tex`

**Step 1: Map mathematical objects to implementation files**

Document the relevant `src/common`, `src/gs/dc`, and `src/gs/main_dft.f90` procedures, callback communicator contracts, physical-point redistribution, fingerprints, and checkpoint publication.

**Step 2: Document validation honestly**

List focused MPI tests and commit milestones. Explain each retained Si64 failure (window-input ownership, local/global row action, generated bulk `memset`, and host-frame descriptor pressure), its evidence, and its fix. Mark the newest Si64 run as pending unless it has completed successfully by verification time.

**Step 3: Run the textual contract**

Run: `python3 tests/dg/check_wpw_lcfo_implementation_note.py`

Expected: PASS.

### Task 4: Render and inspect

**Files:**
- Modify if needed: `docs/notes/wpw_lcfo_divided_scf_implementation.tex`

**Step 1: Detect an available LaTeX engine**

Run: `command -v latexmk || command -v lualatex || command -v pdflatex`

Expected: an engine path, or a recorded environment limitation.

**Step 2: Build without modifying tracked source indirectly**

Run `latexmk -lualatex -interaction=nonstopmode -halt-on-error -outdir=<temporary-directory> docs/notes/wpw_lcfo_divided_scf_implementation.tex`, or the closest available equivalent.

Expected: exit 0 and a PDF in the temporary directory.

**Step 3: Audit the log**

Check for undefined control sequences, undefined references, missing glyphs, and material overfull boxes. Patch the source and rebuild until clean or document a toolchain limitation.

**Step 4: Inspect exact staged content before commit**

Run: `git diff --cached --check` and `git diff --cached` after explicitly staging only the note and its contract test.

**Step 5: Commit**

Commit message: `docs: add divided LCFO implementation note`.
