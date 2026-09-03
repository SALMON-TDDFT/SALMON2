#!/usr/bin/env python3
"""Static safety contract for the Task 5 generalized Hybrid pipeline.

The numerical MPI fixtures intentionally use small matrices and grids.  This
check protects the production source against count overflow and against
per-element/per-row collective patterns that those fixtures cannot exercise at
production scale.
"""

from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[2]
PIPELINE_PATH = ROOT / "src/gs/dc/dg_hybrid_projected_fragment_pipeline.f90"
SOURCE = PIPELINE_PATH.read_text(errors="replace")


def without_comments(text: str) -> str:
    """Remove Fortran comments while preserving statement and loop boundaries."""

    cleaned = []
    for line in text.splitlines():
        # The reviewed procedures contain no quoted exclamation marks.  Keeping
        # this intentionally small makes the source contract easy to audit.
        cleaned.append(line.split("!", 1)[0])
    return "\n".join(cleaned)


def extract_procedure(source: str, name: str, kind: str) -> str:
    """Return one named free-form Fortran procedure, excluding interfaces."""

    start_pattern = re.compile(
        rf"^[ \t]*(?!end\b)(?:[a-z][a-z0-9_]*(?:\([^\n]*?\))?[ \t]+)*"
        rf"{kind}[ \t]+{re.escape(name)}\b",
        re.IGNORECASE | re.MULTILINE,
    )
    start = start_pattern.search(source)
    if start is None:
        raise RuntimeError(f"missing {kind} {name} in {PIPELINE_PATH}")
    end_pattern = re.compile(
        rf"^[ \t]*end[ \t]+{kind}(?:[ \t]+{re.escape(name)})?[ \t]*$",
        re.IGNORECASE | re.MULTILINE,
    )
    end = end_pattern.search(source, start.end())
    if end is None:
        raise RuntimeError(f"unterminated {kind} {name} in {PIPELINE_PATH}")
    return source[start.start() : end.end()]


def logical_lines(text: str) -> list[str]:
    """Normalize continuations/semicolons without erasing loop nesting."""

    text = without_comments(text).lower()
    text = re.sub(r"&[ \t]*\n[ \t]*", "", text)
    text = text.replace(";", "\n")
    return [re.sub(r"\s+", " ", line).strip() for line in text.splitlines() if line.strip()]


def compact(text: str) -> str:
    return re.sub(r"\s+", "", without_comments(text).lower().replace("&", ""))


def loop_regions(text: str) -> list[tuple[str, list[str]]]:
    """Extract nested DO regions from a procedure using normalized statements."""

    regions: list[tuple[str, list[str]]] = []
    stack: list[tuple[str, int]] = []
    lines = logical_lines(text)
    for index, line in enumerate(lines):
        if re.match(r"^(?:[a-z][a-z0-9_]*\s*:\s*)?do\b", line):
            stack.append((line, index))
        elif re.match(r"^end\s*do\b", line):
            if not stack:
                raise RuntimeError(f"unbalanced END DO while parsing {PIPELINE_PATH}")
            header, begin = stack.pop()
            regions.append((header, lines[begin + 1 : index]))
    if stack:
        raise RuntimeError(f"unbalanced DO while parsing {PIPELINE_PATH}")
    return regions


def assignment_expressions(text: str) -> dict[str, str]:
    assignments: dict[str, str] = {}
    for statement in logical_lines(text):
        match = re.match(r"^([a-z][a-z0-9_]*)\s*=\s*(.*)$", statement)
        if match:
            assignments[match.group(1)] = re.sub(r"\s+", "", match.group(2))
    return assignments


failures: list[str] = []


def require(condition: bool, message: str) -> None:
    if not condition:
        failures.append(message)


build = extract_procedure(SOURCE, "build_dg_hybrid_projected_fragment_basis", "subroutine")
finalize = extract_procedure(SOURCE, "finalize_dg_hybrid_dual_basis_catalog", "subroutine")
coverage = extract_procedure(SOURCE, "validate_coverage_pairs", "subroutine")

build_compact = compact(build)
finalize_compact = compact(finalize)

# Count and receipt arithmetic must not first overflow a default integer.
require(
    re.search(r"\bnpw=npw\+", build_compact) is None,
    "unsafe default-integer accumulation remains in the projected PW count",
)
require(
    re.search(r"int\(size\([^)]*\)\+size\(", build_compact) is None,
    "workspace receipt sums SIZE values before widening them",
)
require(
    "sum(expected_fragment_ranks)" not in finalize_compact,
    "expected fragment ranks use an unchecked default-integer SUM",
)
require(
    re.search(
        r"mpi_allreduce\([^\n]*,nuncompressed\*nuncompressed,", finalize_compact
    )
    is None,
    "dual Gram MPI count multiplies default integers without a checked product",
)
require(
    "tail_tolerance**2" not in finalize_compact,
    "tail acceptance squares an unrestricted finite tolerance",
)
require(
    "ncomplete<=nuncompressed" in finalize_compact
    and "nseed<=ncomplete" in finalize_compact,
    "terminal and seed dimensions are not bounded before square workspaces",
)
require(
    "expected_fragment_wannier_ranks>=0" in finalize_compact
    and "expected_fragment_wannier_ranks>=1" not in finalize_compact,
    "the optional WF-sector cardinality wrongly excludes a PW-only fragment",
)

# Replicated-input validation may use a bounded collective hash/chunk, but it
# must never issue a collective for every scalar element.
for helper_name, helper_kind in (
    ("replicated_integer_array", "function"),
    ("replicated_int64_array", "function"),
    ("replicated_complex_matrix", "function"),
):
    helper = extract_procedure(SOURCE, helper_name, helper_kind)
    for header, body in loop_regions(helper):
        require(
            not any("mpi_allreduce" in statement for statement in body),
            f"{helper_name} performs MPI_Allreduce inside element loop `{header}`",
        )

# The terminal fingerprint and sparse publisher checks must use a bounded
# number of collectives, independent of global rows and canonical basis rank.
for header, body in loop_regions(finalize):
    body_text = "".join(body)
    if "global_row_count" in header:
        require(
            "mpi_bcast" not in body_text,
            "terminal map fingerprint broadcasts once per global spatial row",
        )
    canonical_column_loop = (
        "nuncompressed" in header
        or "size(uncompressed_global_basis_ids)" in header.replace(" ", "")
    )
    if canonical_column_loop:
        require(
            "mpi_bcast" not in body_text,
            "sparse publisher payload broadcasts once per canonical basis column",
        )

# Geometric/operator coverage is established by raw buffer row IDs.  A valid
# row can be a node of every retained basis function, so nonzero amplitudes are
# not an admissible coverage gate.
for line in logical_lines(coverage):
    if "buffer_values" in line:
        require(
            not any(token in line for token in ("abs(", "any(", ">0", "/=0", "nonzero")),
            "coverage requires a nonzero basis amplitude instead of raw row presence",
        )

# A deficient terminal map must span the coefficient-space range of G, not
# merely have orthonormal columns and a nonsingular T^H G T.  Require an
# explicit comparison of T T^H with the Moore--Penrose range projector G G+.
statements = logical_lines(finalize)
assignments = assignment_expressions(finalize)
map_expression = "matmul(union_to_complete,conjg(transpose(union_to_complete)))"
metric_expressions = ("matmul(gram,gram_inverse)", "matmul(gram_inverse,gram)")
map_projectors = {
    name
    for name, expression in assignments.items()
    if map_expression in expression
}
metric_projectors = {
    name
    for name, expression in assignments.items()
    if any(metric_expression in expression for metric_expression in metric_expressions)
}
require(map_expression in finalize_compact, "terminal-map validation does not form T T^H")
require(
    any(metric_expression in finalize_compact for metric_expression in metric_expressions),
    "terminal-map validation does not form the G G+ retained projector",
)

span_comparison = False
comparison_variable = ""
for statement in statements:
    match = re.match(r"^([a-z][a-z0-9_]*)\s*=\s*(.*)$", statement)
    if match is None:
        continue
    name = match.group(1)
    expression = re.sub(r"\s+", "", match.group(2))
    mentions_map = map_expression in expression or any(variable in expression for variable in map_projectors)
    mentions_metric = any(metric_expression in expression for metric_expression in metric_expressions) or any(
        variable in expression for variable in metric_projectors
    )
    if mentions_map and mentions_metric and ("abs(" in expression or "norm" in expression):
        span_comparison = True
        comparison_variable = name
        break
require(span_comparison, "terminal-map validation never compares T T^H with G G+")
if span_comparison:
    guarded = any(
        comparison_variable in statement
        and re.search(r"(?:>|>=|\.gt\.|\.ge\.)", statement)
        for statement in statements
        if statement.startswith("if")
    )
    require(guarded, "terminal retained-span defect is computed but is not an acceptance gate")


if failures:
    print("Task 5 pipeline source-safety contract: RED")
    for failure in failures:
        print(f"- {failure}")
    raise SystemExit(1)

print("Task 5 pipeline source-safety contract: PASS")
