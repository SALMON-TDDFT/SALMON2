#!/usr/bin/env python3
import argparse
import hashlib
import math
from pathlib import Path
import re
import shlex
import tempfile


ROOT = Path(__file__).resolve().parents[2]
FORMAL_SAMPLE = ROOT / "samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac"
FORMAL_INPUT_NAME = "inputfile_gs_w90_pseudo_sawf_aligned_nw576_nb664"
FORMAL_ATOM_NAME = "atom_sawf_aligned.dat"
FORMAL_SYMMETRY_NAME = "sym_sawf_aligned.dat"
FORMAL_SYSNAME = "c64_dc_pseudo_sawf_aligned_nw576_nb664"
FORMAL_ATOM_SHA256 = "d38a4a6fa0e772816610fbe3e7c603ee4a8bc061876883ee9b295e479db71df2"


class ValidationError(RuntimeError):
    pass


FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"


def _number(text: str, label: str) -> float:
    try:
        value = float(text.replace("D", "E").replace("d", "e"))
    except ValueError as exc:
        raise ValidationError(f"invalid {label}: {text}") from exc
    if not math.isfinite(value):
        raise ValidationError(f"non-finite {label}: {text}")
    return value


def _read(path: Path, label: str) -> str:
    try:
        return path.read_text(errors="replace")
    except OSError as exc:
        raise ValidationError(f"cannot read {label} {path}: {exc}") from exc


def _seed_location(log_text: str) -> tuple[Path, str, str]:
    run_lines = [line for line in log_text.splitlines()
                 if "[DC-LCFO-WANNIER] run:" in line]
    if len(run_lines) != 1:
        raise ValidationError("expected exactly one external Wannier90 run line")
    run_line = run_lines[0]
    match = re.search(r"\bcd\s+'([^']+)'\s+&&.*'([^']+)'\s*$", run_line)
    if match is None:
        raise ValidationError("cannot recover seed directory/name from Wannier90 run line")
    seed_name = match.group(2)
    if re.fullmatch(r"[A-Za-z0-9_-]+", seed_name) is None:
        raise ValidationError("Wannier90 seed must be a plain basename without path/dot syntax")
    lowered = run_line.lower()
    if "wannier90.x" not in lowered:
        raise ValidationError("normal Wannier90 executable was not invoked")
    if re.search(r"\b(?:export[_-]?only|import[_-]?only|seed[_-]?only|reuse|skip)\b", lowered):
        raise ValidationError("stale-seed/export-only Wannier90 command is forbidden")
    return Path(match.group(1)).resolve(), seed_name, run_line


def _validate_provenance(log_text: str, seed_dir: Path, expected_bands: int) -> str:
    records = re.findall(
        r"\[DC-LCFO-SAWF-DMN\]\s+strict current-run seed fragment=(\d+)\s+token=(\S+)",
        log_text,
    )
    if not records:
        raise ValidationError("missing strict current-run seed provenance records")
    fragments = {int(fragment) for fragment, _ in records}
    tokens = {token for _, token in records}
    if fragments != set(range(1, 9)):
        raise ValidationError(f"current-run provenance does not cover fragments 1..8: {sorted(fragments)}")
    if len(tokens) != 1 or not next(iter(tokens)).strip():
        raise ValidationError("current-run seed records do not share one non-empty token")
    token = next(iter(tokens))

    fragment_root = seed_dir.parent / "fragments"
    for fragment in range(1, 9):
        sidecar = fragment_root / f"{fragment:06d}" / "wavefunctions_wannier_seed.provenance"
        lines = _read(sidecar, "provenance sidecar").splitlines()
        if len(lines) < 4:
            raise ValidationError(f"truncated provenance sidecar: {sidecar}")
        try:
            magic_version = [int(value) for value in lines[0].split()]
            metadata = [int(value) for value in lines[2].split()]
            basis = [int(value) for value in lines[3].split()]
        except ValueError as exc:
            raise ValidationError(f"malformed provenance sidecar: {sidecar}") from exc
        if magic_version != [-22022220, 1] or len(metadata) != 5:
            raise ValidationError(f"invalid provenance metadata shape: {sidecar}")
        if lines[1].strip() != token:
            raise ValidationError(f"stale/mixed provenance token: {sidecar}")
        if metadata[0] != 8 or metadata[1] != fragment or min(metadata[2:]) <= 0:
            raise ValidationError(f"wrong provenance fragment metadata: {sidecar}")
        if metadata[4] != expected_bands:
            raise ValidationError(f"provenance total-state dimension is not {expected_bands}: {sidecar}")
        if len(basis) != metadata[2] or any(
            count <= 0 or count > metadata[3] for count in basis
        ):
            raise ValidationError(f"wrong provenance basis dimensions: {sidecar}")
    return token


def _validate_operations(log_text: str, expected_operations: int) -> tuple[int, float]:
    pattern = re.compile(
        rf"\[DC-LCFO-SAWF-DMN\]\s+operation=(\d+)\s+"
        rf"singular_min=\s*({FLOAT})\s+singular_max=\s*({FLOAT})\s+"
        rf"closure_residual=\s*({FLOAT})\s+tolerance=\s*({FLOAT})"
    )
    records = pattern.findall(log_text)
    if len(records) != expected_operations:
        raise ValidationError(
            f"expected {expected_operations} validated symmetry operations, found {len(records)}"
        )
    if [int(record[0]) for record in records] != list(range(1, expected_operations + 1)):
        raise ValidationError(
            f"SAWF operation records are not exactly sequential 1..{expected_operations}"
        )
    tolerances = []
    for operation, singular_min, singular_max, closure, tolerance in records:
        smin = _number(singular_min, f"operation {operation} singular_min")
        smax = _number(singular_max, f"operation {operation} singular_max")
        residual = _number(closure, f"operation {operation} closure_residual")
        tol = _number(tolerance, f"operation {operation} tolerance")
        if tol <= 0.0 or smin <= 0.0 or smax < smin:
            raise ValidationError(f"invalid singular-value gate for operation {operation}")
        if max(abs(smin - 1.0), abs(smax - 1.0)) > tol:
            raise ValidationError(f"operation {operation} singular values exceed tolerance")
        if residual > 2.0 * tol:
            raise ValidationError(
                f"operation {operation} closure residual {residual:.6e} exceeds {2.0 * tol:.6e}"
            )
        tolerances.append(tol)
    logged_tolerance = tolerances[0]
    formatting_atol = max(1.0e-15, 5.0e-6 * logged_tolerance)
    if any(not math.isclose(value, logged_tolerance, rel_tol=5.0e-6, abs_tol=formatting_atol)
           for value in tolerances[1:]):
        raise ValidationError("operation records do not share one logged tolerance")
    if not math.isclose(logged_tolerance, 1.0e-6, rel_tol=5.0e-6, abs_tol=1.0e-15):
        raise ValidationError("formal C64 SAWF tolerance is not 1e-6")
    return len(records), logged_tolerance


def _validate_publication(
    log_text: str,
    operation_count: int,
    expected_bands: int,
    tolerance: float,
) -> None:
    pattern = re.compile(
        rf"\[DC-LCFO-SAWF-DMN\]\s+published operations=(\d+)\s+bands=(\d+)\s+"
        rf"unitarity_max=\s*({FLOAT})\s+hamiltonian_max=\s*({FLOAT})\s+"
        rf"amn_max=\s*({FLOAT})\s+group_wann_max=\s*({FLOAT})\s+"
        rf"group_band_max=\s*({FLOAT})"
    )
    records = pattern.findall(log_text)
    if len(records) != 1:
        raise ValidationError("expected exactly one DMN publication summary")
    operations, bands, *metric_text = records[0]
    if int(operations) != operation_count or int(bands) != expected_bands:
        raise ValidationError(
            f"published DMN dimensions are not {operation_count} operations / {expected_bands} bands"
        )
    limits = {
        "unitarity": 2.0 * tolerance,
        "hamiltonian": tolerance,
        "amn": 2.0 * tolerance,
        "group_wann": 4.0 * tolerance,
        "group_band": 4.0 * tolerance,
    }
    for (label, limit), text in zip(limits.items(), metric_text):
        value = _number(text, f"published {label}_max")
        if value > limit:
            raise ValidationError(
                f"published {label} residual {value:.6e} exceeds {limit:.6e}"
            )


def _require_nonempty(path: Path, label: str) -> None:
    try:
        if path.stat().st_size <= 0:
            raise ValidationError(f"{label} is empty: {path}")
    except OSError as exc:
        raise ValidationError(f"{label} is missing: {path}") from exc


def _run_root(seed_dir: Path) -> Path:
    if seed_dir.name != "total" or seed_dir.parent.name != "data_dcdft":
        raise ValidationError("seed directory is not RUN_ROOT/data_dcdft/total")
    return seed_dir.parent.parent.resolve()


def _input_raw_assignment(input_file: Path, name: str) -> str:
    text = _read(input_file, "SALMON input file")
    values = []
    pattern = re.compile(rf"(?i)^\s*{re.escape(name)}\s*=\s*(.*?)\s*$")
    for line in text.splitlines():
        active = line.split("!", 1)[0].strip()
        match = pattern.match(active)
        if match is not None:
            value = match.group(1).strip()
            if value.endswith(","):
                value = value[:-1].rstrip()
            values.append(value)
    if len(values) != 1 or not values[0]:
        raise ValidationError(f"input must define exactly one non-empty {name}")
    return values[0]


def _input_string_assignment(input_file: Path, name: str) -> str:
    raw = _input_raw_assignment(input_file, name)
    match = re.fullmatch(r"(['\"])(.*?)\1", raw)
    if match is None or not match.group(2).strip():
        raise ValidationError(f"input {name} must be one non-empty quoted string")
    return match.group(2).strip()


def _input_integer_assignment(input_file: Path, name: str) -> int:
    raw = _input_raw_assignment(input_file, name)
    if re.fullmatch(r"[+-]?\d+", raw) is None:
        raise ValidationError(f"input {name} must be one integer")
    return int(raw)


def _input_integer_vector(input_file: Path, name: str) -> tuple[int, int, int]:
    raw = _input_raw_assignment(input_file, name)
    fields = [field for field in re.split(r"[\s,]+", raw) if field]
    if len(fields) != 3 or any(re.fullmatch(r"[+-]?\d+", field) is None for field in fields):
        raise ValidationError(f"input {name} must contain exactly three integers")
    return tuple(int(field) for field in fields)


def _resolve_run_root_file(
    run_root: Path,
    configured_name: str,
    expected_name: str,
    legacy_name: str,
    label: str,
) -> Path:
    if configured_name != expected_name:
        raise ValidationError(f"formal aligned input must name {expected_name} exactly")
    candidate = Path(configured_name)
    resolved = (candidate if candidate.is_absolute() else run_root / candidate).resolve()
    if resolved.parent != run_root or resolved.name != expected_name:
        raise ValidationError(f"formal aligned {label} escapes the recovered run root")
    _require_nonempty(resolved, f"aligned {label}")
    legacy = run_root / legacy_name
    if legacy.exists() and resolved.samefile(legacy):
        raise ValidationError(f"aligned {label} aliases legacy run-root {legacy_name}")
    return resolved


def _parse_atom_rows(
    path: Path,
) -> tuple[tuple[str, tuple[float, float, float], str, tuple[str, ...]], ...]:
    rows = []
    for line_number, line in enumerate(_read(path, "aligned atom file").splitlines(), 1):
        active = line.split("#", 1)[0].strip()
        if not active:
            continue
        try:
            fields = shlex.split(active, posix=False)
        except ValueError as exc:
            raise ValidationError(f"malformed atom row {line_number}: {exc}") from exc
        if len(fields) < 5:
            raise ValidationError(f"atom row {line_number} has fewer than five fields")
        coordinates = tuple(
            _number(fields[index], f"atom row {line_number} coordinate")
            for index in range(1, 4)
        )
        rows.append((fields[0], coordinates, fields[4], tuple(fields[5:])))
    return tuple(rows)


def _validate_formal_atom(actual_path: Path) -> None:
    reference_path = FORMAL_SAMPLE / FORMAL_ATOM_NAME
    try:
        reference_digest = hashlib.sha256(reference_path.read_bytes()).hexdigest()
    except OSError as exc:
        raise ValidationError(f"cannot read committed aligned atom reference: {exc}") from exc
    if reference_digest != FORMAL_ATOM_SHA256:
        raise ValidationError("committed aligned atom reference SHA-256 is not trusted")
    actual = _parse_atom_rows(actual_path)
    reference = _parse_atom_rows(reference_path)
    if len(actual) != 64 or len(reference) != 64:
        raise ValidationError("formal aligned atom geometry must contain exactly 64 rows")
    cell = 13.44
    tolerance = 1.0e-12
    for index, (actual_row, reference_row) in enumerate(zip(actual, reference), 1):
        if (actual_row[0], actual_row[2], actual_row[3]) != (
            reference_row[0], reference_row[2], reference_row[3]
        ):
            raise ValidationError(f"aligned atom metadata differs at row {index}")
        for actual_value, reference_value in zip(actual_row[1], reference_row[1]):
            difference = abs((actual_value - reference_value) % cell)
            if min(difference, cell - difference) > tolerance:
                raise ValidationError(f"aligned atom coordinate differs at row {index}")


def _resolve_symmetry_provenance(
    seed_dir: Path,
    input_file: Path | None,
    symmetry_file: Path | None,
) -> tuple[Path, str]:
    run_root = _run_root(seed_dir)
    if input_file is None:
        raise ValidationError("formal aligned validation requires --input-file")
    input_file = Path(input_file).resolve()
    if input_file.parent != run_root or input_file.name != FORMAL_INPUT_NAME:
        raise ValidationError(
            f"formal input must be RUN_ROOT/{FORMAL_INPUT_NAME}"
        )
    _require_nonempty(input_file, "formal aligned input file")

    expected_strings = {
        "sysname": FORMAL_SYSNAME,
        "wannier_site_symmetry": "file",
        "wannier_symmetry_file": FORMAL_SYMMETRY_NAME,
        "file_atom_coor": FORMAL_ATOM_NAME,
        "yn_eigenexa": "y",
    }
    strings = {
        name: _input_string_assignment(input_file, name)
        for name in expected_strings
    }
    for name, expected in expected_strings.items():
        if strings[name] != expected:
            raise ValidationError(f"formal input {name} is not {expected!r}")

    expected_integers = {
        "wannier_num_wann": 576,
        "wannier_num_bands": 664,
        "nstate": 664,
    }
    for name, expected in expected_integers.items():
        if _input_integer_assignment(input_file, name) != expected:
            raise ValidationError(f"formal input {name} is not {expected}")
    if any(value == 256 for value in (
        _input_integer_assignment(input_file, "wannier_num_wann"),
        _input_integer_assignment(input_file, "wannier_num_bands"),
        _input_integer_assignment(input_file, "nstate"),
    )):
        raise ValidationError("formal aligned input regressed to a 256-state basis")

    expected_vectors = {
        "num_fragment(1:3)": (2, 2, 2),
        "num_rgrid_buffer(1:3)": (6, 6, 6),
        "num_rgrid(1:3)": (32, 32, 32),
    }
    for name, expected in expected_vectors.items():
        if _input_integer_vector(input_file, name) != expected:
            raise ValidationError(f"formal input {name} is not {expected}")

    named_path = _resolve_run_root_file(
        run_root, strings["wannier_symmetry_file"], FORMAL_SYMMETRY_NAME,
        "sym.dat", "symmetry file",
    )
    atom_path = _resolve_run_root_file(
        run_root, strings["file_atom_coor"], FORMAL_ATOM_NAME,
        "atom.dat", "atom file",
    )
    _validate_formal_atom(atom_path)

    explicit_path = None
    if symmetry_file is not None:
        candidate = Path(symmetry_file)
        explicit_path = (candidate if candidate.is_absolute() else run_root / candidate).resolve()

    if explicit_path is not None and named_path != explicit_path:
        raise ValidationError("explicit symmetry file differs from the input-file setting")
    resolved = named_path
    if resolved.parent != run_root or resolved.name != "sym_sawf_aligned.dat":
        raise ValidationError(
            "formal aligned symmetry must be the dedicated run-root sym_sawf_aligned.dat"
        )
    legacy = run_root / "sym.dat"
    if legacy.exists() and resolved.exists() and resolved.samefile(legacy):
        raise ValidationError("aligned symmetry file aliases legacy run-root sym.dat")
    _require_nonempty(resolved, "aligned symmetry file")
    return resolved, strings["sysname"]


def _validate_symmetry_file(sym_path: Path, expected_operations: int) -> None:
    if expected_operations != 2:
        raise ValidationError("formal aligned C64 symmetry gate requires exactly 2 operations")
    rows = []
    for line_number, line in enumerate(_read(sym_path, "aligned symmetry file").splitlines(), 1):
        comment_positions = [position for position in (line.find("#"), line.find("!"))
                             if position >= 0]
        if comment_positions:
            line = line[:min(comment_positions)]
        if not line.strip():
            continue
        fields = line.split()
        if len(fields) != 4:
            raise ValidationError(f"aligned symmetry row {line_number} does not contain four values")
        rows.append([_number(field, f"aligned symmetry row {line_number}") for field in fields])
    if len(rows) != 6:
        raise ValidationError(
            f"formal aligned C64 symmetry must contain 6 active rows, found {len(rows)}"
        )
    operations = [rows[offset:offset + 3] for offset in range(0, len(rows), 3)]
    identity = [[1.0, 0.0, 0.0, 0.0],
                [0.0, 1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0]]
    inversion = [[-1.0, 0.0, 0.0, 15.0 / 32.0],
                 [0.0, -1.0, 0.0, 15.0 / 32.0],
                 [0.0, 0.0, -1.0, 15.0 / 32.0]]

    def same(left: list[list[float]], right: list[list[float]]) -> bool:
        return all(math.isclose(left[row][column], right[row][column],
                                rel_tol=0.0, abs_tol=1.0e-12)
                   for row in range(3) for column in range(4))

    if sum(same(operation, identity) for operation in operations) != 1:
        raise ValidationError("formal aligned C64 symmetry lacks exactly one identity operation")
    if sum(same(operation, inversion) for operation in operations) != 1:
        raise ValidationError(
            "formal aligned C64 symmetry lacks inversion W=-I, tau=(15/32)^3"
        )


def _validate_alignment_records(
    log_text: str, expected_operations: int, tolerance: float
) -> None:
    tagged_lines = [
        line for line in log_text.splitlines()
        if "[DC-LCFO-SAWF-ALIGN]" in line
    ]
    if len(tagged_lines) != expected_operations:
        raise ValidationError(
            f"expected exactly {expected_operations} SAWF alignment lines, "
            f"found {len(tagged_lines)}"
        )
    pattern = re.compile(
        rf"(?m)^\s*\[DC-LCFO-SAWF-ALIGN\]\s+operation=(\d+)\s+"
        rf"grid_map_ok=([TF])\s+fragment_map_ok=([TF])\s+"
        rf"max_targets_per_source=(\d+)\s+max_grid_residual=\s*({FLOAT})"
        rf"(?:\s+center_grid=\s*({FLOAT})\s+({FLOAT})\s+({FLOAT}))?\s*$"
    )
    records = list(pattern.finditer(log_text))
    if len(records) != expected_operations:
        raise ValidationError(
            f"expected {expected_operations} bounded SAWF alignment records, found {len(records)}"
        )
    if [int(record.group(1)) for record in records] != list(
        range(1, expected_operations + 1)
    ):
        raise ValidationError("SAWF alignment records are not exactly sequential 1..2")
    for record in records:
        operation = int(record.group(1))
        if record.group(2) != "T" or record.group(3) != "T":
            raise ValidationError(f"operation {operation} failed grid/fragment alignment")
        if int(record.group(4)) != 1:
            raise ValidationError(f"operation {operation} splits a fragment core")
        residual = _number(record.group(5), f"operation {operation} grid residual")
        if residual < 0.0 or residual > tolerance:
            raise ValidationError(f"operation {operation} grid residual exceeds tolerance")
        center = record.groups()[5:8]
        if operation == 1 and any(value is not None for value in center):
            raise ValidationError("identity alignment record unexpectedly reports a center")
        if operation == 2:
            if any(value is None for value in center):
                raise ValidationError("inversion alignment record lacks center_grid")
            center_values = tuple(
                _number(value, "inversion center_grid") for value in center
            )
            if any(not math.isclose(value, 7.5, rel_tol=0.0, abs_tol=tolerance)
                   for value in center_values):
                raise ValidationError(
                    "formal aligned inversion center_grid is not half-grid (7.5,7.5,7.5)"
                )

    d_band = list(re.finditer(r"\[DC-LCFO-SAWF-DMN\]\s+operation=\d+", log_text))
    if len(d_band) != expected_operations:
        raise ValidationError("cannot establish alignment-before-D_band ordering")
    if records[-1].start() >= d_band[0].start():
        raise ValidationError("SAWF alignment validation occurs after D_band operation logs")


def _validate_artifacts(
    seed_dir: Path,
    seed_name: str,
    expected_operations: int,
    expected_bands: int,
    logged_tolerance: float,
) -> None:
    dmn = seed_dir / f"{seed_name}.dmn"
    _require_nonempty(dmn, "published DMN")
    try:
        with dmn.open(errors="replace") as handle:
            header = []
            while len(header) < 2:
                line = handle.readline()
                if line == "":
                    break
                if line.strip():
                    header.append(line.strip())
    except OSError as exc:
        raise ValidationError(f"cannot read published DMN header: {dmn}") from exc
    try:
        dimensions = [int(value) for value in header[1].split()]
    except (IndexError, ValueError) as exc:
        raise ValidationError(f"malformed published DMN header: {dmn}") from exc
    expected_dimensions = [expected_bands, expected_operations, 1, 1]
    if dimensions != expected_dimensions:
        raise ValidationError(
            f"published DMN header is not {expected_bands}/{expected_operations}/Gamma-only: {dimensions}"
        )

    win = _read(seed_dir / f"{seed_name}.win", "activated win file")
    active = [line.strip().lower() for line in win.splitlines()
              if line.lstrip() and not line.lstrip().startswith(("!", "#"))]
    if active.count("site_symmetry = true") != 1:
        raise ValidationError("activated .win must contain site_symmetry = true exactly once")
    eps = [line for line in active if re.match(r"^symmetrize_eps\s*=", line)]
    if len(eps) != 1:
        raise ValidationError("activated .win must contain one finite positive symmetrize_eps")
    symmetrize_eps = _number(eps[0].split("=", 1)[1].strip(), "symmetrize_eps")
    if symmetrize_eps <= 0.0 or not math.isclose(
        symmetrize_eps, logged_tolerance, rel_tol=5.0e-6, abs_tol=1.0e-15
    ):
        raise ValidationError("activated .win symmetrize_eps differs from logged tolerance")
    if not math.isclose(symmetrize_eps, 1.0e-6, rel_tol=5.0e-6, abs_tol=1.0e-15):
        raise ValidationError("formal C64 .win symmetrize_eps is not 1e-6")
    num_wann = [line for line in active if re.match(r"^num_wann\s*=", line)]
    num_bands = [line for line in active if re.match(r"^num_bands\s*=", line)]
    if num_wann != ["num_wann = 576"] or num_bands != [f"num_bands = {expected_bands}"]:
        raise ValidationError(
            f"activated .win dimensions are not 576 Wannier / {expected_bands} bands"
        )

    wout = _read(seed_dir / f"{seed_name}.wout", "Wannier90 output")
    if not re.search(r"Using symmetry-adapted WF mode\s*:\s*T\b", wout):
        raise ValidationError("Wannier90 did not report symmetry-adapted WF mode")
    if "All done: wannier90 exiting" not in wout:
        raise ValidationError("Wannier90 output has no successful completion marker")
    wout_bands = re.search(r"Number of input Bloch states\s*:\s*(\d+)", wout)
    wout_wann = re.search(r"Number of Wannier Functions\s*:\s*(\d+)", wout)
    if wout_bands is None or int(wout_bands.group(1)) != expected_bands:
        raise ValidationError(f"Wannier90 output does not report {expected_bands} input bands")
    if wout_wann is None or int(wout_wann.group(1)) != 576:
        raise ValidationError("Wannier90 output does not report 576 Wannier functions")
    for artifact, label in (
        (seed_dir / f"{seed_name}.chk", "Wannier90 checkpoint"),
        (seed_dir / "wannier90_global_basis.bin", "global Wannier basis"),
        (seed_dir / "wannier_flux_eigen_seed.bin", "Wannier eigen seed"),
    ):
        _require_nonempty(artifact, label)


def _validate_post_wannier_path(log_text: str) -> None:
    markers = [
        "[DC-LCFO-WANNIER] Wannier90 completed.",
        "[DC-LCFO-W90-GLOBAL] wrote 576 Wannier functions to ",
        "[DC-LCFO-W90-SEED] wrote Flux-LCFO eigen seed in Wannier basis: states=576 wann=576 file=",
    ]
    positions = []
    for marker in markers:
        if log_text.count(marker) != 1:
            raise ValidationError(f"post-Wannier gate must occur exactly once: {marker}")
        positions.append(log_text.index(marker))
    end_positions = [match.start() for match in re.finditer(r"(?m)^\s*end SALMON\s*$", log_text)]
    if len(end_positions) != 1:
        raise ValidationError("expected exactly one final end SALMON marker")
    if not positions[0] < positions[1] < positions[2] < end_positions[0]:
        raise ValidationError("checkpoint import/global basis/eigen seed/end SALMON are out of order")
    later = log_text[positions[0]:]
    fatal_pattern = re.compile(
        r"(?i)(?:\bfatal\b|\berror(?:\s+stop)?\b|\babnormal\b|\bimproperly\b|"
        r"mpi_abort|exited due|segmentation fault|traceback)"
    )
    match = fatal_pattern.search(later)
    if match is not None:
        line = later[match.start():].splitlines()[0]
        raise ValidationError(f"fatal/error text appears after Wannier90 completion: {line}")


def validate_integration(
    log_path: Path,
    expected_operations: int = 2,
    expected_bands: int = 664,
    *,
    input_file: Path | None = None,
    symmetry_file: Path | None = None,
) -> None:
    if expected_operations <= 0 or expected_bands <= 0:
        raise ValidationError("expected operation and band counts must be positive")
    if expected_operations != 2 or expected_bands != 664:
        raise ValidationError("formal aligned C64 gate requires exactly 2 operations / 664 bands")
    log_text = _read(log_path, "SALMON integration log")
    relevant = "\n".join(line for line in log_text.splitlines()
                         if "SAWF" in line or "WANNIER" in line)
    if re.search(r"(?<![A-Za-z])(?:nan|[+-]?inf(?:inity)?)(?![A-Za-z])", relevant, re.I):
        raise ValidationError("SAWF/Wannier integration log contains NaN or Inf")

    seed_dir, seed_name, _ = _seed_location(log_text)
    symmetry_path, input_sysname = _resolve_symmetry_provenance(
        seed_dir, input_file, symmetry_file
    )
    if seed_name != input_sysname:
        raise ValidationError("Wannier90 run seed differs from formal input sysname")
    _validate_symmetry_file(symmetry_path, expected_operations)
    token = _validate_provenance(log_text, seed_dir, expected_bands)
    operation_count, tolerance = _validate_operations(log_text, expected_operations)
    _validate_alignment_records(log_text, operation_count, tolerance)
    _validate_publication(log_text, operation_count, expected_bands, tolerance)

    required = [
        "[DC-LCFO-SAWF] SAWF pseudo_channels ordering: complete atom-major s+p+d shells",
        "[DC-LCFO-SAWF] .amn keeps shared complete atom-major s+p+d ordering",
        "[DC-LCFO-SAWF] activated site_symmetry in",
        "[DC-LCFO-WANNIER] Wannier90 completed.",
    ]
    for marker in required:
        if log_text.count(marker) != 1:
            raise ValidationError(f"required integration gate must occur exactly once: {marker}")
    publication_position = log_text.index("[DC-LCFO-SAWF-DMN] published operations=")
    run_position = log_text.index("[DC-LCFO-WANNIER] run:")
    completion_position = log_text.index("[DC-LCFO-WANNIER] Wannier90 completed.")
    activation_position = log_text.index("[DC-LCFO-SAWF] activated site_symmetry in")
    if not publication_position < activation_position < run_position < completion_position:
        raise ValidationError("DMN publication, .win activation, and Wannier90 execution are out of order")

    _validate_post_wannier_path(log_text)
    _validate_artifacts(
        seed_dir, seed_name, expected_operations, expected_bands, tolerance
    )
    print(
        f"PASS SAWF integration: token={token} operations={operation_count} "
        f"bands={expected_bands} tolerance={tolerance:.3e} seed={seed_name}"
    )


def _write_good_fixture(
    root: Path,
    expected_operations: int,
    expected_bands: int,
) -> Path:
    seed_dir = root / "data_dcdft" / "total"
    fragment_root = root / "data_dcdft" / "fragments"
    seed_dir.mkdir(parents=True)
    token = "2026-7-10-12-34-56-123"
    for fragment in range(1, 9):
        directory = fragment_root / f"{fragment:06d}"
        directory.mkdir(parents=True)
        (directory / "wavefunctions_wannier_seed.provenance").write_text(
            f"-22022220 1\n{token}\n8 {fragment} 1 400 {expected_bands}\n72\n"
        )
    seed = FORMAL_SYSNAME
    (root / FORMAL_SYMMETRY_NAME).write_text(
        "# Generated by align_periodic_structure_to_fragments.py\n"
        "# translation=11/64 11/64 11/64\n"
        "# buffer=6 6 6\n"
        "1 0 0 0\n0 1 0 0\n0 0 1 0\n"
        "-1 0 0 0.46875\n0 -1 0 0.46875\n0 0 -1 0.46875\n"
    )
    (root / FORMAL_ATOM_NAME).write_bytes(
        (FORMAL_SAMPLE / FORMAL_ATOM_NAME).read_bytes()
    )
    (root / FORMAL_INPUT_NAME).write_text(
        "&control\n"
        f"  sysname = '{seed}'\n/\n"
        "&dc\n"
        "  num_fragment(1:3) = 2,2,2\n"
        "  num_rgrid_buffer(1:3) = 6,6,6\n"
        "  wannier_num_wann = 576\n"
        f"  wannier_num_bands = {expected_bands}\n"
        "  wannier_site_symmetry = 'file'\n"
        f"  wannier_symmetry_file = '{FORMAL_SYMMETRY_NAME}',\n"
        "  wannier_symmetry_tolerance = 1.0d-6\n/\n"
        "&parallel\n"
        "  yn_eigenexa = 'y'\n/\n"
        "&system\n"
        f"  nstate = {expected_bands}\n"
        f"  file_atom_coor = '{FORMAL_ATOM_NAME}'\n/\n"
        "&rgrid\n"
        "  num_rgrid(1:3) = 32,32,32\n/\n"
    )
    (seed_dir / f"{seed}.dmn").write_text(
        "SALMON SAWF integration fixture\n"
        f"       {expected_bands}        {expected_operations}         1         1\n"
    )
    (seed_dir / f"{seed}.win").write_text(
        f"num_wann = 576\nnum_bands = {expected_bands}\n"
        "site_symmetry = true\nsymmetrize_eps = 1d-6\n"
    )
    (seed_dir / f"{seed}.wout").write_text(
        "Number of Wannier Functions : 576\n"
        f"Number of input Bloch states : {expected_bands}\n"
        "Using symmetry-adapted WF mode : T\nAll done: wannier90 exiting\n"
    )
    for artifact in (
        f"{seed}.chk",
        "wannier90_global_basis.bin",
        "wannier_flux_eigen_seed.bin",
    ):
        (seed_dir / artifact).write_bytes(b"validated artifact\n")
    lines = [
        "[DC-LCFO-SAWF] SAWF pseudo_channels ordering: complete atom-major s+p+d shells",
        "[DC-LCFO-SAWF] .amn keeps shared complete atom-major s+p+d ordering",
    ]
    lines.extend(
        f"[DC-LCFO-SAWF-DMN] strict current-run seed fragment={fragment} token={token}"
        for fragment in range(1, 9)
    )
    lines.extend([
        "[DC-LCFO-SAWF-ALIGN] operation=1 grid_map_ok=T fragment_map_ok=T "
        "max_targets_per_source=1 max_grid_residual= 0.00000E+00",
        "[DC-LCFO-SAWF-ALIGN] operation=2 grid_map_ok=T fragment_map_ok=T "
        "max_targets_per_source=1 max_grid_residual= 0.00000E+00 "
        "center_grid=     7.50000     7.50000     7.50000",
    ])
    lines.extend(
        f"[DC-LCFO-SAWF-DMN] operation={operation} singular_min= 9.999995E-01 "
        "singular_max= 1.0000005E+00 closure_residual= 1.50000E-06 tolerance= 1.00000E-06"
        for operation in range(1, expected_operations + 1)
    )
    lines.extend([
        f"[DC-LCFO-SAWF-DMN] published operations={expected_operations} bands={expected_bands} "
        "unitarity_max= 1.50000E-06 hamiltonian_max= 9.00000E-07 "
        "amn_max= 1.50000E-06 group_wann_max= 3.50000E-06 group_band_max= 3.50000E-06",
        f"[DC-LCFO-SAWF] activated site_symmetry in {seed_dir / (seed + '.win')}",
        f"[DC-LCFO-WANNIER] run: cd '{seed_dir}/' && env -u PMIX_RANK "
        f"'/opt/salmon/wannier90.x' '{seed}'",
        "[DC-LCFO-WANNIER] Wannier90 completed.",
        f"[DC-LCFO-W90-GLOBAL] wrote 576 Wannier functions to "
        f"{seed_dir / 'wannier90_global_basis.bin'}",
        "[DC-LCFO-W90-SEED] wrote Flux-LCFO eigen seed in Wannier basis: "
        f"states=576 wann=576 file={seed_dir / 'wannier_flux_eigen_seed.bin'}",
        "end SALMON",
    ])
    log = root / "run.log"
    log.write_text("\n".join(lines) + "\n")
    return log


def self_test() -> None:
    with tempfile.TemporaryDirectory() as tmp_name:
        root = Path(tmp_name)
        expected_operations = 2
        expected_bands = 664
        log = _write_good_fixture(root, expected_operations, expected_bands)
        input_file = root / FORMAL_INPUT_NAME
        symmetry_file = root / FORMAL_SYMMETRY_NAME
        atom_file = root / FORMAL_ATOM_NAME
        validate_integration(
            log, expected_operations, expected_bands,
            input_file=input_file, symmetry_file=symmetry_file,
        )
        validate_integration(
            log, expected_operations, expected_bands, input_file=input_file
        )

        try:
            validate_integration(
                log, expected_operations, expected_bands,
                symmetry_file=symmetry_file,
            )
        except ValidationError:
            pass
        else:
            raise AssertionError("self-test accepted explicit symmetry without input")

        original = log.read_text()
        seed_dir = root / "data_dcdft/total"
        seed = f"c64_dc_pseudo_sawf_aligned_nw576_nb{expected_bands}"

        def reject(label: str) -> None:
            try:
                validate_integration(
                    log, expected_operations, expected_bands,
                    input_file=input_file, symmetry_file=symmetry_file,
                )
            except ValidationError:
                return
            raise AssertionError(f"self-test accepted {label}")

        cases = {
            "missing operation": original.replace(
                "[DC-LCFO-SAWF-DMN] operation=2", "[REMOVED] operation=2", 1
            ),
            "non-finite metric": original.replace(
                "unitarity_max= 1.50000E-06", "unitarity_max= NaN", 1
            ),
            "failed closure": original.replace(
                "closure_residual= 1.50000E-06", "closure_residual= 2.10000E-06", 1
            ),
            "failed singular bound": original.replace(
                "singular_min= 9.999995E-01", "singular_min= 9.999980E-01", 1
            ),
            "wrong publication count": original.replace(
                "published operations=2", "published operations=3", 1
            ),
            "wrong publication bands": original.replace(
                "published operations=2 bands=664", "published operations=2 bands=640", 1
            ),
            "unitarity above 2 tolerance": original.replace(
                "unitarity_max= 1.50000E-06", "unitarity_max= 2.10000E-06", 1
            ),
            "Hamiltonian above tolerance": original.replace(
                "hamiltonian_max= 9.00000E-07", "hamiltonian_max= 1.10000E-06", 1
            ),
            "AMN above 2 tolerance": original.replace(
                "amn_max= 1.50000E-06", "amn_max= 2.10000E-06", 1
            ),
            "Wannier group above 4 tolerance": original.replace(
                "group_wann_max= 3.50000E-06", "group_wann_max= 4.10000E-06", 1
            ),
            "band group above 4 tolerance": original.replace(
                "group_band_max= 3.50000E-06", "group_band_max= 4.10000E-06", 1
            ),
            "inconsistent logged tolerance": original.replace(
                "operation=2 singular_min= 9.999995E-01 singular_max= 1.0000005E+00 "
                "closure_residual= 1.50000E-06 tolerance= 1.00000E-06",
                "operation=2 singular_min= 9.999995E-01 singular_max= 1.0000005E+00 "
                "closure_residual= 1.50000E-06 tolerance= 1.10000E-06",
                1,
            ),
            "missing global basis tag": original.replace(
                "[DC-LCFO-W90-GLOBAL] wrote 576", "[REMOVED-W90-GLOBAL] wrote 576", 1
            ),
            "missing eigen seed tag": original.replace(
                "[DC-LCFO-W90-SEED] wrote", "[REMOVED-W90-SEED] wrote", 1
            ),
            "missing final marker": original.replace("end SALMON", "[REMOVED END]", 1),
            "later fatal": original + "FATAL post-Wannier failure\n",
            "stale route": original.replace("'/opt/salmon/wannier90.x'", "'reuse'", 1),
            "missing alignment line": original.replace(
                "[DC-LCFO-SAWF-ALIGN] operation=2", "[REMOVED-SAWF-ALIGN] operation=2", 1
            ),
            "integer inversion center": original.replace(
                "7.50000     7.50000     7.50000",
                "7.00000     7.00000     7.00000",
                1,
            ),
            "q31 inversion center": original.replace(
                "7.50000     7.50000     7.50000",
                "15.50000    15.50000    15.50000",
                1,
            ),
            "split fragment map": original.replace(
                "operation=2 grid_map_ok=T fragment_map_ok=T max_targets_per_source=1",
                "operation=2 grid_map_ok=T fragment_map_ok=F max_targets_per_source=8",
                1,
            ),
            "negative grid residual": original.replace(
                "max_grid_residual= 0.00000E+00",
                "max_grid_residual=-1.00000E-12",
                1,
            ),
            "alignment after D_band": original.replace(
                "[DC-LCFO-SAWF-ALIGN] operation=2 grid_map_ok=T fragment_map_ok=T "
                "max_targets_per_source=1 max_grid_residual= 0.00000E+00 "
                "center_grid=     7.50000     7.50000     7.50000\n",
                "", 1
            ).replace(
                "[DC-LCFO-SAWF-DMN] operation=2 singular_min=",
                "[DC-LCFO-SAWF-DMN] operation=2 singular_min=",
                1,
            ) + (
                "[DC-LCFO-SAWF-ALIGN] operation=2 grid_map_ok=T fragment_map_ok=T "
                "max_targets_per_source=1 max_grid_residual= 0.00000E+00 "
                "center_grid=     7.50000     7.50000     7.50000\n"
            ),
        }
        for label, broken in cases.items():
            if broken == original:
                raise AssertionError(f"self-test mutation did not change fixture: {label}")
            log.write_text(broken)
            reject(label)
        log.write_text(original)

        dmn = seed_dir / f"{seed}.dmn"
        original_dmn = dmn.read_text()
        dmn.write_text(original_dmn.replace("664        2", "664        3"))
        reject("wrong DMN operation dimension")
        dmn.write_text(original_dmn.replace("664        2", "640        2"))
        reject("wrong DMN band dimension")
        dmn.write_text(original_dmn)

        win = seed_dir / f"{seed}.win"
        original_win = win.read_text()
        win.write_text(original_win.replace("symmetrize_eps = 1d-6", "symmetrize_eps = 2d-6"))
        reject("win tolerance mismatch")
        win.write_text(original_win.replace("num_bands = 664", "num_bands = 640"))
        reject("win band dimension mismatch")
        win.write_text(original_win)

        wout = seed_dir / f"{seed}.wout"
        original_wout = wout.read_text()
        wout.write_text(original_wout.replace("input Bloch states : 664", "input Bloch states : 640"))
        reject("Wannier90 input-band mismatch")
        wout.write_text(original_wout)

        for artifact in (
            seed_dir / f"{seed}.chk",
            seed_dir / "wannier90_global_basis.bin",
            seed_dir / "wannier_flux_eigen_seed.bin",
        ):
            contents = artifact.read_bytes()
            artifact.unlink()
            reject(f"missing artifact {artifact.name}")
            artifact.write_bytes(contents)

        sym_path = symmetry_file
        original_sym = sym_path.read_text()
        bad_symmetries = {
            "wrong inversion translation": original_sym.replace("0.46875", "0.125", 1),
            "non-finite symmetry": original_sym.replace("0.46875", "NaN", 1),
            "extra symmetry operation": original_sym + "1 0 0 0\n0 1 0 0\n0 0 1 0\n",
            "malformed symmetry row": original_sym.replace("-1 0 0 0.46875", "-1 0 0"),
        }
        for label, broken in bad_symmetries.items():
            sym_path.write_text(broken)
            reject(label)
        sym_path.write_text(original_sym)

        original_input = input_file.read_text()
        legacy_sym = root / "sym.dat"
        legacy_sym.write_text(original_sym.replace("0.46875", "0.125"))
        input_file.write_text(
            original_input.replace("sym_sawf_aligned.dat", "sym.dat")
        )
        try:
            validate_integration(
                log, expected_operations, expected_bands, input_file=input_file
            )
        except ValidationError:
            pass
        else:
            raise AssertionError("self-test accepted legacy symmetry path")
        input_file.write_text(
            original_input.replace(
                "sym_sawf_aligned.dat", "../sym_sawf_aligned.dat"
            )
        )
        try:
            validate_integration(
                log, expected_operations, expected_bands, input_file=input_file
            )
        except ValidationError:
            pass
        else:
            raise AssertionError("self-test accepted symmetry path traversal")
        input_file.write_text(original_input)
        try:
            validate_integration(
                log, expected_operations, expected_bands,
                input_file=input_file, symmetry_file=legacy_sym,
            )
        except ValidationError:
            pass
        else:
            raise AssertionError("self-test accepted explicit legacy symmetry path")

        legacy_sym.write_text(original_sym)
        sym_path.unlink()
        sym_path.hardlink_to(legacy_sym)
        try:
            validate_integration(
                log, expected_operations, expected_bands,
                input_file=input_file, symmetry_file=sym_path,
            )
        except ValidationError:
            pass
        else:
            raise AssertionError("self-test accepted hard-link alias to legacy sym.dat")
        sym_path.unlink()
        sym_path.write_text(original_sym)

        input_mutations = {
            "sysname mismatch": original_input.replace(FORMAL_SYSNAME, "misleading_seed", 1),
            "site symmetry disabled": original_input.replace(
                "wannier_site_symmetry = 'file'", "wannier_site_symmetry = 'none'", 1
            ),
            "site symmetry case mismatch": original_input.replace(
                "wannier_site_symmetry = 'file'", "wannier_site_symmetry = 'FILE'", 1
            ),
            "legacy atom name": original_input.replace(FORMAL_ATOM_NAME, "atom.dat", 1),
            "atom path traversal": original_input.replace(
                FORMAL_ATOM_NAME, f"../{FORMAL_ATOM_NAME}", 1
            ),
            "256 Wannier functions": original_input.replace(
                "wannier_num_wann = 576", "wannier_num_wann = 256", 1
            ),
            "256 bands": original_input.replace(
                "wannier_num_bands = 664", "wannier_num_bands = 256", 1
            ),
            "256 states": original_input.replace("nstate = 664", "nstate = 256", 1),
            "wrong fragments": original_input.replace(
                "num_fragment(1:3) = 2,2,2", "num_fragment(1:3) = 1,2,2", 1
            ),
            "wrong buffer": original_input.replace(
                "num_rgrid_buffer(1:3) = 6,6,6",
                "num_rgrid_buffer(1:3) = 5,6,6",
                1,
            ),
            "wrong mesh": original_input.replace(
                "num_rgrid(1:3) = 32,32,32", "num_rgrid(1:3) = 24,32,32", 1
            ),
            "EigenExa disabled": original_input.replace(
                "yn_eigenexa = 'y'", "yn_eigenexa = 'n'", 1
            ),
        }
        for label, broken in input_mutations.items():
            input_file.write_text(broken)
            reject(label)
        input_file.write_text(original_input)

        atom_contents = atom_file.read_bytes()
        atom_file.unlink()
        reject("missing aligned atom")
        atom_file.write_bytes((FORMAL_SAMPLE / "atom.dat").read_bytes())
        reject("unshifted aligned atom")
        atom_file.write_bytes(
            atom_contents.replace(b"2.31", b"2.32", 1)
        )
        reject("stale aligned atom coordinate")
        atom_file.write_bytes(atom_contents)

        legacy_atom = root / "atom.dat"
        legacy_atom.write_bytes(atom_contents)
        atom_file.unlink()
        atom_file.hardlink_to(legacy_atom)
        reject("aligned atom hard-link alias to legacy atom.dat")
        atom_file.unlink()
        atom_file.write_bytes(atom_contents)

        global FORMAL_ATOM_SHA256
        expected_atom_sha256 = FORMAL_ATOM_SHA256
        FORMAL_ATOM_SHA256 = "0" * 64
        reject("untrusted committed atom reference")
        FORMAL_ATOM_SHA256 = expected_atom_sha256

        external_directory = root / "external"
        external_directory.mkdir()
        external_input = external_directory / FORMAL_INPUT_NAME
        external_input.write_text(original_input)
        try:
            validate_integration(
                log, expected_operations, expected_bands,
                input_file=external_input, symmetry_file=symmetry_file,
            )
        except ValidationError:
            pass
        else:
            raise AssertionError("self-test accepted external formal input")

        for bad_seed in ("../seed", "/absolute/seed", "seed.name", r"folder\seed"):
            bad_log = original.replace(f"'{seed}'", f"'{bad_seed}'", 1)
            try:
                _seed_location(bad_log)
            except ValidationError:
                pass
            else:
                raise AssertionError(f"self-test accepted unsafe seed {bad_seed!r}")

        renamed_seed = "misleading_seed"
        for suffix in ("dmn", "win", "wout", "chk"):
            (seed_dir / f"{renamed_seed}.{suffix}").write_bytes(
                (seed_dir / f"{seed}.{suffix}").read_bytes()
            )
        log.write_text(original.replace(f"'{seed}'", f"'{renamed_seed}'", 1))
        reject("run seed differs from input sysname")
        log.write_text(original)

        sidecar = root / "data_dcdft/fragments/000008/wavefunctions_wannier_seed.provenance"
        original_sidecar = sidecar.read_text()
        sidecar.write_text(original_sidecar.replace("-22022220 1", "-22022221 1"))
        reject("wrong provenance magic")
        sidecar.write_text(original_sidecar.replace("8 8 1 400 664", "8 8 2 400 664"))
        reject("wrong provenance basis dimensions")
        sidecar.write_text(original_sidecar.replace("8 8 1 400 664", "8 8 1 400 640"))
        reject("wrong provenance total-state dimension")
        sidecar.write_text(original_sidecar)
        sidecar.write_text(sidecar.read_text().replace("2026-7-10-12-34-56-123", "stale"))
        reject("mixed provenance")
        sidecar.write_text(original_sidecar)
        print("PASS SAWF integration checker rejection self-tests")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("log", nargs="?", type=Path)
    parser.add_argument(
        "--expected-operations",
        type=int,
        default=2,
        help="required symmetry operation count (default: 2 for the formal C64 sample)",
    )
    parser.add_argument(
        "--expected-bands",
        type=int,
        default=664,
        help="required full LCFO band count (default: 664 for the formal C64 sample)",
    )
    parser.add_argument(
        "--input-file",
        type=Path,
        help=(
            "required run-root formal SALMON input binding sysname, aligned "
            "atom geometry, and symmetry provenance"
        ),
    )
    parser.add_argument(
        "--symmetry-file",
        type=Path,
        help=(
            "optional cross-check of the input-named aligned symmetry path; "
            "relative paths are resolved against the recovered run root"
        ),
    )
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.log is None:
        parser.error("log is required unless --self-test is used")
    if args.input_file is None:
        parser.error("--input-file is required for formal aligned validation")
    if args.expected_operations <= 0 or args.expected_bands <= 0:
        parser.error("--expected-operations and --expected-bands must be positive")
    try:
        validate_integration(
            args.log.resolve(), args.expected_operations, args.expected_bands,
            input_file=args.input_file,
            symmetry_file=args.symmetry_file,
        )
    except ValidationError as exc:
        raise SystemExit(f"FAIL SAWF integration: {exc}") from exc


if __name__ == "__main__":
    main()
