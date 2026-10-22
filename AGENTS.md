# Mandatory eigensolver policy

PZHEEV is prohibited by the human user's instruction (2026-10-22).

- Do not execute PZHEEV, including workspace queries, reproductions, diagnostics,
  optional branches, or fallback paths. PZHEEVD is a different, permitted driver;
  replacement implementations still require numerical and capacity validation.
- Before building any SALMON source tree, including a historical snapshot or
  extracted archive, run `python3 tools/check_forbidden_eigensolvers.py
  --source-root <that-source-tree>/src` using this workspace's checker. Stop on
  prohibited references, unavailable cloud-only files, or inspection errors.
- Before transferring or submitting a newly built executable, inspect all its
  application objects with the same checker's `--object-tree <object-directory>`.
  For ad-hoc link commands, inspect every linked application object, including
  objects outside the main build directory. Never infer safety from the filename
  or binary hash alone. Vendor libraries exporting unused PZHEEV symbols do not
  by themselves prove the application calls it.
- CMake configure and incremental builds must retain the prohibition checks.
  Never disable or bypass them to build an old source tree.
- Preserve historical source, logs, archives and failed binaries as evidence;
  they are not approved for execution or reuse. Do not mass-rewrite snapshots.
- Do not modify, replace, cancel or resubmit existing remote production jobs or
  fixed binaries without the user's applicable authorization. Audit their actual
  execution paths and report any uncertainty instead of assuming they are safe.
