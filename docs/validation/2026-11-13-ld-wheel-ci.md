# LD wheel CI verification

Run: https://github.com/SALMON-TDDFT/SALMON2/actions/runs/43354726041
Code: 887391751b6918c7edc6ac819e49fc07dc82202a
Branch: codex/salmonpy-ld-wheels-20261113

Windows AMD64, Linux x86_64, macOS Intel and macOS arm64 each built and repaired
CPython 3.10, 3.11, 3.12, 3.13 and 3.14 wheels. Each installed wheel passed 13
numerical/input-lifecycle tests (260 checks across 20 wheels). Linux sdist rebuild
and packaging hash/RECORD checks passed. Artifacts are attached to the run.

Failures fixed during actual execution:
- macos-14 runner retired: use macos-15.
- Homebrew gfortran runtimes required macOS 15: raise deployment target to 15.
- manylinux2014 could not install current NumPy wheels for newer Python, and
  fallback source build failed: use manylinux_2_28 (glibc 2.28+ usage baseline).

No merge or PyPI publication. CI code branch includes the isolated LD packaging
and permitted eigensolver replacement; existing FDTD extraction remains in the
working checkout and is not claimed to be covered by these wheel tests.
