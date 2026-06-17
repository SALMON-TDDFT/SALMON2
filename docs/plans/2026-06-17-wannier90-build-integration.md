# Wannier90 Build Integration Implementation Plan

**Goal:** Add a build-time Wannier90 dependency path so SALMON can fetch, build, and link a Wannier90 static library.

**Architecture:** Reuse SALMON's existing ExternalProject pattern. LAPACK/ScaLAPACK/EigenExa are resolved first, then Wannier90 is built as a serial static library and appended to `EXTERNAL_LIBS`. Fortran runtime use is intentionally stubbed for now; this task only establishes the dependency boundary.

**Tech Stack:** CMake ExternalProject, configure.py options, Fortran/C preprocessor config.

---

### Task 1: Build Option Wiring

**Files:**
- Modify: `CMakeLists.txt`
- Modify: `configure.py`
- Modify: `src/config.h.in`

**Steps:**
1. Add `USE_WANNIER90` beside other third-party library options.
2. Add `--enable-wannier90` and `--disable-wannier90`.
3. Emit `#cmakedefine USE_WANNIER90`.

### Task 2: ExternalProject Builder

**Files:**
- Create: `cmakefiles/Builder/build_wannier90.cmake`
- Modify: `cmakefiles/build_required_packages.cmake`

**Steps:**
1. Download Wannier90 release tarball.
2. Generate a minimal `make.inc` using SALMON's Fortran compiler and flags.
3. Build `libwannier.a`.
4. Import the library target and append it to `EXTERNAL_LIBS`.

### Task 3: Configure Smoke Test

**Commands:**
- `cmake -S . -B /tmp/salmon_wannier90_config -DUSE_WANNIER90=ON`

**Expected:** CMake configures and creates a `wannier90` imported target. Full build may require network access to fetch Wannier90.
