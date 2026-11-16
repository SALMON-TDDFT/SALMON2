# Direct checkpoint coverage beyond PZ

User authorized extension after committing1851c533. Preserve same-topology fresh-job direct restore; no propagation replay.

Supported initial extension: builtin TB-mBJ (automatic and fixed cval), Libxc-PZ and Libxc-PBE. Builtin `pbe` is unavailable in the solver; HSE requires a separate resource serialization effort and is not admitted here.

Audit: s_xc_functional contains immutable configuration, cval and optional opaque Libxc handles. TB-mBJ automatic c is evaluation-local. Checkpoint V2 records XC/x/c names, xctype/ispin, derivative flags, input cval and current scalar cval. Validate against newly initialized configuration, restore scalar only after successful collective payload read, retain fresh Libxc handles. Preserve all existing orbital/potential/Hartree/RT/Poisson history serialization.

Validation: build with explicit supported Libxc Fortran90 installation in a new external directory. Prepare separate Si SCF inputs per functional and use the same XC in RT. Compare continuous128 with checkpoint63→128, including fixed/automatic TB-mBJ; odd split exercises buffer parity. Track restored step and propagation calls. PZ remains a regression. Validate source manifests, packaging, tests and synchronize English/Japanese/AI guides. These are numerical continuation fixtures, not optical convergence or physical accuracy validation.
