#!/usr/bin/env python3
from pathlib import Path

root = Path(__file__).resolve().parents[2]
builder = (root / "cmakefiles/Builder/build_eigenexa.cmake").read_text().lower()
patcher = (root / "cmakefiles/Builder/patches/apply_eigenexa_2_4b_single_thread.cmake").read_text().lower()

assert "apply_eigenexa_2_4b_single_thread.cmake" in builder
assert "patch_command" in builder
assert "if ( local_size > 1 ) then" in patcher
assert "allocate(u0_z" in patcher
assert "deallocate(u0_z" in patcher
assert "failed to patch eigenexa single-thread initialization" in patcher

print("PASS EigenExa 2.4b single-thread work arrays are patched")
