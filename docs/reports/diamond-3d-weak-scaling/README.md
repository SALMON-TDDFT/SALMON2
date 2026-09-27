# Diamond 3D DC-HSE / real-space RT inputs

These active inputs now propagate the wavefunctions on the full real-space mesh
with SALMON Taylor4. DC/LCFO supplies the initial orbitals only. There is no
fixed LCFO basis projection during RT, so the initial LCFO basis does not limit
accessible excited states. This requires the new real-space HSE RT implementation;
the older namelist-only and radius patches are insufficient.

**All earlier LCFO RT timings, memory estimates and weak-scaling results are
historical fixed-basis results. They do not validate this route or its spectra.**
This implementation prioritizes correctness. Full-domain exchange FFTs and
collectives remain, and no new large-system performance claim is made.

| Cells | Atoms | MPI ranks | RT occupied orbitals |
|---|---:|---:|---:|
| 4×4×4 | 512 | 64 | 1024 |
| 6×6×6 | 1728 | 216 | 3456 |
| 8×8×8 | 4096 | 512 | 8192 |
| 10×10×10 | 8000 | 1000 | 16000 |

Core mesh is 16³ per cell; GS fragment mesh including buffers is 32³.
GS inputs are unchanged physically. RT uses dt=0.02 a.u., 16 steps and an x impulse.
Each RT directory reads ../gs/data_dcdft. Check GS/LCFO convergence first.

All algorithm controls are in &functional:

```fortran
 yn_hse_realspace_rt='y'
 yn_hse_wannier='y'
 hse_rt_wf_radius=6d0
 hse_rt_ace_interval=1
 hse_rt_u_interval=1
 hse_rt_fft_batch=1
 yn_hse_rt_fft_measure='n'
 yn_hse_rt_seed_distributed='y'
```

The WF radius is in bohr. Zero means full support; positive values mask exchange
sources and are an additional approximation. 6 bohr is a comparison setting,
not a certified 99.9% radius. Initial coverage below 99.9% produces a warning.
The radius does not truncate the propagated mesh wavefunction. MLWF U transports
only the occupied source representation. ACE is stored on mesh rows.
The old yn_hse_lcfo_rt / yn_hse_lcfo_direct_wf flags must not be enabled.
No SALMON_LCFO environment settings are used; gs-env.sh and rt-env.sh are inert.

Generate and statically validate (no numerical jobs):

```sh
python3 generate.py
python3 validate.py
```

Run GS and RT sequentially, e.g. use 64 MPI ranks for 4×4×4. On Fugaku, use
the allocation's launcher and set OMP_NUM_THREADS to the cores per rank.
Do not assume old fixed-basis memory estimates bound new real-space RT memory.
The manifest retains old LCFO memory quantities only as historical reference.
Bundled older patches are historical; build the complete real-space RT revision.

The matching `patches/realspace-hse-rt.patch` plus before/after SHA256 lists are
bundled here too. They require the namelist-only source (`9d2ebc27`, compatible
with the separate Tofu fix). Use the same hash→dry-run→apply→hash→CMake procedure
as [Si's input guide](../si-3d-weak-scaling/README.md). Fugaku is untested.
