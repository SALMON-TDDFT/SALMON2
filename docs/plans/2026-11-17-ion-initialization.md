# Owned ion initialization inputs

First slice of the approved remaining GS/TDDFT initialization ownership work.
Standard and fast-PR initialization_sub receive the same value-owned
s_ion_initial_config (nion,nelem,reduced,positions,species). Capture is the only
helper that reads legacy globals. initialize_ions validates the complete input,
resolves reduced coordinates against the selected system lattice and owns ion
positions/species/mass/velocity/force arrays; rejected inputs leave the target
unchanged. No FFT/MPI resources or whole solver contexts are copied.

The existing init_dft_system captures once, uses the helper and preserves its
legacy reduced-coordinate Rion normalization for downstream consumers. The
explicit helper never reads/writes global inputs. This does not claim complete
independent GS or remove global cell/k/symmetry/communicator dependencies.

Verify two different captured inputs after poisoning globals, alternating existing
systems, exact Cartesian/reduced references and no cross-state mutations; invalid
sizes/species/nonfinite input must preserve originals. Compile fail-first against
old module, then run standard and fast-PR helper probes with OMP1/2. Build new
external libraries from verified sources and compare ordinary mixed old/new
results (PBE and fast-PR Taylor/PTCN), including protected retry and direct restart
if feasible. Preserve manifests and provenance, test Python and distribution
contents, update paired manuals/AI JSON/checklist, retain local evidence. No Actions.
