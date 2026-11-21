# Optional final observation-info I/O

Bounded R3 continuation after point/integral optional-output handling.

- Check obs0_info.data open/write/close IOSTAT. A failure warns and proceeds
  through normal field/state/MPI cleanup instead of invalidating completed
  physical evolution. Close only owned units and retain the first I/O error.
- Preserve formatting,calculation order and existing result/checkpoint files.
  No extra collective/global state/numerical work or post-finalization retry.
  Missing/incomplete metadata must be checked before visualization/analysis.
- Required restart saves and other output families remain separate policies.
  Keep all added guards inside standard132-column source lines.

Record:`salmonpy-fdtd-info-output-recovery-20261121/summary.json` in canonical
calc2026. Prior runtime fails a real directory collision at finalization. New
source/installed GNU Mac MPI1/2:10 directory/read-only/caller-file faults finish
normally;completed trajectories,integrals and checkpoints stay byte-identical.
Normal scalar/tensor metadata/output bytes match prior controls. Zero-replay
continuation,state/MPI cleanup and caller-file lifetime pass. Both suites679
run,one optional skip,no failures. Final source/object audits,Git-only build
and sdist→wheel pass. No GS collection or optional CI dispatch.

Close-only failure was not injected. This is not a blanket filesystem-failure
certificate or permission to silently accept failed required restart writes.
