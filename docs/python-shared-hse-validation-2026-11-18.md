# Python shared-model HSE verification — 2026-11-18

NativeFastHSEBackend.open_model returns one model/verified work copy/host-MPI
attachment. Explicit state leases dispatch through generation-checked tokens.
The first state initializes native inputs. A separate model allocator supports
creation after zero live states without retaining a hidden full numerical state.
Legacy single-state creation remains unchanged.

The local numerical campaign uses synthetic8³/full8-k original-GS fixtures and
learned Taylor/learned PT-CN. Each runs internal MPI1/OMP1, MPI2/OMP1 and MPI2/OMP2.
Five jobs per case (isolated A, isolated B, shared A/B, preceding singleton A/B)
use two independent solver groups. All current/energy arrays agree exactly in
six cases/30 MPI invocations. A uses dt=.001/capacity32 and z drive; B uses
dt=.002/capacity34 and rotating xy drive.

Shared cases also save atstep8, restore into fresh compatible A/B states and
compare every future current/energy to continuous execution, with no propagation
replay. Original states remain unchanged. A model closed with one state alive
retains its attachment through that state's final step; final close retires the
parent dependency. Isolated controls close all numerical states and recreate
fresh original-GS states twice in the still-open model.

Additional MPI2 probes use public native_step_targets to roll back one state
without truncating the sibling's flushed stream, and simultaneous native step
backups. A truncated fresh-state checkpoint fails and is closed; siblings and
subsequent correct fresh restores remain usable. The actual shared Python example
is executed separately. Five lightweight unit checks cover entry points, leases,
empty-model reservation release, bound ctypes signatures and journal scoping.

No numerical propagation formula or checkpoint encoding changed. The new native
model allocator shares the same checked constructor as live-token creation.
The attached model owns child communicators until final detach; this does not
complete arbitrary heterogeneous/symmetry reentrancy, unlimited object cycling
or automatic multipoint MPI-plan integration. Physical convergence is separate.
Actions are not used. Source/installed tests and sdist-to-wheel checks accompany
this source-layout change; raw archives stay private local evidence.
