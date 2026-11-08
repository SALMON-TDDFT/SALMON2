#!/bin/bash
# Run inside an allocated 1-node, 4-rank, 12-thread job on Fugaku.
# The executable's application-object audit must have passed before deployment.
set -euo pipefail
root=${1:?usage: run-stage.sh WORKDIR EXECUTABLE EXPECTED_SHA256 STAGE}
exe=${2:?missing executable};expected=${3:?missing recorded SHA256};stage=${4:?missing stage}
case "$stage" in gs|updated|snapshot|restart|external) ;; *) echo 'invalid stage' >&2; exit 2;; esac
if command -v sha256sum >/dev/null; then actual=$(sha256sum "$exe"); else actual=$(shasum -a 256 "$exe"); fi
actual=${actual%% *}
if test "$actual" != "$expected"; then echo 'Executable differs from audited SHA256' >&2; exit 2; fi
export OMP_NUM_THREADS=12 PARALLEL=12 OMP_DYNAMIC=FALSE OMP_THREAD_LIMIT=12 FLIB_FASTOMP=FALSE
# Four ranks per node use the established Fugaku software-barrier policy.
export FLIB_NOHARDBARRIER=TRUE FLIB_CNTL_BARRIER_ERR=FALSE
# Sequential vendor BLAS inside q OpenMP avoids nested oversubscription.
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
cd "$root/$stage"
for previous in run.log launcher.log launcher.stderr stdout.rank-* stderr.rank-*; do
  if test -e "$previous"; then echo "$previous exists: preserve evidence; choose a new work directory" >&2; exit 2; fi
done
if test "$stage" != gs; then
  cp "$root/gs/hsex-w.bin" hsex-w.bin
fi
if test "$stage" = restart; then
  test -f "$root/updated/checkpoint_rt_000001/hsex-active.bin"
  ln -s "$root/gs/data_for_restart" restart
fi
if test "$stage" = external; then
  cp "$root/reference/pump.bin" python-pump.bin
fi
# Redirect on each rank: Fujitsu mpiexec need not forward launcher stdin.
status=0
mpiexec -n 4 bash -c '
set -eu
cd "$1"
rank=${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-${PMI_RANK:-}}}
if test -z "$rank"; then echo "MPI rank environment unavailable" >&2; exit 2; fi
exec "$2" < input > "stdout.rank-$rank" 2> "stderr.rank-$rank"
' hsex-si "$PWD" "$exe" > launcher.log 2> launcher.stderr || status=$?
if test -f stdout.rank-0; then cat stdout.rank-0 > run.log; else : > run.log; fi
if test "$status" -ne 0; then
  echo "SALMON MPI launch failed ($status): retain launcher.log, launcher.stderr and rank logs" >&2
fi
exit "$status"
