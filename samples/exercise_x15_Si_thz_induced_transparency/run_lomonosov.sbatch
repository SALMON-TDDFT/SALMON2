#!/bin/bash
#SBATCH --job-name=x15_si_thz
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#
# One run directory per job: submit with the directory as the only argument,
#     for d in runs/*/; do sbatch run_lomonosov.sbatch "$d"; done
# Each directory is self-contained, so the jobs cannot clobber one another.
set -euo pipefail
RUNDIR="${1:?usage: sbatch run_lomonosov.sbatch <run-directory>}"
SALMON="${SALMON:-salmon}"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"
cd "$RUNDIR"
[ -f Si_prim_tm.data ] || "$SALMON" < Si_prim_epm_gs.inp > gs.log 2>&1

# Resume from the B4 checkpoint if this job is a requeue: without the resume input
# the solver reopens the outputs with status='replace' and silently restarts from
# t = 0, so a wall-clock kill would cost the whole run rather than <= 200 steps.
INP=Si_prim_sbe_rt.inp
if [ -s Si_prim_sbe_rt.data ] && ls ./*ckpt* >/dev/null 2>&1; then
    INP=Si_prim_sbe_rt_resume.inp
fi
"$SALMON" < "$INP" >> rt.log 2>&1
