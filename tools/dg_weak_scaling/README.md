# RT-DG weak-scaling check

Use this directory to collect comparable weak-scaling logs for the RT-DG path.

## What to keep fixed

For a weak-scaling check, keep the local work per fragment/subgroup fixed:

- same fragment size and local real-space grid per fragment
- same `nstate_frag`
- same `MPI ranks per fragment subgroup`
- same RT step count `nt`
- same timing environment variables

Scale these together:

- total fragments
- total atoms/electrons
- total MPI ranks
- total cell length/grid along the replicated direction

For the current DG fragment orbital mode, a useful first series is:

| label | fragments | MPI ranks | subgroup ranks |
| --- | ---: | ---: | ---: |
| ws_64 | 64 | 256 | 4 |
| ws_128 | 128 | 512 | 4 |
| ws_256 | 256 | 1024 | 4 |
| ws_512 | 512 | 2048 | 4 |

## Recommended timing flags

Set these in the job script:

```sh
export OMP_NUM_THREADS=1
export SALMON_DG_RK_TIMING=1
export SALMON_DG_STAGE_TIMING=1
export SALMON_DG_DENSITY_TIMING=1
export SALMON_DG_DERIVATIVE_TRACE=1
export SALMON_DG_INITIAL_EIGENSTATE_CHECK=1
```

If the purpose is only timing, keep `nt` modest, for example 10-20, and ignore
the first RT step when reading manually because it includes extra cache setup.

## Parse logs

After runs complete, put each `run.log` under a case directory:

```text
results/
  ws_64/run.log
  ws_128/run.log
  ws_256/run.log
  ws_512/run.log
```

Then run:

```sh
python tools/dg_weak_scaling/parse_dg_weak_scaling.py --format markdown results
python tools/dg_weak_scaling/parse_dg_weak_scaling.py results > dg_weak_scaling.csv
```

Important columns:

- `rk_mean_s`: mean RK step time if `rk timing` is present
- `stage_mean_s`: mean density/H update time from `[DG-RK] stage ... density/H update`
- `deriv_mean_s`: mean derivative time from `[DG-RK] stage ... derivative`
- `full_rel_final`: stationary check after the distributed DG eigensolver
- `dist_eig_solved`: whether the full distributed initial-state solve succeeded

Good weak scaling means `rk_mean_s` stays roughly flat as `mpi_ranks` and
`n_frag` grow together. If `rk_mean_s` grows, compare `stage_mean_s` and
`deriv_mean_s` first; that separates density/H reconstruction from coefficient
propagation and overlap solve costs.
