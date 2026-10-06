"""Matched fixed-work Si HSE06 CPU/cuFFT trials; no convergence/performance claims from a pilot."""
import argparse
import json
import os
from pathlib import Path
import re
import shlex
import signal
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--exe', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--nk', type=int, default=2, help='k points on each axis; use 8 for target benchmark')
    p.add_argument('--grid', type=int, default=12, help='real-space mesh points on each axis')
    p.add_argument('--steps', type=int, default=3)
    p.add_argument('--ranks', type=int, default=1)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--block-rows', type=int, default=8)
    p.add_argument('--repeats', type=int, default=1)
    p.add_argument('--energy-tolerance', type=float, default=1e-6, help='CPU/GPU SCF energy-trajectory tolerance in eV')
    p.add_argument('--backends', nargs='+', choices=('cpu', 'cufft'), default=['cpu', 'cufft'])
    p.add_argument('--mpi-args', default='', help='extra site-specific mpiexec options, shell-quoted as one argument')
    p.add_argument('--timeout', type=int, default=3300, help='per trial seconds; PBS walltime is the job-wide limit')
    a = p.parse_args()
    if min(a.nk, a.grid) < 2 or min(a.steps, a.ranks, a.threads, a.block_rows, a.repeats, a.timeout) < 1:
        p.error('positive controls and at least two mesh points per axis required')
    if not (0 < a.energy_tolerance < float('inf')):
        p.error('finite positive energy tolerance required')
    if a.ranks > a.nk**3:
        p.error('each rank must own at least one k point')
    exe = a.exe.expanduser().resolve()
    if not exe.is_file():
        p.error('SALMON executable missing')
    out = a.output.expanduser().resolve()
    out.mkdir(parents=True, exist_ok=False)
    template = (ROOT / 'testsuites/420_bulk_Si_hse_gs/inputfile').read_text()
    template = re.sub(r'^!.*\n', '', template, count=1)
    template = template.replace('num_kgrid = 4, 4, 4', f'num_kgrid = {a.nk}, {a.nk}, {a.nk}')
    template = template.replace('dl = 0.855d0, 0.855d0, 0.855d0', f'num_rgrid = {a.grid}, {a.grid}, {a.grid}')
    template = template.replace('nscf = 160', f'nscf = {a.steps}').replace('threshold = 1.d-8', 'threshold = 1.d-30')
    template = template.replace('nproc_k=4', f'nproc_k={a.ranks}')
    controls = f"\n yn_hse_wannier='n'\n yn_hse_profile='y'\n hse_block_rows={a.block_rows}\n"
    template = template.replace("xc = 'hse06'", "xc = 'hse06'" + controls)
    metadata = vars(a).copy()
    metadata['exe'] = str(exe); metadata['output'] = str(out)
    metadata['note'] = 'Fixed SCF work, not a converged material result. Wall time includes initialization.'
    (out / 'settings.json').write_text(json.dumps(metadata, indent=2) + '\n')
    records = []
    for repeat in range(1, a.repeats + 1):
        order = a.backends if repeat % 2 else list(reversed(a.backends))
        for backend in order:
            folder = out / f'{backend}-{repeat}'
            folder.mkdir()
            inp = template.replace("xc = 'hse06'", f"xc = 'hse06'\n exx_kpoint_backend='{backend}'")
            (folder / 'inputfile').write_text(inp)
            shutil.copyfile(ROOT / 'testsuites/pseudo/Si_rps.dat', folder / 'Si_rps.dat')
            env = dict(os.environ, SALMON_EXE=str(exe), OMP_NUM_THREADS=str(a.threads),
                       OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1')
            env.pop('NVCOMPILER_ACC_NOTIFY', None)
            env.pop('NVCOMPILER_ACC_TIME', None)
            command = [os.environ.get('MPIEXEC', 'mpiexec'), *shlex.split(a.mpi_args),
                       '-n', str(a.ranks), 'bash', str(HERE / 'rank_run.sh')]
            start = time.monotonic()
            with (folder / 'output.log').open('w') as log:
                try:
                    process = subprocess.Popen(command, cwd=folder, env=env, stdout=log,
                                               stderr=subprocess.STDOUT, start_new_session=True)
                    returncode = process.wait(timeout=a.timeout)
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid, signal.SIGTERM)
                    try:
                        process.wait(timeout=10)
                    except subprocess.TimeoutExpired:
                        os.killpg(process.pid, signal.SIGKILL)
                        process.wait()
                    (folder / 'FAILED').write_text('Trial timeout; launcher process group terminated.\n')
                    raise
            elapsed = time.monotonic() - start
            text = (folder / 'output.log').read_text(errors='replace')
            record = dict(backend=backend, repeat=repeat, wall_seconds=elapsed, returncode=returncode)
            record['gpu_dispatch_seen'] = 'EXX_KPOINT_BACKEND=cufft' in text
            record['profile_rows'] = re.findall(r'HSE_PROFILE rank0 block=.*', text)
            record['iteration_energy_ev'] = [float(v) for v in re.findall(
                r'iter=\s*\d+\s+Total Energy=\s*([\d.Ee+\-]+)', text)]
            timing = re.search(r'^scf iterations\s+(\S+)\s+(\S+)', text, re.M)
            if timing:
                record['scf_wall_max_seconds'] = float(timing[2])
                record['scf_seconds_per_iteration'] = float(timing[2]) / a.steps
            total = re.search(r'total calculation time,\s*([\d.Ee+\-]+)', text)
            if total:
                record['salmon_total_seconds'] = float(total[1])
            infos = sorted(folder.glob('**/*_info.data'))
            if infos:
                match = re.search(r'Total energy \(eV\) =\s*([\d.Ee+\-]+)', infos[-1].read_text())
                if match:
                    record['energy_ev'] = float(match[1])
            records.append(record)
            (out / 'summary.json').write_text(json.dumps(records, indent=2) + '\n')
            if returncode or 'end SALMON' not in text:
                raise RuntimeError(f'{folder}: SALMON failed; inspect output.log')
            if len(record['iteration_energy_ev']) != a.steps:
                raise RuntimeError('SCF iteration count differs from requested fixed work')
            if backend == 'cufft' and not record['gpu_dispatch_seen']:
                raise RuntimeError('GPU route was not exercised')
            print(f'{backend} repeat={repeat} total wall={elapsed:.3f}s', flush=True)
    comparisons = []
    for repeat in range(1, a.repeats + 1):
        pair = {r['backend']: r for r in records if r['repeat'] == repeat}
        if set(pair) == {'cpu', 'cufft'}:
            delta = max(abs(x-y) for x, y in zip(pair['cpu']['iteration_energy_ev'],
                                                pair['cufft']['iteration_energy_ev']))
            comparisons.append(dict(repeat=repeat, max_energy_difference_ev=delta,
                                    energy_match=delta <= a.energy_tolerance,
                                    wall_speedup=pair['cpu']['wall_seconds']/pair['cufft']['wall_seconds']))
    (out / 'comparison.json').write_text(json.dumps(comparisons, indent=2) + '\n')
    if any(not c['energy_match'] for c in comparisons):
        raise RuntimeError('CPU/GPU energy trajectory mismatch: inspect comparison.json')
    print(f'Results: {out / "summary.json"}; rank RSS and GPU memory samples remain in each trial directory.')


if __name__ == '__main__':
    main()
