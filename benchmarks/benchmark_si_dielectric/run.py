"""Sequential Si k4 impulse-only GS/RT comparison; no background subtraction."""
import argparse
import hashlib
import json
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]
FUNCTIONALS = ('pbe', 'hse06', 'pbe0', 'pbeh40')


def make_input(xc, mode, producer, dt=.08, nt=4375):
    if xc not in FUNCTIONALS or mode not in ('gs', 'rt'):
        raise ValueError('unknown functional or calculation')
    gs = mode == 'gs'
    control = "write_gs_restart_data='wfn'" if gs else f"directory_read_data='{producer}/data_for_restart/'"
    extra = ''
    if xc != 'pbe':
        extra = "exx_mlwf_interval=5\n exx_mlwf_maxiter=100\n exx_mlwf_tolerance=1d-7\n yn_hse_profile='y'"
        if gs:
            extra += '\n exx_pre_scf_threshold=1d-4'
    atoms = [(0,0,0),(.25,.25,.25),(.5,0,.5),(0,.5,.5),(.5,.5,0),(.75,.25,.75),(.25,.75,.75),(.75,.75,.25)]
    coords = '\n'.join(" 'Si' " + ' '.join(f'{10.26*x:.10f}' for x in a) + ' 1' for a in atoms)
    text = f"""&calculation
 theory='{'dft' if gs else 'tddft_response'}'
 yn_dc='n'
/
&control
 sysname='Si'
 yn_restart='n'
 {control}
/
&units
 unit_system='a.u.'
/
&system
 yn_periodic='y'
 al=10.26d0,10.26d0,10.26d0
 nstate=16
 nelec=32
 natom=8
 nelem=1
/
&pseudo
 izatom(1)=14
 file_pseudo(1)='Si_rps.dat'
 lmax_ps(1)=2
 lloc_ps(1)=2
/
&functional
 xc='{'libxc_pbe' if xc == 'pbe' else xc}'
 {extra}
/
&rgrid
 num_rgrid=16,16,16
/
&kgrid
 num_kgrid=4,4,4
/
&parallel
 nproc_k=4
 nproc_ob=1
 nproc_rgrid=1,1,1
/
&scf
 nscf=600
 ncg=4
 threshold=1d-8
 nscf_init_no_diagonal=0
 nscf_init_redistribution=0
/
&atomic_coor
{coords}
/
"""
    if not gs:
        text += f"""&tgrid
 dt={dt}
 nt={nt}
/
&emfield
 trans_longi='tr'
 ae_shape1='impulse'
 e_impulse=1d-4
 epdir_re1=1,0,0
/
&analysis
 out_rt_energy_step=1
 nenergy=800
 de=0.0018374661087827496
/
"""
    return text


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    temporary = path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n')
    temporary.replace(path)


def execute(root, exe, mpi, name, inp, steps=None):
    folder = root/name
    folder.mkdir(exist_ok=True)
    status_file = folder/'status.json'
    if status_file.exists():
        previous = json.loads(status_file.read_text())
        if previous.get('complete') and (folder/'inputfile').read_text() == inp and previous['binary_sha256'] == digest(exe):
            return previous
        raise RuntimeError(f'{name}: existing result/input differs or is incomplete; preserve it before retry')
    if (folder/'output').exists():
        raise RuntimeError(f'{name}: existing output without completion metadata')
    (folder/'inputfile').write_text(inp)
    shutil.copy2(ROOT/'testsuites/pseudo/Si_rps.dat', folder/'Si_rps.dat')
    command = [mpi, '--bind-to', 'none', '-n', '4', sys.executable,
               str(Path(__file__).with_name('rank_run.py')), str(exe)]
    row = dict(name=name, complete=False, command=command, started=time.time(),
               binary_sha256=digest(exe), input_sha256=digest(folder/'inputfile'))
    save(status_file, row)
    save(root/'active.json', row)
    print('START', name, flush=True)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    with (folder/'inputfile').open('rb') as fi, (folder/'output').open('wb') as fo:
        process = subprocess.Popen(command, cwd=folder, stdin=fi, stdout=fo, stderr=subprocess.STDOUT, env=env)
        row['launcher_pid'] = process.pid
        save(root/'active.json', row)
        row['returncode'] = process.wait()
    row['wall_seconds'] = time.time()-row['started']
    log = (folder/'output').read_text(errors='replace')
    row['complete'] = row['returncode'] == 0 and 'end SALMON' in log
    if steps is None:
        row['complete'] &= '#GS converged at' in log
    else:
        import numpy as np
        data = np.loadtxt(folder/'Si_rt.data', ndmin=2)
        row['complete'] &= len(data) == steps and bool(np.isfinite(data).all())
        row['wall_seconds_per_step_including_setup'] = row['wall_seconds']/steps
    ranks = [json.loads(p.read_text()) for p in folder.glob('rank-*.json')]
    if len(ranks) != 4:
        row['complete'] = False
    row['rank_peak_rss_bytes'] = [p['peak_rss_bytes'] for p in ranks]
    row['sum_rank_peaks_bytes_not_simultaneous'] = sum(row['rank_peak_rss_bytes'])
    save(status_file, row)
    save(root/'active.json', row)
    print('END', name, row['wall_seconds'], row['complete'], flush=True)
    if not row['complete']:
        raise RuntimeError(f'{name}: failed or unconverged; see output')
    return row


def check_dt(root, xc):
    import numpy as np
    coarse = np.loadtxt(root/f'{xc}-probe-008/Si_rt.data')
    fine = np.loadtxt(root/f'{xc}-probe-004/Si_rt.data')[1::2]
    np.testing.assert_allclose(coarse[:,0], fine[:,0], rtol=0, atol=1e-9)
    error = float(np.max(abs(coarse[:,13:16]-fine[:,13:16])))
    signal = float(np.max(abs(fine[:,13:16])))
    energy_drifts = []
    norm_errors = []
    for suffix in ('008', '004'):
        folder = root/f'{xc}-probe-{suffix}'
        energy = np.loadtxt(folder/'Si_rt_energy.data', ndmin=2)
        if not np.isfinite(energy).all():
            raise RuntimeError('nonfinite RT energy')
        energy_drifts.append(float(np.max(abs(energy[:,1]-energy[0,1])))/8)
        log = (folder/'output').read_text(errors='replace')
        pattern = r'^\s*\d+\s+[\d.]+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+\S+\s*$'
        norms = [float(v) for v in re.findall(pattern, log, re.M)]
        if not norms:
            raise RuntimeError('missing RT electron norm')
        norm_errors.append(max(abs(v-32) for v in norms))
    row = dict(current_max_difference_au=error, current_max_au=signal,
               energy_max_drift_ha_per_atom=energy_drifts, electron_norm_error=norm_errors,
               accepted=(error <= max(1e-9, .01*signal) and max(energy_drifts)<1e-5
                         and max(norm_errors)<1e-5))
    save(root/f'{xc}-dt-check.json', row)
    if not row['accepted']:
        raise RuntimeError(f'{xc}: refine time step before production')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--exe', type=Path, required=True)
    parser.add_argument('--mpi', default='/opt/homebrew/bin/mpiexec')
    parser.add_argument('--phase', choices=['prepare', 'preflight', 'production', 'all'], default='preflight')
    parser.add_argument('--functionals', nargs='+', choices=FUNCTIONALS, default=list(FUNCTIONALS))
    a = parser.parse_args()
    root = a.root.resolve();root.mkdir(parents=True, exist_ok=True)
    exe = root/'salmon'
    if not exe.exists():
        shutil.copy2(a.exe.resolve(), exe)
    elif digest(exe) != digest(a.exe):
        raise RuntimeError('executable snapshot differs; choose a new comparison directory')
    if a.phase == 'prepare':
        for xc in a.functionals:
            producer = root/f'{xc}-gs'
            for suffix, mode, dt, nt in [('gs','gs',.08,4375), ('probe-008','rt',.08,16),
                                         ('probe-004','rt',.04,32), ('impulse','rt',.08,4375)]:
                folder = root/f'{xc}-{suffix}';folder.mkdir(exist_ok=True)
                if (folder/'output').exists():
                    raise RuntimeError('refusing to rewrite a run that has started')
                (folder/'inputfile').write_text(make_input(xc, mode, producer, dt, nt))
        return
    if a.phase in ('preflight', 'all'):
        for xc in a.functionals:
            gs = root/f'{xc}-gs'
            execute(root, exe, a.mpi, f'{xc}-gs', make_input(xc, 'gs', gs))
            for dt, steps, suffix in [(.08,16,'008'), (.04,32,'004')]:
                execute(root, exe, a.mpi, f'{xc}-probe-{suffix}', make_input(xc, 'rt', gs, dt, steps), steps)
            check_dt(root, xc)
    if a.phase in ('production', 'all'):
        for xc in a.functionals:
            if not json.loads((root/f'{xc}-dt-check.json').read_text())['accepted']:
                raise RuntimeError('preflight required')
            gs = root/f'{xc}-gs'
            execute(root, exe, a.mpi, f'{xc}-impulse', make_input(xc, 'rt', gs), 4375)
    save(root/'queue-status.json', dict(complete=True, phase=a.phase, functionals=a.functionals))


if __name__ == '__main__':
    main()
