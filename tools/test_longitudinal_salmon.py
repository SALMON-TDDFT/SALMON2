"""Compare salmonpy.maxwell field components + TDDFT worker with SALMON's own fields.

PYTHONPATH=. python tools/test_longitudinal_salmon.py SALMON_EXECUTABLE NEW_DIRECTORY [--trans-longi lo|2d|ms1d]

The same synthetic, normalized (not converged) Si fixture is propagated twice:
  1. ordinary SALMON with trans_longi='lo' or '2d', or the 1D multiscale solver with
     one or more macroscopic points ('ms1d'): the field is computed inside SALMON;
  2. the SALMONpy TDDFT worker (trans_longi='tr') coupled to salmonpy.maxwell.
Ac_tot and the matter current must agree to SALMON's printed precision.
Requires SALMONPY_RUNTIME_SOURCE (prepared runtime source, for Si_rps.dat).
"""
import argparse
import hashlib
import json
from contextlib import ExitStack
import os
import shutil
import struct
import subprocess
from pathlib import Path
import numpy as np
from salmonpy.fdtd import SalmonExecutable
from salmonpy.fields import AcosPulse
from salmonpy.maxwell import Longitudinal, ThinFilm, Multiscale1D, couple
from salmonpy.tddft import TDDFTInitialState
from salmonpy.analysis import read_table

DT, NT = 0.02, 400
PULSE = dict(amplitude=0.02, omega=0.5, duration=6.0, polarization=(0., 0., 1.))
FILM = dict(thickness=10.0, epsilon=(1.0, 2.25))
LAYER = dict(spacing=5.0, vacuum=(300, 50))  # c*dt = 2.74 < spacing: not the magic step


def fixture(bundle, runtime_source):
    restart = bundle/'data_for_restart'; restart.mkdir(parents=True)
    def record(stream, data):
        stream.write(struct.pack('=i', len(data))); stream.write(data); stream.write(struct.pack('=i', len(data)))
    with (restart/'info.bin').open('wb') as s:
        for value in (1, 2, 0, 1, 0): record(s, struct.pack('=i', value))
    with (restart/'occupation.bin').open('wb') as s: record(s, np.array([2., 2.], dtype='=f8').tobytes())
    x = np.arange(8)[:, None, None]; z = np.arange(8)[None, None, :]
    psi = np.empty((8, 8, 8, 2), dtype='=c16', order='F')
    psi[:, :, :, 0] = 1/np.sqrt(8.**3)
    psi[:, :, :, 1] = np.exp(2j*np.pi*(x+z)/8)/np.sqrt(8.**3)
    (restart/'wfn.bin').write_bytes(psi.tobytes(order='F'))
    shutil.copy2(Path(runtime_source)/'samples/exercise_06_bulkSi_rt/Si_rps.dat', bundle/'Si_rps.dat')


def inp(trans_longi, emfield_extra, maxwell='', theory='tddft_pulse'):
    return f"""&calculation
 theory='{theory}'
/
&control
 sysname='fixture', directory_read_data='./data_for_restart/'
/
&units
 unit_system='au'
/
&system
 yn_periodic='y', al=8,8,8, nelem=1, natom=1, nelec=4, nstate=2
/
&pseudo
 file_pseudo(1)='Si_rps.dat', izatom(1)=14, lloc_ps(1)=2
/
&functional
 xc='PZ'
/
&rgrid
 num_rgrid=8,8,8
/
&kgrid
 num_kgrid=1,1,1
/
&tgrid
 dt={DT}, nt={NT}
/
&propagation
 propagator='middlepoint', yn_predictor_corrector='n'
/
&emfield
 trans_longi='{trans_longi}'{emfield_extra}
/
{maxwell}&atomic_red_coor
 'Si' 0.0 0.0 0.0 1
/
"""


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('executable'); p.add_argument('directory')
    p.add_argument('--trans-longi', choices=('lo', '2d', 'ms1d'), default='lo')
    p.add_argument('--cells', type=int, default=1, help='material points for ms1d (default: 1)')
    p.add_argument('--mpiexec', default='mpiexec', help='MPI launcher for multi-point SALMON reference')
    args = p.parse_args()
    if args.cells < 1:
        p.error('--cells must be positive')
    if args.cells != 1 and args.trans_longi != 'ms1d':
        p.error('--cells applies only to ms1d')
    runtime_source = os.environ.get('SALMONPY_RUNTIME_SOURCE')
    if not runtime_source: raise SystemExit('set SALMONPY_RUNTIME_SOURCE')
    dest = Path(args.directory).resolve(); dest.mkdir(exist_ok=False)
    exe = Path(args.executable).resolve()
    sha = hashlib.sha256(exe.read_bytes()).hexdigest()
    mode = args.trans_longi

    # 1. SALMON's own induced field.
    native = dest/'salmon'; fixture(native, runtime_source)
    extra = (f", ae_shape1='Acos2', epdir_re1=0,0,1, E_amplitude1={PULSE['amplitude']}, "
             f"omega1={PULSE['omega']}, tw1={PULSE['duration']}, phi_cep1=0.0")
    maxwell = ''
    if mode == '2d':  # SALMON reads these from &maxwell
        maxwell = (f"&maxwell\n film_thickness={FILM['thickness']}, epsilon_em(1)={FILM['epsilon'][0]}, "
                   f"epsilon_em(2)={FILM['epsilon'][1]}\n/\n")
    if mode == 'ms1d':
        return multiscale(args, dest, exe, sha, runtime_source, extra)
    (native/'salmon.inp').write_text(inp(mode, extra, maxwell))
    with (native/'stdout.log').open('w') as log:
        subprocess.run([str(exe)], stdin=(native/'salmon.inp').open(), stdout=log, stderr=subprocess.STDOUT,
                       cwd=native, check=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
    table = read_table(native/'fixture_rt.data')
    salmon_A = np.column_stack([table.column(f'Ac_tot_{c}') for c in 'xyz'])
    salmon_J = np.column_stack([table.column(f'Jm_{c}') for c in 'xyz'])

    # 2. SALMONpy: transverse worker + field component, coupled through couple().
    py = dest/'salmonpy'; fixture(py, runtime_source)
    (py/'rt.inp').write_text(inp('tr', ''))
    pulse = AcosPulse(**PULSE)
    field = Longitudinal(pulse) if mode == 'lo' else ThinFilm(external=pulse, **FILM)
    with TDDFTInitialState.from_salmon(py, 'rt.inp', executable=SalmonExecutable(exe, sha),
                                       max_steps=NT, threads=1, timeout=600.) as initial:
        result = couple(initial, field, dt=DT, steps=NT)
    py_A = result.A[:, 0]; py_J = -result.current[:, 0]  # back to SALMON matter current

    # SALMON prints rows for itt = 1..nt; align by the time column.
    index = np.rint(table.column('Time')/DT).astype(int)
    if not np.allclose(index*DT, table.column('Time'), atol=1e-7): raise SystemExit('unexpected time grid')
    py_A = py_A[index]; py_J = py_J[index]
    scaleA = np.abs(salmon_A).max(); scaleJ = np.abs(salmon_J).max()
    errA = np.abs(py_A-salmon_A).max()/scaleA
    errJ = np.abs(py_J-salmon_J).max()/scaleJ
    external = pulse.sample(np.arange(NT+2)*DT).values[index]
    induced = np.abs(salmon_A - external).max()/scaleA
    print(f'mode={mode} rows={len(index)} (t={index[0]*DT:g}..{index[-1]*DT:g}) '
          f'max|Ac_tot|={scaleA:.3e} induced/total={induced:.3e}')
    print(f'relative max difference: Ac_tot {errA:.2e}, matter current {errJ:.2e}')
    if induced < 1e-3: raise SystemExit('induced field too small for a meaningful comparison')
    if errA > 1e-12 or errJ > 1e-12: raise SystemExit('FAIL: SALMON and SALMONpy differ')
    print('PASS')


def multiscale(args, dest, exe, sha, runtime_source, extra):
    native = dest/'salmon'
    maxwell = (f"&multiscale\n fdtddim='1d', nx_m={args.cells}, ny_m=1, nz_m=1, "
               f"hx_m={LAYER['spacing']}, hy_m={LAYER['spacing']}, hz_m={LAYER['spacing']}, "
               # SALMON's checker requires nxvacl_m/nxvacr_m: nxvac = (|nxvacl_m|+1, |nxvacr_m|)
               f"nxvacl_m={LAYER['vacuum'][0]-1}, nxvacr_m={LAYER['vacuum'][1]}\n/\n"
               "&maxwell\n boundary_em(1,1)='abc', boundary_em(1,2)='abc'\n/\n"
               "&analysis\n out_ms_step=1\n/\n")
    (native/'salmon.inp').write_text(inp('tr', extra, maxwell, theory='multi_scale_maxwell_tddft'))
    command = [str(exe)] if args.cells == 1 else [args.mpiexec, '-n', str(args.cells), str(exe)]
    with (native/'stdout.log').open('w') as log:
        subprocess.run(command, stdin=(native/'salmon.inp').open(), stdout=log, stderr=subprocess.STDOUT,
                       cwd=native, check=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
    # SALMON can return status zero after an input/MPI setup error. Validate output before proceeding.
    reference = read_ms1d_reference(native, cells=args.cells, steps=NT)
    py = dest/'salmonpy'; fixture(py, runtime_source)
    (py/'rt.inp').write_text(inp('tr', ''))
    pulse = AcosPulse(**PULSE)
    field = Multiscale1D(cells=args.cells, incident=pulse, **LAYER).new_state(dt=DT, steps=NT)
    grids = [field.grid_A.copy()]
    with TDDFTInitialState.from_salmon(py, 'rt.inp', executable=SalmonExecutable(exe, sha),
                                       max_steps=NT, threads=1, timeout=600.) as initial:
        with ExitStack() as stack:
            states = [stack.enter_context(initial.new_state(dt=DT, initial_A=A)) for A in field.A]
            j = np.array([state.current for state in states])
            charge = [j.copy()]
            for _ in range(NT):
                A = field.step(j)
                grids.append(field.grid_A.copy())
                j = np.array([state.step(A[i]) for i, state in enumerate(states)])
                charge.append(j.copy())
    grids, charge = np.asarray(grids), np.asarray(charge)
    np.savez(dest/'reference.npz', **reference)
    np.savez(dest/'python.npz', grid_A=grids, current=charge, dt=DT)
    report = compare_ms1d_reference(reference, grids, charge)
    report.update(executable_sha256=sha, fixture='synthetic normalized Si; not converged',
                  dt=DT, steps=NT, cells=args.cells)
    (dest/'comparison.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))
    if not report['passed']:
        raise SystemExit('FAIL: SALMON and SALMONpy differ')
    print('PASS')


def read_ms1d_reference(native, *, cells, steps, vacuum=LAYER['vacuum']):
    """Read full-grid A and delayed charge current from SALMON ms output.

    feedback_current[n] is j[n-1] for n >= 2; at n=0,1 SALMON writes zero.
    Arrays retain every material point, including points with identical material.
    """
    indices = np.arange(1-vacuum[0], cells+vacuum[1]+1)
    grids, feedback = [], []
    for n in range(steps+1):
        table = np.loadtxt(Path(native)/'fixture_RT_Ac'/f'fixture_Ac_{n:06d}.data', comments='#')
        if not np.array_equal(table[:, 0], indices):
            raise ValueError('unexpected multiscale grid')
        grids.append(table[:, 3:6])
        feedback.append(table[(indices >= 1) & (indices <= cells), 12:15])
    return dict(grid_indices=indices, grid_A=np.asarray(grids),
                feedback_current=np.asarray(feedback))


def compare_ms1d_reference(reference, grid_A, current, *, tolerance=1e-12):
    """Compare endpoint histories with SALMON's one-step-lagged output.

    current[n, point, xyz] must be physical charge j at n*dt.  Initial j is
    intentionally omitted by SALMON's field update, and the final endpoint j
    is not present in its final field dump. Check these alignments explicitly.
    """
    native_A = np.asarray(reference['grid_A'])
    feedback = np.asarray(reference['feedback_current'])
    grid_A, current = np.asarray(grid_A), np.asarray(current)
    if grid_A.shape != native_A.shape or current.shape != feedback.shape:
        raise ValueError('history shape differs from SALMON reference')
    if len(grid_A) < 3:
        raise ValueError('at least two steps are required to compare current feedback')
    scale_A = float(np.max(np.abs(native_A)))
    scale_j = np.max(np.abs(feedback[2:]), axis=(0, 2))
    error_A = float(np.max(np.abs(native_A-grid_A)))
    error_j = np.max(np.abs(feedback[2:]-current[1:-1]), axis=(0, 2))
    # Zero signals must agree exactly; never hide a mismatch with a large floor.
    relative_A = error_A/scale_A if scale_A else (0. if error_A == 0 else float('inf'))
    relative_j = np.divide(error_j, scale_j, out=np.full_like(error_j, np.inf), where=scale_j != 0)
    relative_j[(scale_j == 0) & (error_j == 0)] = 0.
    initial_zero = bool(np.all(feedback[:2] == 0))
    return dict(relative_grid_A_error=relative_A,
                relative_current_error_per_point=relative_j.tolist(),
                max_abs_grid_A=scale_A, max_abs_current_per_point=scale_j.tolist(),
                initial_feedback_zero=initial_zero, tolerance=tolerance,
                passed=bool(initial_zero and relative_A <= tolerance and np.all(relative_j <= tolerance)))


if __name__ == '__main__':
    main()
