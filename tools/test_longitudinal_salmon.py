"""Compare salmonpy.maxwell field components + TDDFT worker with SALMON's own fields.

PYTHONPATH=. python tools/test_longitudinal_salmon.py SALMON_EXECUTABLE NEW_DIRECTORY [--trans-longi lo|2d|ms1d]

The same synthetic, normalized (not converged) Si fixture is propagated twice:
  1. ordinary SALMON with trans_longi='lo' or '2d', or the 1D multiscale solver with
     one macroscopic point ('ms1d'): the field is computed inside SALMON;
  2. the SALMONpy TDDFT worker (trans_longi='tr') coupled to salmonpy.maxwell.
Ac_tot and the matter current must agree to SALMON's printed precision.
Requires SALMONPY_RUNTIME_SOURCE (prepared runtime source, for Si_rps.dat).
"""
import argparse
import hashlib
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
    args = p.parse_args()
    runtime_source = os.environ.get('SALMONPY_RUNTIME_SOURCE')
    if not runtime_source: raise SystemExit('set SALMONPY_RUNTIME_SOURCE')
    dest = Path(args.directory); dest.mkdir(exist_ok=False)
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
    maxwell = ("&multiscale\n fdtddim='1d', nx_m=1, ny_m=1, nz_m=1, "
               f"hx_m={LAYER['spacing']}, hy_m={LAYER['spacing']}, hz_m={LAYER['spacing']}, "
               # SALMON's checker requires nxvacl_m/nxvacr_m: nxvac = (|nxvacl_m|+1, |nxvacr_m|)
               f"nxvacl_m={LAYER['vacuum'][0]-1}, nxvacr_m={LAYER['vacuum'][1]}\n/\n"
               "&maxwell\n boundary_em(1,1)='abc', boundary_em(1,2)='abc'\n/\n"
               "&analysis\n out_ms_step=1\n/\n")
    (native/'salmon.inp').write_text(inp('tr', extra, maxwell, theory='multi_scale_maxwell_tddft'))
    with (native/'stdout.log').open('w') as log:
        subprocess.run([str(exe)], stdin=(native/'salmon.inp').open(), stdout=log, stderr=subprocess.STDOUT,
                       cwd=native, check=True, env={**os.environ, 'OMP_NUM_THREADS': '1'})
    py = dest/'salmonpy'; fixture(py, runtime_source)
    (py/'rt.inp').write_text(inp('tr', ''))
    pulse = AcosPulse(**PULSE)
    field = Multiscale1D(cells=1, incident=pulse, **LAYER).new_state(dt=DT, steps=NT)
    grids = [field.grid_A.copy()]; charge = []
    with TDDFTInitialState.from_salmon(py, 'rt.inp', executable=SalmonExecutable(exe, sha),
                                       max_steps=NT, threads=1, timeout=600.) as initial:
        with initial.new_state(dt=DT, initial_A=field.A[0]) as state:
            j = np.asarray(state.current)  # TDDFTState returns charge current j
            for n in range(NT):
                charge.append(j)
                A = field.step(j[None])
                grids.append(field.grid_A.copy())
                j = np.asarray(state.step(A[0]))
    errA = errJ = 0.; scaleA = scaleJ = 0.
    for itt in range(NT+1):
        table = np.loadtxt(native/'fixture_RT_Ac'/f'fixture_Ac_{itt:06d}.data', comments='#')
        ix = table[:, 0].astype(int)
        if not np.array_equal(ix, np.arange(1-LAYER['vacuum'][0], 2+LAYER['vacuum'][1])):
            raise SystemExit('unexpected multiscale grid')
        scaleA = max(scaleA, np.abs(table[:, 3:6]).max())
        errA = max(errA, np.abs(table[:, 3:6]-grids[itt]).max())
        if itt >= 2:  # file itt holds j_em = -curr(itt-1); j(0) is never fed back
            cell = np.where(ix == 1)[0][0]
            scaleJ = max(scaleJ, np.abs(table[cell, 12:15]).max())
            errJ = max(errJ, np.abs(table[cell, 12:15]-charge[itt-1]).max())
    print(f'mode=ms1d steps={NT} grid cells={len(ix)} max|Ac|={scaleA:.3e} max|j|={scaleJ:.3e}')
    print(f'relative max difference: Ac on whole grid {errA/scaleA:.2e}, current {errJ/scaleJ:.2e}')
    if errA/scaleA > 1e-12 or errJ/scaleJ > 1e-12: raise SystemExit('FAIL: SALMON and SALMONpy differ')
    print('PASS')


if __name__ == '__main__':
    main()
