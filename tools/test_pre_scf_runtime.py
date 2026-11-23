"""Compare tiny GS warmup transitions with an audited original fast-HSE build.

Control-flow regression only: loose thresholds are NOT physical convergence.
All generated restart payloads are removed after comparing scalar results.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess


def audited(path):
    record = json.loads(path.read_text())
    executable = Path(record['executable'])
    if (record.get('schema') != 'salmonpy.hse-rt-state-build/1'
            or record.get('status') != 'built-and-audited'
            or hashlib.sha256(executable.read_bytes()).hexdigest() != record['sha256']):
        raise ValueError('unchanged audited build required')
    return record, executable


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build', type=Path)
    parser.add_argument('reference', type=Path)
    parser.add_argument('output', type=Path)
    parser.add_argument('--dc-ranks', type=int, choices=(2, 4),
                        help='two DC fragments, with one or two spatial ranks each')
    parser.add_argument('--mpiexec', default='mpiexec')
    parser.add_argument('--require-gs-owner', action='store_true',
                        help='require explicit GS HSE owner in the current build only')
    args = parser.parse_args()
    records = [audited(p) for p in (args.build, args.reference)]
    root = Path(__file__).resolve().parents[1]
    output = args.output.resolve()
    protected = [root]
    for record, executable in records:
        protected.extend((executable.parent, Path(record['source_root'])))
    if any(output.is_relative_to(p.resolve()) for p in protected):
        raise ValueError('new external output directory required')
    pseudo = root / 'salmonpy/examples/si_intro_data/Si_rps.dat'
    if not pseudo.exists():
        pseudo = root.parent / 'examples/si_intro_data/Si_rps.dat'
    pseudo_bytes = pseudo.read_bytes()
    output.mkdir(parents=True, exist_ok=False)
    env = dict(os.environ, OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    results = []
    for steps in (0, 1, 2, 3):
        runs = []
        for label, (_, executable) in zip(('current', 'reference'), records):
            directory = output / f'steps-{steps}' / label
            directory.mkdir(parents=True)
            (directory / 'Si_rps.dat').write_bytes(pseudo_bytes)
            text = f"""&calculation
 theory='dft'
/
&control
 sysname='gs', base_directory='./', yn_restart='n', yn_self_checkpoint='n'
/
&units
 unit_system='au'
/
&system
 yn_periodic='y', spin='unpolarized', nelem=1, natom=1, nelec=4, nstate=2, al=8,8,8
/
&pseudo
 file_pseudo(1)='Si_rps.dat', izatom(1)=14, lloc_ps(1)=2
/
&functional
 xc='hse06', exx_pre_scf_threshold={1e8 if steps else 0}, exx_pre_scf_steps={max(steps, 1)}, yn_hse_wannier='n'
/
&rgrid
 num_rgrid=8,8,8
/
&kgrid
 num_kgrid=2,1,1
/
&scf
 nscf=10, convergence='rho_dne', threshold=1d8, nscf_init_mix_zero=1, iseed_number_change=0
/
&atomic_red_coor
'Si' .3 .4 .45 1
/
"""
            if args.dc_ranks:
                text = text.replace("yn_self_checkpoint='n'", "yn_self_checkpoint='n', write_gs_restart_data='no'")
                text = text.replace("theory='dft'", "theory='dft', yn_dc='y'")
                text = text.replace('natom=1, nelec=4, nstate=2, al=8,8,8',
                                    'natom=2, nelec=8, nstate=4, al=16,8,8, temperature_k=10000')
                text = text.replace('num_rgrid=8,8,8', 'num_rgrid=16,8,8')
                text = text.replace('num_kgrid=2,1,1', 'num_kgrid=1,1,1')
                text = text.replace('&atomic_red_coor', '&atomic_coor')
                text = text.replace("'Si' .3 .4 .45 1", "'Si' 4 3.2 3.6 1\n'Si' 12 3.2 3.6 1")
                spatial = args.dc_ranks // 2
                text += f'''&dc
 num_fragment=2,1,1
 num_rgrid_buffer=4,0,0
 nproc_rgrid_tot=2,{spatial},1
 nstate_frag=4
 yn_dc_lcfo='n'
/
&parallel
 nproc_k=1, nproc_ob=1, nproc_rgrid=1,{spatial},1
/
'''
            (directory / 'gs.inp').write_text(text)
            try:
                with (directory / 'gs.inp').open() as inp, (directory / 'run.log').open('w') as log:
                    command = [str(executable)]
                    if args.dc_ranks:
                        command = [args.mpiexec, '-n', str(args.dc_ranks), *command]
                    subprocess.run(command, cwd=directory, env=env,
                                   stdin=inp, stdout=log, stderr=subprocess.STDOUT,
                                   timeout=120, check=True)
                log = (directory / 'run.log').read_text()
                if label == 'current' and args.require_gs_owner and 'EXX_GS_OWNER explicit' not in log:
                    raise AssertionError('current GS did not select explicit HSE ownership')
                marker = 'end SALMON' if args.dc_ranks else '#GS converged'
                if marker not in log or re.search(r'\b(?:NaN|Infinity)\b', log, re.I):
                    raise AssertionError(f'invalid GS result: {directory}')
                switch = re.findall(r'EXX_PRE_SCF switch to target hybrid before iteration (\d+)', log)
                if switch != ([str(steps + 1)] if steps else []):
                    raise AssertionError(f'unexpected phase transition: {switch}')
                diagnostics = [line.strip() for line in log.splitlines()
                               if line.startswith('EXX_PRE_SCF') or 'Total Energy=' in line
                               or '#GS converged' in line or line.startswith('DC #SCF')]
                folders = [directory]
                if args.dc_ranks:
                    folders = [directory / 'data_dcdft/fragments' / f'{i:06d}' for i in (1, 2)]
                    iterations = re.findall(r'DC #SCF =\s*(\d+)', log)
                    if not iterations or int(iterations[-1]) >= 10:
                        raise AssertionError('DC did not exit SCF early')
                scalar = [(folder / name).read_bytes() for folder in folders
                          for name in ('gs_info.data', 'gs_eigen.data', 'gs_k.data')]
                if any(re.search(rb'\b(?:NaN|Infinity)\b', data, re.I) for data in scalar):
                    raise AssertionError('nonfinite scalar output')
                runs.append((scalar, diagnostics))
            finally:
                for payload in directory.rglob('data_for_restart'):
                    shutil.rmtree(payload)

                (directory / 'Si_rps.dat').unlink(missing_ok=True)
        if runs[0] != runs[1]:
            raise AssertionError(f'GS scalar output or phase trace differs for steps={steps}')
        results.append({'pre_scf_steps': steps, 'scalar_outputs_byte_equal': True,
                        'phase_and_energy_trace_equal': True})
    report = {'schema': 'salmonpy.pre-scf-regression/1', 'status': 'passed',
              'scope': 'tiny GS/DC control-flow; not physical convergence or concurrent-state isolation',
              'mpi_ranks': args.dc_ranks or 0, 'require_gs_owner': args.require_gs_owner,
              'builds': [record for record, _ in records], 'cases': results}
    (output / 'result.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'GS pre-SCF comparison PASS (4 cases, MPI ranks={args.dc_ranks or 0})')


if __name__ == '__main__':
    main()
