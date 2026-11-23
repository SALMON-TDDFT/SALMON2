"""Short Ni SOC test using an existing audited Libxc-PBE build; no GS export.

BUILD_JSON NEW_OUTPUT --potential PINNED_NI_UPF
Use --allow-download instead only with explicit authorization on a connected node.
The test performs three SCF iterations; it does not certify physical convergence.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import time
import warnings

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,_read_status
from salmonpy.preparation import prepare_dft
from salmonpy.pseudopotentials import catalog_reference,download


def nickel_model(potential):
    # One-atom fcc primitive cell. Lattice constant is a test geometry.
    a=3.52/0.52917721067
    cell=[[0,a/2,a/2],[a/2,0,a/2],[a/2,a/2,0]]
    return DFT([28],cell,[[0,0,0]],[potential],18,36,(16,16,16),(2,2,2),
               xc='libxc_pbe',spin_orbit=True,electronic_temperature_k=300)


def run(build_json,output,*,potential=None,allow_download=False):
    build=json.loads(Path(build_json).read_text())
    executable=Path(build['executable']).resolve()
    if build.get('status')!='built-and-audited' or hashlib.sha256(executable.read_bytes()).hexdigest()!=build.get('sha256'):
        raise ValueError('unchanged audited executable required')
    cache=(executable.parent/'CMakeCache.txt').read_text()
    import re
    if not re.search(r'^USE_LIBXC:\w+=(ON|TRUE|1)$',cache,re.M):
        raise ValueError('Ni PBE requires a Libxc-enabled build; no XC substitution')
    if potential is None and allow_download is not True:
        raise ValueError('provide --potential or explicitly authorize --allow-download')
    output=Path(output).expanduser().resolve();output.mkdir(parents=True,exist_ok=False)
    reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-ni')
    if potential is None:
        potential=download(reference,output/'potential',allow_download=True)
    pseudo=Pseudopotential.from_file(potential,atomic_number=28,local_channel=2,reference=reference)
    model=nickel_model(pseudo)
    scf=SCFSettings(max_iterations=3,tolerance=1e-12,mixing_method='simple_dm',mixing_rate=.2,
                    initial_wavefunction='random',random_seed=0)
    prepared=prepare_dft(model,output/'input',scf=scf,save_restart=False)
    command=[sys.executable,str(ROOT/'tools/run_nlcc_native_tests.py'),str(executable.parent),
             str(output/'execution'),'--test','test_nickel_soc','--input',str(prepared/'gs.inp'),
             '--asset','pseudo1.upf',str(prepared/'pseudo1.upf')]
    started=time.perf_counter()
    subprocess.run(command,check=True)
    evidence=json.loads((output/'execution/test-record.json').read_text())
    cases=[]
    for case in evidence['runs']:
        if case['returncode'] or any(marker not in case['stdout'] for marker in
            ('NI_SOC_SHORT_SCF_PASS','NI_SOC_RESTORE_EXACT_PASS')):
            raise AssertionError('Ni SOC native test failed')
        work=output/'execution'/f"ranks{case['ranks']}-threads{case['threads']}"
        with warnings.catch_warnings(record=True) as messages:
            warnings.simplefilter('always')
            status=_read_status(work/'gs_status.dat')
        if not math.isfinite(status['energy_hartree']) or status['iterations']!=3:
            raise AssertionError('invalid short SCF result')
        exported=list((work/'data_for_restart').rglob('*')) if (work/'data_for_restart').exists() else []
        if any(p.is_file() for p in exported):
            raise AssertionError('short CI fixture exported GS data')
        cases.append({'ranks':case['ranks'],'threads':case['threads'],**status,
                      'nonconvergence':'expected for a three-iteration smoke test; no RT continuation',
                      'warnings':[str(m.message) for m in messages],
                      'soc_state_restore':'bitwise exact','gs_exported':False})
    report={'schema':'salmonpy.nickel-soc-smoke/1','status':'passed','kmesh':[2,2,2],
            'k_generation':'SALMON standard','reference':reference.catalog_id,
            'semicore_shells':list(reference.semicore_shells),'valence_configuration':reference.valence_configuration,
            'potential_sha256':pseudo.sha256,'runtime_sha256':build['sha256'],
            'test_executable_sha256':evidence['executable_sha256'],'cases':cases,
            'elapsed_seconds':time.perf_counter()-started,
            'initialization':'seeded random spinors, seed0; no prescribed magnetization vector',
            'electronic_temperature_k':300,'scope':'initialization and noncollinear-array ownership only; not full SCF/TDDFT checkpoint or converged Ni magnetic properties'}
    (output/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print('Ni SOC short test PASS:',output/'result.json')
    return report


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('build',type=Path);parser.add_argument('output',type=Path)
    source=parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--potential',type=Path)
    source.add_argument('--allow-download',action='store_true')
    args=parser.parse_args()
    run(args.build,args.output,potential=args.potential,allow_download=args.allow_download)


if __name__=='__main__':
    main()
