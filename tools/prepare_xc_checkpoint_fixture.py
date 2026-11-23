"""Explicit small Si XC checkpoint regression preparation, not optical convergence.

Usage: python -I tools/prepare_xc_checkpoint_fixture.py AUDITED_EXECUTABLE NEW_ROOT [INITIAL_GUESS_RESTART]
Runs each XC's own SCF with unchanged small starter discretization. TB-mBJ is
prepared directly through the solver input because DFT's public wrapper is narrower.
"""
import hashlib,json,os,shutil,subprocess,sys
import numpy as np
from pathlib import Path
from salmonpy.examples.si_intro import _model
from salmonpy.bands import _input_groups

exe=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]).resolve();out.mkdir(exist_ok=False)
assets=Path(__file__).resolve().parents[1]/'salmonpy/examples/si_intro_data'
config=json.loads((assets/'input.json').read_text())
for tag,xc,cval in [('tbmbj-auto','tbmbj',None),('tbmbj-fixed','tbmbj',1.2),
                    ('libxc-pz','libxc_pz',None),('libxc-pbe','libxc_pbe',None)]:
    d=out/tag;d.mkdir();shutil.copy2(assets/'Si_rps.dat',d/'Si_rps.dat')
    model,scf=_model(config,d)
    shutil.copy2(d/'Si_rps.dat',d/'pseudo1_rps.dat')
    text=model.input_text(scf).replace("xc='pz'",f"xc='{xc}'"+(f', cval={cval}' if cval is not None else ''))
    if len(sys.argv) > 3:
        guess=Path(sys.argv[3]).resolve()
        (d/'initial_guess').mkdir()
        for name in ('info.bin','occupation.bin','wfn.bin'):
            shutil.copy2(guess/name,d/'initial_guess'/name)
        text=text.replace("yn_restart='n'", "yn_restart='y', read_gs_restart_data='wfn', yn_reset_step_restart='y', directory_read_data='./initial_guess/'",1)
    (d/'gs.inp').write_text(text)
    result=subprocess.run([str(exe)],input=text,text=True,cwd=d,
        env={**os.environ,'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'},
        stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=600)
    (d/'gs.log').write_text(result.stdout)
    if result.returncode or 'end SALMON' not in result.stdout or '#GS converged' not in result.stdout or not (d/'data_for_restart/wfn.bin').exists():
        raise RuntimeError(f'{tag}: SCF preparation failed; inspect {d}/gs.log')
    if not np.isfinite(np.fromfile(d/'data_for_restart/wfn.bin',dtype=np.complex128)).all() or 'NaN' in result.stdout:
        raise RuntimeError(f'{tag}: nonfinite SCF result; inspect {d}/gs.log')
    groups=_input_groups(text);groups['calculation']="\n theory='tddft_pulse'\n"
    groups['control']="\n sysname='rt', directory_read_data='./data_for_restart/', yn_restart='n'\n"
    groups.pop('scf',None);groups['tgrid']='\n dt=.003, nt=128\n'
    groups['emfield']="\n ae_shape1='none', trans_longi='tr'\n"
    groups['propagation']="\n propagator='middlepoint', yn_predictor_corrector='n'\n"
    (d/'rt.inp').write_text('\n'.join('&'+key+body+'/' for key,body in groups.items())+'\n')
    (d/'fixture.json').write_text(json.dumps({'xc':xc,'cval':cval,
        'executable_sha256':hashlib.sha256(exe.read_bytes()).hexdigest(),
        'initial_guess': {name:hashlib.sha256((d/'initial_guess'/name).read_bytes()).hexdigest()
                          for name in ('info.bin','occupation.bin','wfn.bin')} if len(sys.argv)>3 else None,
        'scope':'small Si numerical checkpoint regression, no optical convergence claim'},indent=2)+'\n')
    print(tag,flush=True)
