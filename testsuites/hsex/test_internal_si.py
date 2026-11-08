#!/usr/bin/env python3
"""Small real-Si SALMON integration validation; Python is only a test harness.

Keep WORKDIR on scratch. No GS/checkpoint payload is collected by this harness.
The separate shell runner makes MPI launch explicit on both Mac and Fugaku.
"""
from pathlib import Path
import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[2]
FIXTURES=Path(__file__).with_name('internal_si')
ATOL=2e-12
RTOL=2e-11
CURRENT_ATOL=1e-12


def prepare(args):
    work=args.work.resolve()
    if work.exists():
        raise FileExistsError('Use a new scratch directory; previous evidence is retained')
    pseudo=args.pseudo.resolve()
    if not pseudo.is_file():
        raise FileNotFoundError(pseudo)
    work.mkdir(parents=True)
    gs=(FIXTURES/'gs.inp').read_text()
    rt=(FIXTURES/'rt.inp').read_text().replace('hsex_w_eta=0.01',f'hsex_w_eta={args.eta}')
    restart=str(work/'gs/data_for_restart')+'/'
    variants={
        'updated':("yn_hsex_update_w='y', hsex_screen_step=1",restart),
        'snapshot':("yn_hsex_update_w='n', hsex_screen_step=1",restart),
        'restart':("yn_hsex_update_w='y', hsex_screen_step=1",str(work/'updated/checkpoint_rt_000001')+'/'),
        'external':("yn_hsex_update_w='n', hsex_screen_step=1, hsex_probe_w='python-pump.bin'",restart),
    }
    for stage in ('gs',*variants):
        directory=work/stage;directory.mkdir()
        shutil.copyfile(pseudo,directory/'Si_rps.dat')
        if stage=='gs':
            text=gs
        else:
            update,data=variants[stage]
            text=rt.replace('@RESTART@',data).replace('@UPDATE@',update)
            if stage=='restart':text=text.replace('&control','&control\n yn_restart=\'y\'',1)
        (directory/'input').write_text(text)
    (work/'reference').mkdir()
    shutil.copyfile(FIXTURES/'run-stage.sh',work/'run-stage.sh')
    (work/'run-stage.sh').chmod(0o755)
    (work/'fixture.json').write_text(json.dumps(dict(grid=[8]*3,kmesh=[2]*3,qmesh=[2]*3,
        boxes=2,cutoff_Ha=1.,eta_input_Ha=args.eta,mpi=4,omp=12,
        pseudo_sha256=hashlib.sha256(pseudo.read_bytes()).hexdigest()),indent=2)+'\n')
    print(work)


def run(args):
    subprocess.run(['bash',str(args.work/'run-stage.sh'),str(args.work.resolve()),
                    str(args.executable.resolve()),args.sha256,args.stage],check=True)


def reference(args):
    sys.path.insert(0,str(ROOT/'tools/hsex'))
    from prepare_w import prepare as build_w
    directory=args.work/'reference';directory.mkdir(exist_ok=True)
    eta=json.loads((args.work/'fixture.json').read_text())['eta_input_Ha']
    if args.kind=='gs':
        build_w(args.work/'gs',directory/'gs.bin',[2]*3,2,1.,qmesh=[2]*3,hamiltonians=args.work/'gs')
    else:
        build_w(args.work/'gs',directory/'pump.bin',[2]*3,2,1.,qmesh=[2]*3,
                pump=args.work/args.snapshot,eta=eta)


def read_w(path):
    data=path.read_bytes()
    if data[:4] not in (b'HXR1',b'HXT1',b'HXP1'):
        raise ValueError(f'Unexpected W format: {path}')
    dims=np.frombuffer(data,'=i4',count=9,offset=4)
    ng=int(np.prod(dims[:3]));nk=int(np.prod(dims[3:6]));na=int(dims[6])**3;nion=int(dims[8])
    offset=40+8*(9+3*nk+3*nion+ng+1+3*nk)
    end=offset+8*nk*na*na
    return data[:offset],np.frombuffer(data,'=f8',count=nk*na*na,offset=offset),data[end:]


def compare_w(actual,reference):
    ah,a,at=read_w(actual);bh,b,bt=read_w(reference)
    if ah!=bh:raise AssertionError(f'W grid/reference metadata differs: {actual}')
    if not np.isfinite(a).all() or not np.isfinite(b).all():raise AssertionError('Nonfinite W')
    # Tail contains pump step, projection size, eigenvalues and occupations.
    # Independent OMP projection reductions may differ by roundoff.
    if at or bt:
        if at[:8]!=bt[:8] or len(at)!=len(bt):raise AssertionError('Pump projection metadata mismatch')
        ta=np.frombuffer(at,'=f8',offset=8);tb=np.frombuffer(bt,'=f8',offset=8)
        if not np.isfinite(ta).all() or not np.isfinite(tb).all():raise AssertionError('Nonfinite pump metadata')
        np.testing.assert_allclose(ta,tb,rtol=0,atol=1e-7,equal_nan=False)
    np.testing.assert_allclose(a,b,atol=ATOL,rtol=RTOL,equal_nan=False)
    return dict(max_absolute_error=float(np.max(abs(a-b))),
                maximum_tolerance_ratio=float(np.max(abs(a-b)/(ATOL+RTOL*abs(b)))))


def current_rows(path,steps):
    values=np.atleast_2d(np.loadtxt(path))
    if values.shape!=(len(steps),16) or not np.isfinite(values).all():
        raise AssertionError(f'Incomplete/nonfinite current trajectory: {path}')
    np.testing.assert_allclose(values[:,0],np.asarray(steps)*.01,rtol=0,atol=1e-10)
    return values


def completed_log(path,require_team=False):
    text=path.read_text()
    if not re.search(r'^end SALMON\s*$',text,re.M):
        raise AssertionError(f'SALMON completion marker missing: {path}')
    if require_team and not re.search(r'HSEX_W_OMP team=\s*12(?:\s|$)',text):
        raise AssertionError(f'Actual W OpenMP team 12 missing: {path}')


def compare_current(actual,reference,steps):
    a=current_rows(actual,steps);b=current_rows(reference,[1,2,3])
    error=0.
    for row in a:
        selected=b[np.isclose(b[:,0],row[0],rtol=0,atol=1e-10)]
        if len(selected)!=1:raise AssertionError('Missing or duplicate current time')
        np.testing.assert_allclose(row[13:16],selected[0,13:16],atol=CURRENT_ATOL,rtol=0)
        error=max(error,float(np.max(abs(row[13:16]-selected[0,13:16]))))
    return dict(rows=len(a),max_absolute_error=error,absolute_tolerance=CURRENT_ATOL)


def check(args):
    p=args.work
    manifest=json.loads((p/'fixture.json').read_text())
    expected=dict(grid=[8]*3,kmesh=[2]*3,qmesh=[2]*3,boxes=2,cutoff_Ha=1.,eta_input_Ha=args.eta,mpi=4,omp=12)
    if any(manifest.get(key)!=value for key,value in expected.items()):
        raise AssertionError('Fixture metadata differs from the required 4 MPI x 12 OMP integration case')
    for directory,log,team in [('gs','run.log',False),(args.updated,args.updated_log,True),
                               (args.snapshot,'run.log',False),(args.restart,'run.log',False)]:
        completed_log(p/directory/log,require_team=team)
    if not re.search(r'HSEX_INTERNAL_W_UPDATE_PASS step=\s*1\b',(p/args.updated/args.updated_log).read_text()):
        raise AssertionError('Internal W update completion marker missing')
    if not re.search(r'HSEX_ACTIVE_W_RESTORED step=\s*1\b',(p/args.restart/'run.log').read_text()):
        raise AssertionError('Updated W restore marker missing')
    updated=current_rows(p/args.updated/'Si_rt.data',[1,2,3])
    frozen=current_rows(p/args.snapshot/'Si_rt.data',[1,2,3])
    effect=float(np.max(abs(updated[:,13:16]-frozen[:,13:16])))
    if effect<=5*CURRENT_ATOL:raise AssertionError('W update has no resolved effect on current')
    gsref=args.gs_reference or p/'reference/gs.bin'
    pumpref=args.pump_reference or p/'reference/pump.bin'
    report=dict(gs_w=compare_w(p/'gs/hsex-w.bin',gsref),
        pump_w=compare_w(p/args.updated/'hsex-updated-000000001.bin',pumpref),
        restart_current=compare_current(p/args.restart/'Si_rt.data',p/args.updated/'Si_rt.data',[2,3]))
    if args.skip_external:
        report['external_w_current']='NOT RUN'
    else:
        completed_log(p/args.external/'run.log')
        report['external_w_current']=compare_current(p/args.external/'Si_rt.data',p/args.updated/'Si_rt.data',[1,2,3])
    report['update_current_effect']=dict(max_absolute_difference=effect,required_minimum=5*CURRENT_ATOL)
    report['w_tolerance']=dict(atol=ATOL,rtol=RTOL)
    report['eta_input_Ha']=args.eta
    report['scope']='Small integration fixture; not production convergence or physical Si spectrum validation'
    (p/'integration-check.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));print('INTERNAL_SI_PASS')


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    commands=parser.add_subparsers(dest='command',required=True)
    p=commands.add_parser('prepare');p.add_argument('work',type=Path);p.add_argument('--pseudo',type=Path,required=True)
    p.add_argument('--eta',type=float,choices=[.001,.01],default=.001)
    p.set_defaults(action=prepare)
    p=commands.add_parser('run');p.add_argument('work',type=Path);p.add_argument('stage',choices=['gs','updated','snapshot','restart','external'])
    p.add_argument('--executable',type=Path,required=True);p.add_argument('--sha256',required=True);p.set_defaults(action=run)
    p=commands.add_parser('reference');p.add_argument('work',type=Path);p.add_argument('kind',choices=['gs','pump'])
    p.add_argument('--snapshot',default='snapshot');p.set_defaults(action=reference)
    p=commands.add_parser('check');p.add_argument('work',type=Path)
    p.add_argument('--snapshot',default='snapshot');p.add_argument('--updated-log',default='run.log')
    p.add_argument('--updated',default='updated');p.add_argument('--restart',default='restart');p.add_argument('--external',default='external')
    p.add_argument('--gs-reference',type=Path);p.add_argument('--pump-reference',type=Path);p.add_argument('--skip-external',action='store_true')
    p.add_argument('--eta',type=float,choices=[.001,.01],default=.001)
    p.set_defaults(action=check)
    args=parser.parse_args();args.action(args)

if __name__=='__main__':main()
