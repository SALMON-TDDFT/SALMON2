#!/usr/bin/env python3
"""Sequential fragment x orbital MPI native integration test. Needs a USE_HSE binary and H_rps.dat.
Usage: python3 test_native.py --binary /abs/salmon --pseudo /abs/H_rps.dat
Each run owns a fresh temporary directory; no production data is modified.
"""
import argparse
import json
import math
import re
import os
from pathlib import Path
import shutil
import subprocess
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--binary", type=Path, required=True)
parser.add_argument("--pseudo", type=Path, required=True)
parser.add_argument("--mpirun", default="mpirun")
parser.add_argument("--work-root", type=Path, default=Path(tempfile.gettempdir()))
parser.add_argument("--axis", choices=["x","y","z"], default="x")
parser.add_argument("--xc", choices=["pbeh40","pbeh40_rvv10","hse06"], default="pbeh40_rvv10")
parser.add_argument("--omp", type=int, default=1)
args = parser.parse_args()
assert args.omp > 0
root = Path(tempfile.mkdtemp(prefix="lcfo-native-test-", dir=args.work_root))
repo = Path(__file__).resolve().parents[2]
env = dict(os.environ, OMP_NUM_THREADS=str(args.omp), OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
command = [args.mpirun, "-np", "2", str(args.binary.resolve())]

def run(name, text, rt=False, reject=False, settings=None, orbital_groups=1):
    folder = root / name
    folder.mkdir()
    shutil.copy(args.pseudo, folder / "H_rps.dat")
    if rt:
        (folder / "data_dcdft").symlink_to(root / "gs/data_dcdft", target_is_directory=True)
    local_env = dict(env)
    controls = dict(settings or {})
    if rt:
        controls.update(yn_hse_lcfo_rt='y', yn_hse_lcfo_continuity='y')
    for key, value in controls.items():
        literal = "'" + value + "'" if key.startswith('yn_') else str(value)
        pattern = r"(?i)(\b" + re.escape(key) + r"\s*=\s*)(?:'[^']*'|[^,\s/]+)"
        if re.search(pattern, text):
            text = re.sub(pattern, lambda match: match[1] + literal, text)
        else:
            text = text.replace('&functional', '&functional\n ' + key + '=' + literal, 1)
    permutation={'x':(0,1,2),'y':(1,0,2),'z':(2,1,0)}[args.axis]
    def rotate(match):
        values=match[2].split(',')
        return match[1]+','.join(values[j] for j in permutation)
    for key in ('num_fragment','num_rgrid_buffer','nproc_rgrid_tot','nproc_rgrid','al','num_rgrid','epdir_re1'):
        text=re.sub(r'(?im)(^\s*'+key+r'\s*=\s*)([^\n/]+)',rotate,text)
    def rotate_atom(match):
        values=match[2].split()
        return match[1]+' '.join(values[j] for j in permutation)+' '+values[3]
    text=re.sub(r"(?m)(^\s*'H'\s+)([^\n]+)",rotate_atom,text)
    (folder / "inputfile").write_text(text)
    local_command=command.copy();local_command[2]=str(2*orbital_groups)
    with (folder / "inputfile").open("rb") as inp, (folder / "run.log").open("wb") as log:
        status = subprocess.run(local_command, cwd=folder, env=local_env, stdin=inp,
                                stdout=log, stderr=subprocess.STDOUT, timeout=120)
    log = (folder / "run.log").read_text()
    if reject:
        assert (reject if isinstance(reject,str) else "LCFO RT: restart/checkpoint output is not supported yet") in log, folder
        if reject is True:
            assert "Native LCFO RT active" not in log, folder
    else:
        assert status.returncode == 0 and "end SALMON" in log, folder
        if rt:
            assert "Native LCFO RT active" in log and "LCFO HSE ACE build" in log, folder
            storage=re.findall(r'LCFO distributed storage rank/local/global/halo rows:\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)',log)
            assert len(storage)==2, (folder,storage)
            assert all(0<int(local)<int(total) and 0<int(halo)<=int(total) for rank,local,total,halo in storage)
            assert 'LCFO timing pack/source/exchange/ACE build' in log, folder
    return folder

def rows(file):
    data = [[float(x) for x in line.split()] for line in file.read_text().splitlines()
            if line.strip() and not line.startswith("#")]
    assert all(math.isfinite(x) for row in data for x in row), file
    return data

source = (repo / "testsuites/422_H_dcdft_hse/inputfile").read_text()
source = source.replace("nproc_rgrid_tot=4,1,1", "nproc_rgrid_tot=2,1,1")
source = source.replace("nproc_k=2", "nproc_k=1").replace("num_kgrid=1,2,1", "num_kgrid=1,1,1")
source = source.replace("xc='hse06'", "xc='"+args.xc+"'")
run("gs", source)
rt = source.replace("theory='dft'", "theory='tddft_response'")
rt = rt.replace("yn_dc='y'", "yn_dc='n'\n yn_conventional_from_dcdft='y'")
rt = rt.replace("nproc_rgrid=1,1,1", "nproc_rgrid=2,1,1")
rt = rt.replace(" nstate=4", " nstate=2").replace(" temperature_k=300d0\n", "")
rt += """
&tgrid
 nt=4
 dt=0.02d0
/
&emfield
 ae_shape1='impulse'
 e_impulse=0.0001d0
 epdir_re1=1d0,0d0,0d0
/
&analysis
 yn_out_dns_rt='y'
 out_dns_rt_step=4
 out_rt_energy_step=1
 nenergy=20
 de=0.01d0
/
"""
coarse = run("rt", rt, rt=True)
fine = run("rt_half_dt", rt.replace("nt=4", "nt=8").replace("dt=0.02d0", "dt=0.01d0"), rt=True)
a = rows(coarse / "H_dc_hse_rt.data")
b = rows(fine / "H_dc_hse_rt.data")
assert len(a) == 4 and len(b) == 8
assert abs(a[-1][0] - b[-1][0]) < 1e-12
current_error = max(abs(a[-1][j] - b[-1][j]) for j in range(13, 16))
assert current_error < 1e-10, current_error
energy = rows(coarse / "H_dc_hse_rt_energy.data")
energy_width = max(row[1] for row in energy[1:]) - min(row[1] for row in energy[1:])
assert energy_width < 1e-9, energy_width
run("reject_restart", rt.replace("&control", "&control\n yn_restart='y'"), rt=True, reject="PBEh40: checkpoint parameter validation" if args.xc != "hse06" else True)
print(json.dumps(dict(work=str(root), steps=4, half_dt_current_error=current_error,
                      post_impulse_energy_width=energy_width, restart_guard="passed"), indent=2))


split = run('rt_orbital2',rt.replace('nproc_ob=1','nproc_ob=2'),rt=True,orbital_groups=2)
split_rows=rows(split/'H_dc_hse_rt.data')
assert len(split_rows)==len(a)
assert max(abs(x[j]-y[j]) for x,y in zip(a,split_rows) for j in range(13,16)) < 1e-11
split_energy=rows(split/'H_dc_hse_rt_energy.data')
assert len(split_energy)==len(energy)
assert max(abs(x[1]-y[1]) for x,y in zip(energy,split_energy)) < 1e-9
if args.xc!='hse06':
    for key,val in [('pbeh_coulomb_radius','3'),('rvv10_b','6.0')]:
        if key=='rvv10_b' and args.xc!='pbeh40_rvv10':continue
        run('reject_'+key,rt,rt=True,settings={key:val},reject='LCFO functional metadata missing or mismatched')
    files=sorted((root/'gs/data_dcdft/fragments').glob('*/functional.txt'))
    assert len(files)==2
    saved=files[-1].read_text()
    for mode in ('missing','run_id','finite_source'):
        if mode=='missing':files[-1].unlink()
        else:
            lines=saved.splitlines()
            if mode=='run_id':lines[1]='stale-run'
            else:
                vals=lines[3].split();vals[-1]='1.0';lines[3]=' '.join(vals)
            files[-1].write_text('\n'.join(lines)+'\n')
        try:
            run('reject_metadata_'+mode,rt,rt=True,reject='LCFO functional metadata missing or mismatched')
        finally:files[-1].write_text(saved)
    other='pbeh40' if args.xc=='pbeh40_rvv10' else 'pbeh40_rvv10'
    run('reject_functional',rt.replace(args.xc,other),rt=True,reject='LCFO functional metadata missing or mismatched')
print('Orbital decomposition and functional provenance passed')

# Grid spacing input leaves the global num_rgrid namelist at zero. Derive the
# distributed extent from owned grid indices, as for spatial stencil layouts.
spacing=run('rt_grid_spacing',rt.replace('num_rgrid=16,8,8','dl=1d0,1d0,1d0'),rt=True)
spacing_rows=rows(spacing/'H_dc_hse_rt.data')
assert len(spacing_rows)==len(a)
assert max(abs(x[j]-y[j]) for x,y in zip(a,spacing_rows) for j in range(13,16)) < 1e-12
spacing_energy=rows(spacing/'H_dc_hse_rt_energy.data')
assert len(spacing_energy)==len(energy)
assert max(abs(x[1]-y[1]) for x,y in zip(energy,spacing_energy)) < 1e-11
print('Grid-spacing input parity passed')

for folder in (coarse,split,spacing):
    cube=(folder/'H_dc_hse_dns_000004.cube').read_text().splitlines()
    natom=abs(int(cube[2].split()[0]))
    charge=sum(float(v) for line in cube[6+natom:] for v in line.split()) # unit grid volume
    assert abs(charge-4)<1e-10, (folder,charge)
print('Electron number conserved')
if args.xc=='hse06':
    metadata={p:p.read_bytes() for p in (root/'gs/data_dcdft/fragments').glob('*/functional.txt')}
    assert len(metadata)==2
    try:
        for p in metadata:p.unlink()
        legacy=run('rt_legacy_metadata',rt,rt=True)
        legacy_rows=rows(legacy/'H_dc_hse_rt.data')
        assert len(legacy_rows)==len(a)
        assert max(abs(x[j]-y[j]) for x,y in zip(a,legacy_rows) for j in range(13,16))<1e-12
    finally:
        for p,data in metadata.items():p.write_bytes(data)
    print('Legacy HSE metadata compatibility passed')
