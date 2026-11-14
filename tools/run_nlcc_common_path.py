"""Staged fixed-physical-path NLCC grid/quadrature diagnostic.

Each point has disposable native scratch; only input/energy/force/convergence
reports are exported. No solver build and no GS state/eigen/k/info collection.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import shutil
import time
import numpy as np
from salmonpy.dft import DFT,Pseudopotential,SCFSettings,_read_status,_read_forces
from measure_nlcc_force_work import SCALAR_FILES,run_scalar_only,ordered_samples,write_json,summarize

PP_SHA='39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282'
LENGTH=.4275
INTERVALS=32


def model_for(grid,shift,pseudo):
    base=np.array([[0,0,0],[0,.5,.5],[.5,0,.5],[.5,.5,0],
                   [.25,.25,.25],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]],float)
    base[0,0]+=(.12+shift)/10.26
    return DFT((14,)*8,np.eye(3)*10.26,base,(pseudo,),32,32,(grid,)*3,(2,)*3,'libxc_pbe')


def evaluate(exe,pseudo,root,grid,index):
    shift=float(np.linspace(0,LENGTH,INTERVALS+1)[index])
    model=model_for(grid,shift,pseudo)
    text=model.input_text(SCFSettings(max_iterations=600,tolerance=1e-11)).replace("write_gs_restart_data='all'","write_gs_restart_data='no'")
    dest=root/f'grid{grid}'/f'point{index:02d}'
    if dest.exists():
        backup=dest/f'interrupted-{time.time_ns()}';backup.mkdir()
        for name in SCALAR_FILES:
            if (dest/name).exists():shutil.move(dest/name,backup/name)
    env={key:value for key,value in os.environ.items() if key!='SALMONPY_TDDFT_WORKER'}
    env.update(OMP_NUM_THREADS='2',OPENBLAS_NUM_THREADS='1')
    start=time.monotonic()
    code=run_scalar_only(exe,text,pseudo.data,dest,env)
    if code or 'end SALMON' not in (dest/'stdout.log').read_text():raise RuntimeError(f'failed SCF {grid}/{index}')
    status=_read_status(dest/'gs_status.dat');forces=_read_forces(dest/'gs_forces.dat',8)
    return dict(index=index,shift_bohr=shift,**status,seconds=time.monotonic()-start,
                forces_hartree_per_bohr=forces.tolist(),path_force_hartree_per_bohr=float(forces[0,0]))


def verify_reused_row(source, copied, row):
    hashes={}
    for name in SCALAR_FILES:
        original=(source/name).read_bytes();actual=(copied/name).read_bytes()
        hashes[name]=hashlib.sha256(original).hexdigest()
        if hashlib.sha256(actual).hexdigest()!=hashes[name]:
            raise ValueError('copied scalar report differs: '+name)
    status=_read_status(source/'gs_status.dat')
    force=_read_forces(source/'gs_forces.dat',8)
    for key in ('energy_hartree','iterations','residual','tolerance'):
        if status[key]!=row.get(key):raise ValueError('reused '+key+' does not match scalar native status')
    if status['tolerance']!=1e-11:raise ValueError('reused SCF tolerance differs from common path conditions')
    if (not np.array_equal(force,np.asarray(row['forces_hartree_per_bohr']))
            or force[0,0]!=row['path_force_hartree_per_bohr']):
        raise ValueError('reused force does not match scalar native report')
    return hashes


def verify_reuse_audit(root, source):
    provenance_hash=hashlib.sha256((source/'provenance.json').read_bytes()).hexdigest()
    rows=ordered_samples(json.loads((root/'grid24/samples.json').read_text()),INTERVALS)
    records=[]
    for row in rows:
        if 'reuse_source' not in row: continue
        src=Path(row['reuse_source']);dest=root/'grid24'/f"point{row['index']:02d}"
        if src.parent!=(source/'grid24') or row['index']%2:
            raise ValueError('invalid reuse provenance')
        hashes=verify_reused_row(src,dest,row)
        records.append({'index':row['index'],'source':str(src),'scalar_hashes':hashes})
    result={'source_provenance_sha256':provenance_hash,'reused_scalar_points':len(records),'records':records}
    write_json(root/'reuse-audit.json',result)
    return result


def seed_existing(root,source,meta):
    old=json.loads((source/'provenance.json').read_text())
    checks={'executable_sha256':meta['executable_sha256'],'effective_pseudo_sha256':PP_SHA,
            'xc':'libxc_pbe','mode':'nlcc','rigid':False,'threads':2,'intervals':16,
            'tolerance':1e-11,'lattice_bohr':10.26,'states':32,'kmesh':[2,2,2]}
    if any(old.get(k)!=v for k,v in checks.items()):raise ValueError('reuse source settings mismatch')
    rows=ordered_samples(json.loads((source/'grid24/samples.json').read_text()),16)
    if [r['index'] for r in rows]!=list(range(17)):raise ValueError('incomplete reuse path')
    shifts=np.linspace(0,LENGTH,17)
    if any(abs(r['shift_bohr']-x)>1e-15 for r,x in zip(rows,shifts)):raise ValueError('reuse physical path differs')
    seeded=[]
    for row in rows:
        new=dict(row);new['index']=2*row['index'];new['reuse_source']=str(source/'grid24'/f"point{row['index']:02d}")
        src=source/'grid24'/f"point{row['index']:02d}";dest=root/'grid24'/f"point{new['index']:02d}";dest.mkdir(parents=True)
        # Whitelist, never copy an entire native result directory.
        for name in SCALAR_FILES:
            if not (src/name).is_file():raise ValueError('missing reused scalar report '+name)
            shutil.copy2(src/name,dest/name)
        new['reused_scalar_hashes']=verify_reused_row(src,dest,new)
        new['source_provenance_sha256']=hashlib.sha256((source/'provenance.json').read_bytes()).hexdigest()
        seeded.append(new)
    write_json(root/'grid24/samples.json',seeded)


def grid_comparisons(stats, profiles, intervals):
    changes={}
    for lo,hi in ((16,24),(24,32)):
        key=f'intervals{intervals}'
        if key not in stats[str(lo)] or key not in stats[str(hi)]: continue
        sx=stats[str(lo)][key];sy=stats[str(hi)][key]
        changes[f'{lo}-to-{hi}']=dict(quadrature_intervals=intervals,
            delta_energy_change_hartree=sy['delta_energy_hartree']-sx['delta_energy_hartree'],
            residual_change_mHa_per_atom=sy['consistency_residual_mHa_per_atom']-sx['consistency_residual_mHa_per_atom'],
            max_force_profile_change_hartree_per_bohr=float(np.max(np.abs(profiles[hi]-profiles[lo]))))
    return changes


def analyze(root,stage):
    stats={};profiles={};fine_profiles={}
    for grid in (16,24,32):
        rows=ordered_samples(json.loads((root/f'grid{grid}/samples.json').read_text()),32)
        lookup={r['index']:r for r in rows}
        coarse=[lookup[i] for i in range(0,33,2)]
        profile=np.array([r['forces_hartree_per_bohr'] for r in coarse]);profiles[grid]=profile
        entry={'intervals16':summarize(coarse,LENGTH),'accepted_points':len(rows)}
        if len(rows)==33:
            entry['intervals32']=summarize(rows,LENGTH)
            fine_profiles[grid]=np.array([r['forces_hartree_per_bohr'] for r in rows])
        stats[str(grid)]=entry
    changes=grid_comparisons(stats,profiles,16)
    fine_changes=grid_comparisons(stats,fine_profiles,32)
    result=dict(stage=stage,path_length_bohr=LENGTH,start_atom_x_bohr=.12,end_atom_x_bohr=.12+LENGTH,
                grids=stats,grid_changes=changes,grid_changes_32_intervals=fine_changes,
                scope='Measured grid/quadrature changes, not automatic production convergence certification; k-point/material validation separate.')
    write_json(root/f'{stage}-analysis.json',result)
    lines=[f'# NLCC共通経路の収束診断（{stage}）','',
           '全格子で同じ物理経路：原子0のx=0.12→0.5475 bohr（長さ0.4275）。Si8、PBE、2³ k点、32状態、SCF閾値1e-11、OMP2。','',
           '| 格子 | 区間数 | ΔE (Ha) | R=ΔE+∫Fdx (mHa/原子) | 積分点数差 (mHa/原子) |',
           '|---|---:|---:|---:|---:|']
    for grid,entry in stats.items():
        for label in ('intervals16','intervals32'):
            if label in entry:
                s=entry[label];lines.append(f"| {grid}³ | {label.removeprefix('intervals')} | {s['delta_energy_hartree']:.12g} | {s['consistency_residual_mHa_per_atom']:+.8f} | {s['quadrature_change_hartree']*1000/8:.8f} |")
    lines+=['','積分点数差は当該区間数と半分の区間数との差で、厳密な誤差上限ではありません。格子変更は同じ経路で比較しています。','',
            '| 格子変更（全て16区間） | ΔEの変化 (Ha) | Rの変化 (mHa/原子) | 共通17点の力の最大変化 (Ha/bohr) |','|---|---:|---:|---:|']
    for name,c in changes.items():lines.append(f"| {name} | {c['delta_energy_change_hartree']:+.8g} | {c['residual_change_mHa_per_atom']:+.8f} | {c['max_force_profile_change_hartree_per_bohr']:.8g} |")
    if fine_changes:
        lines+=['','32区間で揃えた24³→32³の比較（共通33点）：']
        c=fine_changes['24-to-32']
        lines.append(f"- ΔE変化 {c['delta_energy_change_hartree']:+.8g} Ha、R変化 {c['residual_change_mHa_per_atom']:+.8f} mHa/原子、力の最大変化 {c['max_force_profile_change_hartree_per_bohr']:.8g} Ha/bohr。")
    lines+=['','数値差を示す診断で、本番の収束を自動認定しません。k点・実材料の検証、radial補間・cutoff境界、SCF全経路の閾値収束は別です。',
            'GS状態・固有値・k点・GS情報は保存せず、必要な4ファイルだけを一時領域から出しています。OneDriveローカル保存で、クラウド同期は未確認です。']
    (root/f'REPORT-{stage}.md').write_text('\n'.join(lines)+'\n')
    print(json.dumps({'stage':stage,'grid_changes':changes}),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('build_record',type=Path);p.add_argument('pseudo',type=Path);p.add_argument('reuse',type=Path);p.add_argument('output',type=Path)
    p.add_argument('--workers',type=int,default=4);p.add_argument('--resume',action='store_true');a=p.parse_args()
    if not 1<=a.workers<=4:raise ValueError('at most four workers, each uses two OMP threads')
    b=json.loads(a.build_record.read_text());exe=Path(b['executable'])
    if b.get('status')!='built-and-audited' or hashlib.sha256(exe.read_bytes()).hexdigest()!=b['sha256']:raise ValueError('audited binary checksum mismatch')
    raw=a.pseudo.read_bytes()
    if hashlib.sha256(raw).hexdigest()!=PP_SHA:raise ValueError('UPF checksum mismatch')
    pseudo=Pseudopotential(14,2,'.upf',raw)
    meta=dict(executable_sha256=b['sha256'],original_pseudo_sha256=PP_SHA,path_length_bohr=LENGTH,
              intervals=32,grids=[16,24,32],tolerance=1e-11,threads=2,workers=a.workers,
              lattice_bohr=10.26,states=32,kmesh=[2,2,2],xc='libxc_pbe',start_atom_x_bohr=.12,
              script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if a.resume:
        old=json.loads((a.output/'provenance.json').read_text())
        for k in ('executable_sha256','original_pseudo_sha256','path_length_bohr','intervals','grids','tolerance','threads','lattice_bohr','states','kmesh','xc','start_atom_x_bohr'):
            if old.get(k)!=meta[k]:raise ValueError('resume settings differ: '+k)
        write_json(a.output/f'resume-{time.time_ns()}.json',meta)
    else:
        a.output.mkdir(parents=True,exist_ok=False);write_json(a.output/'provenance.json',meta)
        seed_existing(a.output,a.reuse,meta)
    rows={}
    for g in (16,24,32):
        path=a.output/f'grid{g}/samples.json'
        rows[g]=ordered_samples(json.loads(path.read_text()),32) if path.exists() else []
    for stage,tasks in [('coarse',[(g,i) for i in range(0,33,2) for g in (32,16,24)]),
                        ('refined',[(g,i) for i in range(1,32,2) for g in (32,24)])]:
        tasks=[(g,i) for g,i in tasks if i not in {r['index'] for r in rows[g]}]
        with ThreadPoolExecutor(max_workers=a.workers) as pool:
            futures={pool.submit(evaluate,exe,pseudo,a.output,g,i):(g,i) for g,i in tasks}
            for future in as_completed(futures):
                g,i=futures[future]
                row=future.result();rows[g].append(row);rows[g]=ordered_samples(rows[g],32)
                write_json(a.output/f'grid{g}/samples.json',rows[g])
                print(json.dumps(dict(stage=stage,grid=g,index=i,seconds=row['seconds'],iterations=row['iterations'])),flush=True)
        analyze(a.output,stage)
    audit=verify_reuse_audit(a.output,a.reuse)
    forbidden={'gs_eigen.data','gs_k.data','gs_info.data','wfn.bin','rho.bin','rho_inout.bin','occupation.bin','info.bin'}
    bad=[str(p) for p in a.output.rglob('*') if p.is_file() and p.name in forbidden]
    if bad:raise RuntimeError('unexpected retained GS files: '+str(bad))
    write_json(a.output/'retention-audit.json',{'forbidden_GS_files_remaining':bad,'accepted_scalar_points':sum(map(len,rows.values())),'reused_scalar_points':audit['reused_scalar_points']})

if __name__=='__main__':main()
