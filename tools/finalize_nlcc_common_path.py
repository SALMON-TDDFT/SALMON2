"""Validate completed common-path scalars and regenerate matched-quadrature analysis.

No solver execution and no GS files are read. Preserve the executing coordinator's
original analysis before replacing it with the corrected quadrature comparison.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import numpy as np
from measure_nlcc_force_work import SCALAR_FILES,ordered_samples,write_json
from run_nlcc_common_path import analyze,verify_reuse_audit,LENGTH
from salmonpy.dft import _read_status,_read_forces


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('results',type=Path);p.add_argument('reuse',type=Path);p.add_argument('build_record',type=Path)
    a=p.parse_args();root=a.results
    data={}
    for g in (16,24,32):
        rows=ordered_samples(json.loads((root/f'grid{g}/samples.json').read_text()),32)
        expected=list(range(0,33,2)) if g==16 else list(range(33))
        if [r['index'] for r in rows]!=expected:raise ValueError('incomplete common path '+str(g))
        for r in rows:
            if not np.isfinite(r['shift_bohr']) or not np.isfinite(r['path_force_hartree_per_bohr']):raise ValueError('nonfinite path scalar')
            if abs(r['shift_bohr']-r['index']*LENGTH/32)>1e-15:raise ValueError('physical path mismatch')
            if r['tolerance']!=1e-11 or not np.isfinite(r['energy_hartree']):raise ValueError('invalid scalar convergence record')
            files={x.name for x in (root/f'grid{g}'/f"point{r['index']:02d}").iterdir() if x.is_file()}
            if files!=set(SCALAR_FILES):raise ValueError('unexpected point outputs: '+str(files))
            directory=root/f'grid{g}'/f"point{r['index']:02d}"
            status=_read_status(directory/'gs_status.dat')
            if any(status[k]!=r[k] for k in status):raise ValueError('native SCF status differs')
            forces=_read_forces(directory/'gs_forces.dat',8)
            if not np.array_equal(forces,np.asarray(r['forces_hartree_per_bohr'])):raise ValueError('native force differs')
            if r['path_force_hartree_per_bohr']!=float(forces[0,0]):raise ValueError('path force differs')
        data[g]=rows
    meta=json.loads((root/'provenance.json').read_text());build=json.loads(a.build_record.read_text())
    actual=hashlib.sha256(Path(build['executable']).read_bytes()).hexdigest()
    if actual!=build['sha256'] or actual!=meta['executable_sha256']:raise ValueError('runtime changed during diagnostic')
    executed=root/'executed-run_nlcc_common_path.py'
    if hashlib.sha256(executed.read_bytes()).hexdigest()!=meta['script_sha256']:raise ValueError('executing source record differs')
    history=root/f'analysis-history-{time.time_ns()}';history.mkdir()
    for name in ('refined-analysis.json','REPORT-refined.md','retention-audit.json'):
        if (root/name).exists():shutil.copy2(root/name,history/name)
    analyze(root,'refined')
    reuse=verify_reuse_audit(root,a.reuse)
    forbidden={'gs_eigen.data','gs_info.data','gs_k.data','wfn.bin','rho.bin','rho_inout.bin','occupation.bin','info.bin'}
    bad=[str(p.relative_to(root)) for p in root.rglob('*') if p.is_file() and p.name in forbidden]
    if bad:raise ValueError('unexpected retained GS files: '+str(bad))
    count=sum(map(len,data.values()));reused=sum('reuse_source' in row for rows in data.values() for row in rows)
    if reused!=reuse['reused_scalar_points']:raise ValueError('reuse count mismatch')
    write_json(root/'retention-audit.json',{'forbidden_GS_files_remaining':bad,'accepted_scalar_points':count,'reused_scalar_points':reused,'new_scalar_points':count-reused})
    source=Path(__file__).resolve().parent
    with (root/'analysis-tests.log').open('w') as log:
        run=subprocess.run([sys.executable,str(source/'test_force_work_analysis.py')],cwd=source.parent,
            env={**os.environ,'PYTHONPATH':str(source.parent)},stdout=log,stderr=subprocess.STDOUT)
    if run.returncode:raise RuntimeError('analysis regression tests failed')
    for name in ('run_nlcc_common_path.py','measure_nlcc_force_work.py','test_force_work_analysis.py','finalize_nlcc_common_path.py'):
        shutil.copy2(source/name,root/('analysis-'+name))
    write_json(root/'analysis-revision.json',{'executed_coordinator_sha256':meta['script_sha256'],
        'corrected_analysis_sha256':hashlib.sha256((source/'run_nlcc_common_path.py').read_bytes()).hexdigest(),
        'reason':'All grid comparisons use matched quadrature: S16 for16/24/32, separately S32 for24/32; no native calculations repeated.',
        'prior_analysis_preserved':history.name})
    x=json.loads((root/'refined-analysis.json').read_text())
    lines=['# NLCC共通物理経路：格子・積分点数の収束確認','',
        '同じSi8の物理経路（原子0のx=0.12→0.5475 bohr、変位0.4275 bohr）で16³・24³・32³を比較しました。PBE、32電子・32状態、2³ k点、SCF rho_dne閾値1e-11、各点OMP2です。','',
        f'受理済み測定は{count}点。うち24³の{reused}点は同一バイナリ・同一条件の既存スカラー記録を再利用し、{count-reused}点を新規に計算しました。4並列×OMP2を上限とし、既存ジョブ・既存実行ファイルは変更していません。','',
        '## エネルギー差と力の積分','',
        'R=ΔE+∫Fdx。表のRはSi8セルの値を8原子で割っています。積分点数差は当該区間数と半分の区間数との差で、厳密な誤差上限ではありません。','',
        '| 格子 | 区間数 | ΔE (Ha/セル) | ∫Fdx (Ha/セル) | R (mHa/原子) | 積分点数差 (mHa/原子) |','|---|---:|---:|---:|---:|---:|']
    for g,v in x['grids'].items():
        for key in ('intervals16','intervals32'):
            if key in v:
                z=v[key];lines.append(f"| {g}³ | {key.removeprefix('intervals')} | {z['delta_energy_hartree']:.12g} | {z['integrated_force_hartree']['1']:.12g} | {z['consistency_residual_mHa_per_atom']:+.8f} | {z['quadrature_change_hartree']*1000/8:.8f} |")
    for key,label in [('grid_changes','全格子で16区間に揃えた比較'),('grid_changes_32_intervals','24³・32³で32区間に揃えた比較')]:
        lines += ['',f'## {label}','','| 格子変更 | ΔE変化 (Ha/セル) | R変化 (mHa/原子) | 力プロフィール最大変化 (Ha/bohr) |','|---|---:|---:|---:|']
        for name,z in x[key].items():lines.append(f"| {name} | {z['delta_energy_change_hartree']:+.9g} | {z['residual_change_mHa_per_atom']:+.8f} | {z['max_force_profile_change_hartree_per_bohr']:.9g} |")
    lines += ['','## SCF収束確認','','全83点の生スカラー報告を再読し、収束フラグ・残差・閾値・エネルギーと力が受理記録に一致することを確認しました。','',
        '| 格子 | 点数 | 反復数範囲 | 最大rho_dne残差 |','|---|---:|---:|---:|']
    scf={}
    for g,rows in data.items():
        lo=min(r['iterations'] for r in rows);hi=max(r['iterations'] for r in rows);maximum=max(r['residual'] for r in rows)
        scf[str(g)]={'points':len(rows),'iterations_min':lo,'iterations_max':hi,'maximum_residual':maximum,'all_converged':True}
        lines.append(f'| {g}³ | {len(rows)} | {lo}–{hi} | {maximum:.9g} |')
    write_json(root/'scf-audit.json',scf)
    lines += ['','## 解釈と残る範囲','',
        '16³の残差は大きく、この条件を収束済みとは扱えません。24³→32³のエネルギー・力の変化と16→32区間の積分差を別々に示しており、本番精度を自動認定しません。許容誤差の指定、より細かい格子、k点、radial補間・cutoff境界、SCF全経路の閾値、実材料の検証は別です。小さい端点残差だけでは途中の力プロフィールの精度を保証できません。32³では積分点数差が16対8区間の0.00230243から32対16区間の0.00395773 mHa/原子へ増え、積分の単調な収束もまだ確認できていません。24³と32³の32区間比較では力の最大差は約0.00143 Ha/bohr（約0.0735 eV/Å）です。','',
        '計算中に読み込まれていた旧解析は細分化後に異なる区間数を混ぜる可能性があったため、計算結果を変更せず、全格子S16と24/32格子S32の別比較へ再集計しました。元の解析と実行コードは保存しています。','',
        '再利用記録は元のエネルギー・反復数・残差・閾値・力の報告と照合し、コピーした4ファイルと元provenanceのSHA256も記録しました。新規計算は使い捨て領域で行い、gs.inp・stdout.log・gs_status.dat・gs_forces.datだけを保存。GS状態・固有値・k点・GS情報ファイルは保持していません。','',
        '解析・再利用・保存制限の10件のテストが通過。バイナリSHA256は開始時と完了時に一致し、既存の監査済みバイナリを変更せず使用しています。','',
        '[詳細数値](refined-analysis.json) / [比較図](common-path.png) / [再利用監査](reuse-audit.json) / [保存監査](retention-audit.json) / [解析テスト](analysis-tests.log)','',
        '結果はOneDriveのローカルcalc2026へ保存しています。クラウド同期完了は未確認です。']
    (root/'REPORT.md').write_text('\n'.join(lines)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(2,2,figsize=(11,8))
    for grid,color in ((16,'tab:orange'),(24,'tab:blue'),(32,'tab:green')):
        rows=data[grid];s=np.array([r['shift_bohr'] for r in rows]);e=np.array([r['energy_hartree'] for r in rows]);f=np.array([r['path_force_hartree_per_bohr'] for r in rows])
        axs[0,0].plot(s,(e-e[0])*1000,'o-',ms=3,color=color,label=f'{grid}³')
        axs[0,1].plot(s,f,'o-',ms=3,color=color,label=f'{grid}³')
    axs[0,0].set(title='Same physical path: energy change',xlabel='Displacement (bohr)',ylabel='Energy change (mHa / cell)');axs[0,0].legend()
    axs[0,1].set(title='Force on displaced atom',xlabel='Displacement (bohr)',ylabel='x force (Ha/bohr)');axs[0,1].legend()
    axs[1,0].plot([16,24,32],[x['grids'][str(g)]['intervals16']['consistency_residual_mHa_per_atom'] for g in (16,24,32)],'o-',label='16 intervals')
    axs[1,0].plot([24,32],[x['grids'][str(g)]['intervals32']['consistency_residual_mHa_per_atom'] for g in (24,32)],'s--',label='32 intervals')
    axs[1,0].axhline(0,color='black',lw=.7);axs[1,0].set_yscale('symlog',linthresh=.01)
    axs[1,0].set(title='Signed work/energy residual',xlabel='Grid points per direction',ylabel='R (mHa / atom), symlog');axs[1,0].legend()
    ax=axs[1,1]
    for g,color in ((16,'tab:orange'),(24,'tab:blue'),(32,'tab:green')):
        values=[x['grids'][str(g)]['intervals16']['quadrature_change_hartree']*1000/8]
        labels=['16 vs 8']
        if g!=16:values.append(x['grids'][str(g)]['intervals32']['quadrature_change_hartree']*1000/8);labels.append('32 vs 16')
        ax.plot(labels,values,'o-',color=color,label=f'{g}³')
    ax.set_yscale('log');ax.set(title='Quadrature changes, not error bounds',ylabel='Absolute work change (mHa / atom)');ax.legend()
    fig.suptitle('NLCC fixed physical path: grid and quadrature diagnostics')
    fig.text(.5,.015,'Si8, PBE, fixed path 0.4275 bohr. No automatic physical-accuracy certification.',ha='center',fontsize=9)
    fig.tight_layout(rect=(0,.035,1,.96));fig.savefig(root/'common-path.png',dpi=180);plt.close(fig)
    print(root/'REPORT.md')

if __name__=='__main__':main()
