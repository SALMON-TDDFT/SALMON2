"""Apply the explicit RT step-config overlay to a NEW external runtime snapshot.

python tools/prepare_rt_step_config.py PREPARED_RUNTIME NEW_DIRECTORY
Opt-in experiment: does not enable full native multi-context TDDFT. Never edits
its input tree, downloads, builds, or replaces an executable. Failed snapshots
are preserved with their preparation status for inspection, not approved reuse.
"""
import argparse,hashlib,json,os,shutil,subprocess,sys
from pathlib import Path


def prepare(source,directory):
    root=Path(__file__).resolve().parents[1]
    spec=json.loads((root/'runtime/rt-step-config.json').read_text())
    patch=root/'runtime'/spec['patch']
    if hashlib.sha256(patch.read_bytes()).hexdigest()!=spec['patch_sha256']:raise ValueError('overlay hash changed')
    source=Path(source).resolve();dest=Path(directory).resolve()
    if dest.is_relative_to(source) or dest.is_relative_to(root):raise ValueError('new external output required')
    base=json.loads((source/'salmonpy-runtime-source.json').read_text())
    if base.get('upstream_commit')!=spec['upstream_commit'] or base.get('patch_sha256')!=spec['base_patch_sha256']:
        raise ValueError('explicit overlay requires its pinned base profile')
    if (source/'salmonpy-rt-step-config.json').exists():raise ValueError('source already has step-config overlay')
    checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    subprocess.run([sys.executable,str(checker),'--source-root',str(source/'src')],check=True)
    # Claim only a new output; record failure even if copying is interrupted.
    dest.mkdir(parents=True,exist_ok=False)
    status={**spec,'status':'preparing'}
    marker=dest/'salmonpy-rt-step-config.json'
    marker.write_text(json.dumps(status,indent=2)+'\n')
    try:
        shutil.copytree(source,dest,dirs_exist_ok=True)
        env={k:v for k,v in os.environ.items() if k not in ('GIT_DIR','GIT_WORK_TREE','GIT_INDEX_FILE','GIT_COMMON_DIR')}
        env['GIT_CEILING_DIRECTORIES']=str(dest.parent)
        subprocess.run(['git','apply','--check',str(patch)],cwd=dest,env=env,check=True)
        subprocess.run(['git','apply','--whitespace=nowarn',str(patch)],cwd=dest,env=env,check=True)
        subprocess.run([sys.executable,str(checker),'--source-root',str(dest/'src')],check=True)
        status['status']='prepared-and-source-audited'
        status['modified_source_hashes']={name:hashlib.sha256((dest/name).read_bytes()).hexdigest() for name in spec['modified_sources']}
        marker.write_text(json.dumps(status,indent=2)+'\n')
        return dest
    except BaseException:
        status['status']='failed-not-approved';marker.write_text(json.dumps(status,indent=2)+'\n');raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('output')
    a=p.parse_args();print(prepare(a.source,a.output))
