"""Prepare the optional full SALMON runtime from a pinned external checkout.

python tools/prepare_salmon_runtime.py UPSTREAM_CHECKOUT NEW_DIRECTORY
No download, build, branch change, or modification of the input checkout occurs.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import tempfile


def prepare(source, directory):
    root=Path(__file__).resolve().parents[1]
    spec=json.loads((root/'runtime/source.json').read_text())
    patch=root/'runtime'/spec['patch']
    if hashlib.sha256(patch.read_bytes()).hexdigest()!=spec['patch_sha256']:
        raise ValueError('runtime patch checksum mismatch')
    source=Path(source).expanduser().resolve()
    commit=spec['upstream_commit']
    found=subprocess.check_output(['git','-C',str(source),'rev-parse',commit+'^{commit}'],text=True).strip()
    if found!=commit:raise ValueError('upstream revision mismatch')
    dest=Path(directory).expanduser().resolve()
    if dest.is_relative_to(source):raise ValueError('destination must be outside the upstream checkout')
    if dest.is_relative_to(root):raise ValueError('destination must be outside the SALMONpy checkout')
    dest.mkdir(parents=True,exist_ok=False)
    try:
        with tempfile.TemporaryFile() as archive:
            subprocess.run(['git','-C',str(source),'archive','--format=tar',commit],stdout=archive,check=True)
            archive.seek(0)
            with tarfile.open(fileobj=archive,mode='r:') as tar:
                for member in tar:
                    path=PurePosixPath(member.name)
                    if path.is_absolute() or '..' in path.parts or '\\' in member.name:
                        raise ValueError('unsafe upstream archive path')
                    target=dest.joinpath(*path.parts)
                    if member.isdir():target.mkdir(parents=True,exist_ok=True)
                    elif member.isfile():
                        target.parent.mkdir(parents=True,exist_ok=True)
                        with tar.extractfile(member) as src,target.open('wb') as out:shutil.copyfileobj(src,out)
                        target.chmod(member.mode & 0o777)
                    else:raise ValueError('unsupported upstream archive member: '+member.name)
        # Apply to this exported tree, never to a containing checkout/index.
        env={k:v for k,v in os.environ.items() if k not in
             ('GIT_DIR','GIT_WORK_TREE','GIT_INDEX_FILE','GIT_COMMON_DIR')}
        env['GIT_CEILING_DIRECTORIES']=str(dest.parent)
        subprocess.run(['git','apply','--check',str(patch)],cwd=dest,env=env,check=True)
        subprocess.run(['git','apply','--whitespace=nowarn',str(patch)],cwd=dest,env=env,check=True)
        for local,external in spec['shared_kernels'].items():
            if (root/local).read_bytes()!=(dest/external).read_bytes():
                raise ValueError('shared kernel differs from runtime patch: '+local)
        subprocess.run([sys.executable,
                        str(root/'tools/check_forbidden_eigensolvers.py'),'--source-root',str(dest/'src')],check=True)
        (dest/'salmonpy-runtime-source.json').write_text(json.dumps(spec,indent=2)+'\n')
        return dest
    except BaseException:
        shutil.rmtree(dest)
        raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('upstream');parser.add_argument('directory')
    args=parser.parse_args()
    print(prepare(args.upstream,args.directory))
