"""Verify/copy repository-owned solver sources; never apply patches or use Git."""
import hashlib
from importlib.util import source_from_cache
import json
from pathlib import Path,PurePosixPath
import shutil

def digest(path):
    path=Path(path)
    if path.is_symlink() or not path.is_file():raise ValueError('regular source file required')
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):h.update(block)
    return h.hexdigest()

def path_in(root,name):
    if not isinstance(name,str) or not name or '\\' in name or ':' in name:
        raise ValueError('invalid source path')
    parts=PurePosixPath(name)
    if parts.is_absolute() or '..' in parts.parts or '.'==name:
        raise ValueError('source path escapes root')
    root=Path(root).resolve();path=root
    for part in parts.parts:
        path=path/part
        if path.is_symlink():raise ValueError('source symlinks are not supported')
    if not path.resolve().is_relative_to(root):raise ValueError('source path escapes root')
    return path

def manifest(root,spec,key,hash_key,schema):
    path=path_in(root,spec[key])
    if digest(path)!=spec[hash_key]:raise ValueError('source manifest hash changed')
    data=json.loads(path.read_text())
    if data.get('schema')!=schema or not isinstance(data.get('files'),dict) or not data['files']:
        raise ValueError('invalid source inventory')
    return data['files']

def check_tree(tree,files,*,extra=()):
    actual=set()
    for path in tree.rglob('*'):
        if path.is_symlink():raise ValueError('source symlinks are not supported')
        if path.is_file():
            # pip compiles bundled .py resources. Accept caches only for listed
            # source files; never copy or execute their bytecode in prepared trees.
            if path.suffix=='.pyc' and path.parent.name=='__pycache__':
                try: original=Path(source_from_cache(str(path))).relative_to(tree).as_posix()
                except (ValueError,TypeError): original=None
                if original in files:continue
            actual.add(path.relative_to(tree).as_posix())
        elif not path.is_dir():raise ValueError('special source entries are not supported')
    if actual-set(extra)!=set(files):raise ValueError('source inventory differs: missing/extra files')
    for name,expected in files.items():
        if digest(path_in(tree,name))!=expected:raise ValueError('source hash changed: '+name)

def copy_tree(root,spec,dest):
    files=manifest(root,spec,'source_tree_manifest','source_tree_sha256','salmonpy.source-tree/1')
    source=path_in(root,spec['source_tree']);check_tree(source,files)
    dest=Path(dest).resolve()
    if dest.is_relative_to(source) or source.is_relative_to(dest):raise ValueError('overlapping source/destination')
    # Destination may contain only the caller's preparation marker.
    for name in files:
        target=path_in(dest,name)
        if target.exists():raise ValueError('destination source already exists: '+name)
    copy_files(source,dest,files)
    for name,expected in files.items():
        if digest(path_in(dest,name))!=expected:raise ValueError('source changed while copying: '+name)
    return files

def replace_files(root,spec,dest):
    files=manifest(root,spec,'replacements','replacements_sha256','salmonpy.source-replacements/1')
    pending=[]
    for name,entry in files.items():
        source=path_in(root,entry['source']);target=path_in(dest,name)
        if digest(source)!=entry['after_sha256']:raise ValueError('replacement hash changed: '+name)
        if entry['before_sha256'] is None:
            if target.exists():raise ValueError('new source already exists: '+name)
        elif not target.is_file() or digest(target)!=entry['before_sha256']:
            raise ValueError('replacement base differs: '+name)
        pending.append((source,target,entry['after_sha256']))
    # Preflight all files before changing any; caller retains failed output on I/O failure.
    for source,target,expected in pending:
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,target)
        if digest(target)!=expected:raise ValueError('replacement changed while copying: '+str(target))


def copy_files(source,dest,files):
    for name,expected in files.items():
        origin=path_in(source,name);target=path_in(dest,name)
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(origin,target)
        if digest(target)!=expected:raise ValueError('source changed while copying: '+name)
