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


def effective_files(root,spec,*,overlays=()):
    """Resolve verified complete-file preparations without copying or mutations.

    This describes source bytes, not semantic feature compatibility. Each overlay
    is checked against the preceding inventory, including additions and order.
    """
    root=Path(root).resolve()
    files=manifest(root,spec,'source_tree_manifest','source_tree_sha256','salmonpy.source-tree/1')
    tree=path_in(root,spec['source_tree']);check_tree(tree,files)
    resolved={name:{'sha256':expected,'source':path_in(tree,name).relative_to(root).as_posix()}
              for name,expected in files.items()}
    for overlay in overlays:
        replacements=manifest(root,overlay,'replacements','replacements_sha256','salmonpy.source-replacements/1')
        for name,entry in replacements.items():
            path_in(root,name)  # Validate destination names even though no files are written.
            current=resolved.get(name)
            if entry['before_sha256'] is None:
                if current is not None:raise ValueError('new source already exists: '+name)
            elif current is None or current['sha256']!=entry['before_sha256']:
                raise ValueError('replacement base differs: '+name)
            origin=path_in(root,entry['source'])
            if digest(origin)!=entry['after_sha256']:raise ValueError('replacement hash changed: '+name)
            resolved[name]={'sha256':entry['after_sha256'],'source':origin.relative_to(root).as_posix()}
    return resolved


def inventory(root):
    """Current standard and fast-PR effective source matrix, with provenance."""
    root=Path(root).resolve()
    standard=json.loads((root/'source.json').read_text())
    hse=json.loads((root/'hse-rt-state.json').read_text())
    routing=json.loads((root/'hse-explicit-routing.json').read_text())
    if routing['parent_patch_sha256']!=hse['patch_sha256'] or routing['base_commit']!=hse['upstream_commit']:
        raise ValueError('fast HSE routing parent differs')
    profiles={'standard':effective_files(root,standard),
              'fast_hse_reference':effective_files(root,hse),
              'fast_hse':effective_files(root,hse,overlays=(hse,routing))}
    a=profiles['standard'];b=profiles['fast_hse'];common=set(a)&set(b)
    same=sorted(name for name in common if a[name]['sha256']==b[name]['sha256'])
    different=sorted(common-set(same))
    critical=('src/gs/embedded_gs.f90','src/gs/gs_model.f90','src/gs/gs_seed.f90',
              'src/gs/gs_storage.f90','src/io/model_inputs.f90','src/rt/embedded_tddft.f90',
              'src/maxwell/embedded_eh.f90','src/common/force.f90','src/xc/salmon_xc.f90',
              'src/io/checkpoint_restart.f90','src/gs/scf_iteration_dft.f90',
              'src/parallel/communication.f90','src/xc/exx_native.f90','src/xc/exx_symmetry.f90')
    return {'schema':'salmonpy.source-inventory/1','scope':'source byte provenance, not execution/physics certification',
            'counts':{name:len(files) for name,files in profiles.items()},
            'comparison':{'byte_identical':same,'different':different,
                          'standard_only':sorted(set(a)-set(b)),'fast_hse_only':sorted(set(b)-set(a))},
            'critical_paths':{name:{profile:files.get(name) for profile,files in profiles.items()} for name in critical},
            'profiles':profiles}


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description='Inspect verified effective solver source inventories without copying/building.')
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]/'runtime')
    args=parser.parse_args()
    print(json.dumps(inventory(args.root),indent=2))
