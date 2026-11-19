"""Prepare a fresh pinned fast-PR HSE snapshot with optional RT-state overlay.

python tools/prepare_hse_rt_state.py NEW_OUTPUT [--reference-only]
No download/build/PR edits. Source audit precedes every preparation/build; local
CMake guards retain mandatory checks even though the PR excludes those guards.
Failed evidence is preserved and explicitly not approved for execution.
"""
import argparse,hashlib,json,subprocess,sys
from pathlib import Path

import importlib.util
_store_spec=importlib.util.spec_from_file_location('salmonpy_source_store',Path(__file__).with_name('runtime_sources.py'))
store=importlib.util.module_from_spec(_store_spec);_store_spec.loader.exec_module(store)


def prepare(source,directory=None,*,reference_only=False,explicit_routing=False):
    if directory is None:directory=source;source=None
    root=Path(__file__).resolve().parents[1];spec=json.loads((root/'runtime/hse-rt-state.json').read_text())
    checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    routing=None
    if explicit_routing:
        if reference_only:raise ValueError('explicit routing cannot use reference-only preparation')
        routing=json.loads((root/'runtime/hse-explicit-routing.json').read_text())
        if routing['parent_patch_sha256']!=spec['patch_sha256'] or routing['base_commit']!=spec['upstream_commit']:
            raise ValueError('explicit routing requires its pinned parent overlay')
    dest=Path(directory).resolve()
    if (source is not None and dest.is_relative_to(Path(source).resolve())) or dest.is_relative_to(root):raise ValueError('new external directory required')
    dest.mkdir(parents=True,exist_ok=False);marker=dest/'salmonpy-hse-rt-state.json'
    record={**spec,'status':'preparing','overlay_applied':not reference_only,'explicit_rt_routing':explicit_routing}
    marker.write_text(json.dumps(record,indent=2)+'\n')
    try:
        store.copy_tree(root/'runtime',spec,dest)
        subprocess.run([sys.executable,str(checker),'--source-root',str(dest/'src')],check=True)
        for name,expected in spec['reference_source_hashes'].items():
            if hashlib.sha256((dest/name).read_bytes()).hexdigest()!=expected:raise ValueError('pinned source changed: '+name)
        if not reference_only:
            store.replace_files(root/'runtime',spec,dest)
            if routing is not None:
                store.replace_files(root/'runtime',routing,dest)
                record['routing_overlay']=routing
                store.check_shared_kernels(root,dest,routing.get('shared_kernels',{}))
        guard=dest/'cmakefiles/salmonpy_prohibition.cmake'
        subprocess.run([sys.executable,str(checker),'--source-root',str(dest/'src')],check=True)
        record.update(status='prepared-and-source-audited',checker_sha256=hashlib.sha256(checker.read_bytes()).hexdigest(),
            modified_source_hashes={name:hashlib.sha256((dest/name).read_bytes()).hexdigest() for name in spec['modified_sources'] if (dest/name).is_file()},
            cmake_guard_sha256=hashlib.sha256(guard.read_bytes()).hexdigest())
        marker.write_text(json.dumps(record,indent=2)+'\n');return dest
    except BaseException:
        record['status']='failed-not-approved';marker.write_text(json.dumps(record,indent=2)+'\n');raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('output',nargs='?');p.add_argument('--reference-only',action='store_true');p.add_argument('--explicit-routing',action='store_true')
    a=p.parse_args();print(prepare(a.source,a.output,reference_only=a.reference_only,explicit_routing=a.explicit_routing))
