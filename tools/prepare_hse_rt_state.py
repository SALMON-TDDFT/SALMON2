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


def prepare(directory,*,reference_only=False,explicit_routing=False):
    root=Path(__file__).resolve().parents[1]
    if not reference_only:
        entry=importlib.util.spec_from_file_location('salmonpy_current_prepare',Path(__file__).with_name('prepare_salmon_runtime.py'))
        module=importlib.util.module_from_spec(entry);entry.loader.exec_module(module)
        module.store=store  # Shared preparation boundary, including failure injection/audits.
        return module.prepare(directory)
    if explicit_routing:raise ValueError('explicit routing cannot use reference-only preparation')
    spec=json.loads((root/'runtime/source-history/20261119/hse-fast-reference.json').read_text())
    if not (root/'runtime'/spec['source_tree']).is_file():
        raise FileNotFoundError('historical reference archive is not bundled; use the source checkout for developer comparisons')
    dest=Path(directory).resolve()
    if dest.is_relative_to(root):
        raise ValueError('new external directory required')
    dest.mkdir(parents=True,exist_ok=False);marker=dest/'salmonpy-hse-rt-state.json'
    record={**spec,'status':'preparing','overlay_applied':False,'explicit_rt_routing':False}
    try:
        store.copy_tree(root/'runtime',spec,dest)
        subprocess.run([sys.executable,str(Path(__file__).with_name('check_forbidden_eigensolvers.py')),'--source-root',str(dest/'src')],check=True)
        record['status']='prepared-and-source-audited'
        marker.write_text(json.dumps(record,indent=2)+'\n');return dest
    except BaseException:
        record['status']='failed-not-approved';marker.write_text(json.dumps(record,indent=2)+'\n');raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory');p.add_argument('--reference-only',action='store_true');p.add_argument('--explicit-routing',action='store_true')
    a=p.parse_args();print(prepare(a.directory,reference_only=a.reference_only,explicit_routing=a.explicit_routing))
