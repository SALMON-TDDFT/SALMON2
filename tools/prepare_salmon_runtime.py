"""Prepare the bundled solver in a new external directory, without Git or patches.

The single positional argument is the output directory.
The pinned repository-owned source is used.
"""
import argparse,json,subprocess,sys,hashlib
from pathlib import Path

import importlib.util
_store_spec=importlib.util.spec_from_file_location('salmonpy_source_store',Path(__file__).with_name('runtime_sources.py'))
store=importlib.util.module_from_spec(_store_spec);_store_spec.loader.exec_module(store)

def prepare(directory):
    root=Path(__file__).resolve().parents[1]
    spec=json.loads((root/'runtime/source.json').read_text())
    dest=Path(directory).expanduser().resolve()
    if dest.is_relative_to(root):
        raise ValueError('new external directory outside source checkouts required')
    dest.mkdir(parents=True,exist_ok=False)
    marker=dest/'salmonpy-runtime-source.json'
    try:
        store.copy_tree(root/'runtime',spec,dest)
        store.check_shared_kernels(root,dest,spec['shared_kernels'])
        subprocess.run([sys.executable,str(root/'tools/check_forbidden_eigensolvers.py'),'--source-root',str(dest/'src')],check=True)
        marker.write_text(json.dumps(spec,indent=2)+'\n')
        if spec.get('source_profile')=='fast-pr':
            fast=json.loads((root/'runtime/hse-rt-state.json').read_text())
            if (fast['source_tree'],fast['source_tree_sha256'])!=(spec['source_tree'],spec['source_tree_sha256']):
                raise ValueError('ordinary and fast HSE source manifests differ')
            fast.update(status='prepared-and-source-audited',routing_overlay=json.loads((root/'runtime/hse-explicit-routing.json').read_text()),
                        checker_sha256=hashlib.sha256((root/'tools/check_forbidden_eigensolvers.py').read_bytes()).hexdigest())
            (dest/'salmonpy-hse-rt-state.json').write_text(json.dumps(fast,indent=2)+'\n')
        return dest
    except BaseException:
        marker.write_text(json.dumps({**spec,'status':'failed-not-approved'},indent=2)+'\n')
        fast=json.loads((root/'runtime/hse-rt-state.json').read_text())
        fast['status']='failed-not-approved'
        (dest/'salmonpy-hse-rt-state.json').write_text(json.dumps(fast,indent=2)+'\n')
        raise

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory')
    args=parser.parse_args();print(prepare(args.directory))
