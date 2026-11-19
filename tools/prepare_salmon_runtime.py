"""Prepare the bundled solver in a new external directory, without Git or patches.

One positional argument is the output. The historical UPSTREAM OUTPUT form is
accepted for compatibility, but the pinned repository-owned source is used.
"""
import argparse,json,subprocess,sys
from pathlib import Path

import importlib.util
_store_spec=importlib.util.spec_from_file_location('salmonpy_source_store',Path(__file__).with_name('runtime_sources.py'))
store=importlib.util.module_from_spec(_store_spec);_store_spec.loader.exec_module(store)

def prepare(source, directory=None):
    if directory is None:directory=source;source=None
    root=Path(__file__).resolve().parents[1]
    spec=json.loads((root/'runtime/source.json').read_text())
    dest=Path(directory).expanduser().resolve()
    if dest.is_relative_to(root) or (source is not None and dest.is_relative_to(Path(source).expanduser().resolve())):
        raise ValueError('destination must be outside source checkouts')
    dest.mkdir(parents=True,exist_ok=False)
    marker=dest/'salmonpy-runtime-source.json'
    try:
        store.copy_tree(root/'runtime',spec,dest)
        store.check_shared_kernels(root,dest,spec['shared_kernels'])
        subprocess.run([sys.executable,str(root/'tools/check_forbidden_eigensolvers.py'),'--source-root',str(dest/'src')],check=True)
        marker.write_text(json.dumps(spec,indent=2)+'\n')
        return dest
    except BaseException:
        marker.write_text(json.dumps({**spec,'status':'failed-not-approved'},indent=2)+'\n')
        raise

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source_or_output');parser.add_argument('directory',nargs='?')
    args=parser.parse_args();print(prepare(args.source_or_output,args.directory))
