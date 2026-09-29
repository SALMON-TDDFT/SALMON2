"""Link identical build objects with native EXX routing variants for a controlled benchmark.

No production source or build object is modified. 'replicated' uses the existing
replicated paths in the current source; it is not a historical executable.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--build', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
b = a.build.resolve()
out = a.output.resolve()
root = Path(__file__).resolve().parents[2]
source = (root / 'src/xc/exx_native.f90').read_text()
# Production uses distributed ACE with replicated MLWF. Construct the optional
# all-distributed comparison route only in these isolated benchmark copies.
old = '    if(info%isize_o>1)orbital_comm=info%icomm_o'
assert source.count(old) == 2
source = source.replace(old, '    orbital_comm=info%icomm_o')
old = 'maxiter,exx_mlwf_tolerance,status,occupation=system%rocc(info%io_s:info%io_e,:,1),comm_o=orbital_comm)'
assert source.count(old) == 1
source = source.replace(old, old[:-1] + ', &\n      comm_matrix=info%icomm_ro)')
flags = dict(re.findall(r'^(Fortran_\w+) = (.*)$', (b / 'src/CMakeFiles/salmon.dir/flags.make').read_text(), re.M))
link = shlex.split((b / 'src/CMakeFiles/salmon.dir/link.txt').read_text())
compiler = link[0]
out.mkdir(parents=True, exist_ok=True)
manifest = dict(kind='current-source routing comparison', build=str(b), variants={},
                source_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=root, text=True).strip())
(out / 'working-tree.patch').write_bytes(subprocess.check_output(['git', 'diff'], cwd=root))
# Archive new, as-yet-untracked production modules as well as their hashes.
for name in ['exx_distributed_metric.f90', 'exx_distributed_gauge.f90']:
    shutil.copy2(root / 'src/xc' / name, out / name)
for mode in ['replicated', 'ace', 'all']:
    folder = out / mode
    folder.mkdir()
    text = source
    if mode != 'all':
        text = text.replace('    orbital_comm=info%icomm_o', '    if(info%isize_o>1)orbital_comm=info%icomm_o')
        old = 'comm_o=orbital_comm, &\n      comm_matrix=info%icomm_ro)'
        assert text.count(old) == 1
        text = text.replace(old, 'comm_o=orbital_comm)')
    if mode == 'replicated':
        old = 'status,packed=.true., &\n        comm_matrix=info%icomm_ro)'
        assert text.count(old) == 1
        text = text.replace(old, 'status,packed=.true.)')
        old = 'call orbital_ace_build(ace,local,w,system%hvol,info%icomm_r,info%icomm_o,status,comm_matrix=info%icomm_ro)'
        assert text.count(old) == 1
        text = text.replace(old, 'call orbital_ace_build(ace,local,w,system%hvol,info%icomm_r,info%icomm_o,status)')
        old = 'subroutine build_exchange_ace()\n      implicit none\n      if(info%isize_o>1.or.info%isize_r>1)then'
        assert text.count(old) == 1
        text = text.replace(old, old.replace('.or.info%isize_r>1', ''))
    f = folder / 'exx_native.f90'
    f.write_text(text)
    compile_flags = shlex.split(flags['Fortran_FLAGS'])
    compile_flags = [v for v in compile_flags if not v.startswith('-J')]
    cmd = [compiler, *compile_flags, '-J'+str(folder), *shlex.split(flags['Fortran_INCLUDES']),
           '-c', str(f), '-o', str(folder / 'exx_native.o')]
    subprocess.run(cmd, cwd=folder, check=True)
    cmd_link = link.copy()
    index = cmd_link.index('CMakeFiles/salmon.dir/xc/exx_native.f90.o')
    cmd_link[index] = str(folder / 'exx_native.o')
    cmd_link[cmd_link.index('-o')+1] = str(folder / 'salmon')
    subprocess.run(cmd_link, cwd=b / 'src', check=True)
    manifest['variants'][mode] = dict(executable=str(folder / 'salmon'),
        binary_sha256=hashlib.sha256((folder / 'salmon').read_bytes()).hexdigest(),
        native_source_sha256=hashlib.sha256(f.read_bytes()).hexdigest(), compile=cmd, link=cmd_link)
manifest['shared_objects_sha256'] = {
    item: hashlib.sha256((b / 'src' / item).read_bytes()).hexdigest()
    for item in link if item.endswith('.o') and item != 'CMakeFiles/salmon.dir/xc/exx_native.f90.o'}
(out / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
print('Built replicated / ACE-only / ACE+MLWF variants:', out)
