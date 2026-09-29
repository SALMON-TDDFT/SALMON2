"""Build isolated dense-source, sparse-source and one-column routes; keep production unchanged."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import subprocess

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--build', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
root = Path(__file__).resolve().parents[2]
build = a.build.resolve()
out = a.output.resolve()
out.mkdir(parents=True, exist_ok=True)
source = (root / 'src/xc/exx_native.f90').read_text()
flag = 'logical,parameter :: use_sparse_source=.false.'
assert source.count(flag) == 1
flags = dict(re.findall(r'^(Fortran_\w+) = (.*)$', (build / 'src/CMakeFiles/salmon.dir/flags.make').read_text(), re.M))
link = shlex.split((build / 'src/CMakeFiles/salmon.dir/link.txt').read_text())
common = {item: hashlib.sha256((build / 'src' / item).read_bytes()).hexdigest()
          for item in link if item.endswith('.o') and item != 'CMakeFiles/salmon.dir/xc/exx_native.f90.o'}
for mode in ['before', 'after', 'boundary']:
    folder = out / mode
    folder.mkdir()
    text = source if mode == 'before' else source.replace(flag, flag.replace('.false.', '.true.'))
    if mode == 'boundary':
        text = text.replace('block_size=32', 'block_size=1')
    src = folder / 'exx_native.f90'
    src.write_text(text)
    compile_flags = [x for x in shlex.split(flags['Fortran_FLAGS']) if not x.startswith('-J')]
    compile_cmd = [link[0], *compile_flags, '-J'+str(folder), *shlex.split(flags['Fortran_INCLUDES']),
                   '-c', str(src), '-o', str(folder / 'exx_native.o')]
    subprocess.run(compile_cmd, cwd=folder, check=True)
    link_cmd = link.copy()
    link_cmd[link_cmd.index('CMakeFiles/salmon.dir/xc/exx_native.f90.o')] = str(folder / 'exx_native.o')
    link_cmd[link_cmd.index('-o')+1] = str(folder / 'salmon')
    subprocess.run(link_cmd, cwd=build / 'src', check=True)
    variant = dict(executable=str(folder / 'salmon'), compile=compile_cmd, link=link_cmd,
                   binary_sha256=hashlib.sha256((folder / 'salmon').read_bytes()).hexdigest(),
                   native_source_sha256=hashlib.sha256(src.read_bytes()).hexdigest())
    (folder / 'manifest.json').write_text(json.dumps(dict(kind='sparse-source routing '+mode,
        variants={'ace': variant}, shared_objects_sha256=common), indent=2)+'\n')
