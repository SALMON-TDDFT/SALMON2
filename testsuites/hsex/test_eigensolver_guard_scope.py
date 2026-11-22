"""Ensure generated application link guards resolve the mandatory checker."""
from pathlib import Path
import argparse

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--build',type=Path,required=True)
args=parser.parse_args()
files=list(args.build.rglob('build.make'))
if not files and (args.build/'build.ninja').exists():
    files=[args.build/'build.ninja']
guards=[line for path in files for line in path.read_text().splitlines() if '--object-tree' in line]
assert guards, 'No application-object eigensolver guards found in the generated build'
for line in guards:
    assert '"" --object-tree' not in line, 'Empty checker path in an application link guard'
    assert '/tools/check_forbidden_eigensolvers.py' in line, 'Unresolved application checker path'
print('APPLICATION_OBJECT_GUARD_SCOPE_PASS targets=',len(guards))
