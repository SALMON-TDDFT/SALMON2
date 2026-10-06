"""One functional regression entry point; --detailed expands the numerical matrix."""
import argparse
from pathlib import Path
import subprocess
import sys
p = argparse.ArgumentParser()
p.add_argument('--build', type=Path, required=True)
p.add_argument('--detailed', action='store_true')
a = p.parse_args()
here = Path(__file__).resolve().parent
cases = [('finite/run.py', [], False), ('mesh/test_tile.py', [], False), ('mesh/test_stream.py', [], False), ('scalapack/validate.py', ['--build', '{build}'], True), ('rt/check.py', [], False), ('rt/test_transport.py', [], False), ('rt/test_occupation_guard.py', [], False)]
extras = []
for script, arguments, quick in cases + (extras if a.detailed else []):
    command = [sys.executable, str(here / script)]
    command += [str(a.build.resolve()) if value == '{build}' else value for value in arguments]
    if quick and not a.detailed:
        command.append('--quick')
    subprocess.run(command, check=True)
print('Function checks passed:', here.name)
