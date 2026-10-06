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
cases = [('threads/run.py', [], False), ('ace/test_conditioning.py', ['--build', '{build}'], False), ('ace/validate_exchange.py', ['--build', '{build}'], True), ('fft/run.py', ['--build', '{build}', '--action'], True), ('sr/run.py', ['--build', '{build}'], True), ('wannier/test_wannier.py', [], False)]
extras = [('k_exchange/run.py', [], False), ('cpu/run.py', [], False), ('wannier/test_gauge.py', [], False), ('wannier/test_locality.py', [], False), ('gauge/run.py', ['--build', '{build}'], False)]
for script, arguments, quick in cases + (extras if a.detailed else []):
    command = [sys.executable, str(here / script)]
    command += [str(a.build.resolve()) if value == '{build}' else value for value in arguments]
    if quick and not a.detailed:
        command.append('--quick')
    subprocess.run(command, check=True)
subprocess.run(['ctest', '--test-dir', str(a.build.resolve()), '-R', '^exx_history3(_coefficients)?$', '--output-on-failure', '--no-tests=error'], check=True)
print('Function checks passed:', here.name)
