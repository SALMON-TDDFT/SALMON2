#!/usr/bin/env python3
"""Compile-only ACE crash isolation; never modifies the production build.

Works with both the original and scalar-finiteness versions of exx_ace.f90.
Each probe has a private working directory. Do not link or run stub probes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import tempfile

MARKERS = {
    'finite_orbitals': ('finite=.false.', 'finite=.false.'),
    'finite_matrix': ('finite=.false.', 'finite=.false.'),
    'exx_ace_clear': ('ace=empty', 'return'),
    'exx_ace_ready': ('exx_ace_ready=allocated(', 'exx_ace_ready=.false.'),
    'exx_ace_average': ('ierr=1', 'ierr=1'),
    'exx_ace_build': ('ierr=1', 'ierr=1'),
    'exx_ace_apply': ('ierr=1', 'ierr=1'),
}


def variants(source):
    lines = source.splitlines(keepends=True)
    spans = {}
    current = None
    begin = None
    in_module_body = False
    for i, line in enumerate(lines):
        if line.strip().lower() == 'contains':
            in_module_body = True
            continue
        if not in_module_body:
            continue
        match = re.match(r'\s*(?:logical\s+)?(?:subroutine|function)\s+(\w+)\(', line, re.I)
        if match:
            current = match.group(1).lower()
            if current not in MARKERS:
                raise ValueError('Unrecognized procedure: ' + current)
            begin = None
        if current and line.strip().startswith(MARKERS[current][0]) and begin is None:
            begin = i
        if current and re.match(r'\s*end\s+(?:subroutine|function)\s*$', line, re.I):
            if begin is None:
                raise ValueError('Missing body marker: ' + current)
            spans[current] = (begin, i)
            current = None
    if current or not {'exx_ace_clear', 'exx_ace_ready', 'exx_ace_average',
                       'exx_ace_build', 'exx_ace_apply'} <= set(spans):
        raise ValueError('Source structure changed; refusing unsafe stub generation')

    def retain(names):
        result = list(lines)
        for name, (first, last) in sorted(spans.items(), key=lambda x: -x[1][0]):
            if name not in names:
                result[first:last] = ['    ' + MARKERS[name][1] + '\n']
        return ''.join(result)

    return [('all_stubs', retain(set()))] + [
        ('only_' + name, retain({name})) for name in spans
    ] + [('without_' + name, retain(set(spans) - {name})) for name in spans]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--build', type=Path, required=True)
    parser.add_argument('--source', type=Path)
    parser.add_argument('--compiler', default='mpifrtpx')
    parser.add_argument('--gnu-check', action='store_true')
    parser.add_argument('--timeout', type=int, default=120)
    a = parser.parse_args()
    build = a.build.resolve()
    source = (a.source or Path(__file__).resolve().parents[1]/'src/xc/exx_ace.f90').resolve()
    text = source.read_text()
    probes = variants(text)
    compiler = shutil.which(a.compiler)
    if not compiler:
        parser.error('Compiler not found: ' + a.compiler)
    flags_file = build/'src/CMakeFiles/salmon.dir/flags.make'
    flag_text = flags_file.read_text()
    flags = []
    for key in ('Fortran_DEFINES', 'Fortran_FLAGS'):
        match = re.search(r'^' + key + r'[ \t]*=[ \t]*(.*)$', flag_text, re.M)
        if not match:
            parser.error('Missing ' + key + ' in flags.make')
        flags += shlex.split(match.group(1))
    # Keep module output in the private case directory, including GNU local checks.
    safe_flags = []
    skip = False
    for flag in flags:
        if skip:
            skip = False
            continue
        if flag in ('-J', '-module'):
            skip = True
        elif not flag.startswith('-J'):
            safe_flags.append(flag)
    flags = safe_flags
    if any('$' in flag for flag in flags):
        parser.error('Unexpanded make variable in flags; inspect flags.make')
    out = Path(tempfile.mkdtemp(prefix='frtpx-ace-', dir=str(build.parent)))
    results = dict(source=str(source), source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                   scalar_finite_helpers='finite_orbitals' in text, compiler=compiler,
                   flags_make=flag_text, cases=[])
    version = subprocess.run([compiler, '--version' if a.gnu_check else '-V'],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True, timeout=30)
    results['version'] = version.stdout
    print('Output: ' + str(out), flush=True)
    print('Source SHA256: ' + results['source_sha256'], flush=True)
    print('Scalar finite helpers: ' + str(results['scalar_finite_helpers']), flush=True)
    print(version.stdout, flush=True)
    controls = [('full_build_flags', text, flags),
                ('full_O0', text, [f for f in flags if not re.match(r'^-O|^-Kfast', f)] + ['-O0']),
                ('full_no_openmp', text, [f for f in flags if f not in ('-Kopenmp', '-Nfjomplib', '-fopenmp')])]
    cases = controls + [(name, body, flags) for name, body in probes]
    for label, body, options in cases:
        folder = out/label
        folder.mkdir()
        probe = folder/'exx_ace.f90'
        probe.write_text(body)
        command = [compiler] + options + ['-c', str(probe), '-o', str(folder/'probe.o')]
        with (folder/'compile.log').open('w') as log:
            try:
                status = subprocess.run(command, cwd=folder, stdout=log,
                                        stderr=subprocess.STDOUT, timeout=a.timeout).returncode
            except subprocess.TimeoutExpired:
                status = 'timeout'
        log = (folder/'compile.log').read_text(errors='replace')
        row = dict(case=label, status=status, command=command,
                   sigsegv=status != 0 and ('SIGSEGV' in log or status == -11))
        results['cases'].append(row)
        (out/'summary.json').write_text(json.dumps(results, indent=2) + '\n')
        print('{}: status={} SIGSEGV={}'.format(label, status, row['sigsegv']), flush=True)
        if status != 0:
            print('\n'.join(log.splitlines()[-6:]), flush=True)
    print('Summary: ' + str(out/'summary.json'), flush=True)


if __name__ == '__main__':
    main()
