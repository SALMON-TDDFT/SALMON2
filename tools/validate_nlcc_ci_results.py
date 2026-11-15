"""Check completed NLCC CI records; this is not runtime capability certification.

The native tests use prescribed density/core fields only, never retained GS data.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import re


def number(value):
    result = float(value.replace('D', 'E').replace('d', 'e'))
    if not math.isfinite(result):
        raise ValueError('nonfinite native diagnostic')
    return result


def load_build(path, expected_source, mpi):
    path = Path(path).resolve()
    record = json.loads(path.read_text())
    if record.get('status') != 'built-and-audited':
        raise ValueError('completed audited build required')
    if record.get('build_configuration', {}).get('mpi_enabled') is not mpi:
        raise ValueError('build MPI capability does not match role')
    for key in ('upstream_commit', 'patch_sha256') + (('source_tree_sha256',) if 'source_tree_sha256' in expected_source else ()):
        if record.get('source', {}).get(key) != expected_source[key]:
            raise ValueError('build source revision mismatch: '+key)
    exe = Path(record['executable']).resolve()
    if exe != path.parent/'build/salmon':
        raise ValueError('unexpected runtime executable location')
    if hashlib.sha256(exe.read_bytes()).hexdigest() != record.get('sha256'):
        raise ValueError('runtime executable checksum mismatch')
    return record, (path.parent/'source').resolve()


def validate_test(path, source, kind, expected_pairs, *, runtime_sha256, patch_sha256):
    record = json.loads(Path(path).read_text())
    if Path(record['source_root']).resolve() != Path(source).resolve():
        raise ValueError('native test used another source/build tree')
    if (record.get('runtime_executable_sha256') != runtime_sha256
            or record.get('runtime_source_patch_sha256') != patch_sha256):
        raise ValueError('native test runtime binding mismatch')
    expected_test = {'xc': 'test_nlcc_xc_variation', 'core': 'test_nlcc_core_force',
                     'mpi': 'test_nlcc_mpi_force'}[kind]
    if record.get('test') != expected_test or not re.fullmatch('[0-9a-f]{64}', record.get('executable_sha256', '')):
        raise ValueError('native test identity missing or mismatched')
    runs = record.get('runs', [])
    if any(type(r.get('ranks')) is not int or type(r.get('threads')) is not int for r in runs):
        raise ValueError('invalid native rank/thread types')
    pairs = [(r.get('ranks'), r.get('threads')) for r in runs]
    if len(pairs) != len(expected_pairs) or set(pairs) != set(expected_pairs):
        raise ValueError('missing or duplicate native rank/thread cases')
    maxima = []
    for run in runs:
        if type(run.get('returncode')) is not int or run['returncode'] != 0:
            raise ValueError('native test failed')
        stdout = run.get('stdout', '')
        if kind == 'xc':
            if 'NLCC_XC_VARIATION_PASS' not in stdout:
                raise ValueError('XC test success marker missing')
            cases = re.findall(r'^XC_VARIATION\s+(pz|libxc_pz|libxc_pbe)\s+([01])\s+(\S+)\s+(\S+)\s+(\S+)\s*$', stdout, re.M)
            if len(cases) != 6 or {(c[0], c[1]) for c in cases} != {(f, str(n)) for f in ('pz', 'libxc_pz', 'libxc_pbe') for n in (0, 1)}:
                raise ValueError('XC density/functionals coverage incomplete')
            for _, _, fd, inner, error in cases:
                fd, inner, error = map(number, (fd, inner, error))
                if error < 0 or error > 1e-9*max(1., abs(inner)) or abs(abs(fd-inner)-error) > 1e-12:
                    raise ValueError('XC variational check outside tolerance')
                maxima.append(error)
        elif kind == 'core':
            if 'NLCC_CORE_FORCE_PASS' not in stdout:
                raise ValueError('core-force test success marker missing')
            cases = re.findall(r'^\s*CORE_FORCE_VARIATION\s+([123])\s+(\S+)\s*$', stdout, re.M)
            if len(cases) != 3 or {c[0] for c in cases} != {'1', '2', '3'}:
                raise ValueError('core-force geometry coverage incomplete')
            for _, error in cases:
                error = number(error)
                if not 0 <= error <= 1e-7:
                    raise ValueError('core-force variation outside tolerance')
                maxima.append(error)
        else:
            cases = re.findall(r'^\s*MPI_NLCC_FORCE_PASS\s+(\d+)\s+(\S+)\s+(\d+)\s*$', stdout, re.M)
            if len(cases) != 1:
                raise ValueError('MPI test success marker missing or duplicated')
            ranks, error, rejected = cases[0]
            error = number(error)
            if int(ranks) != run['ranks'] or int(rejected) != run['ranks'] or not 0 <= error <= 1e-13:
                raise ValueError('MPI replication or collective rejection failed')
            maxima.append(error)
    return {'test': expected_test, 'test_executable_sha256': record['executable_sha256'],
            'cases': len(runs), 'rank_thread_pairs': [list(pair) for pair in sorted(pairs)],
            'record_sha256': hashlib.sha256(Path(path).read_bytes()).hexdigest(),
            'maximum_reported_error': max(maxima)}


def validate_native_records(serial_build_record, xc_record, core_record, *, source_spec,
                            mpi_build_record=None, mpi_record=None):
    """Validate required serial records and optional MPI records, fail closed.

    Inputs are record paths and the current parsed runtime/source.json. This
    normalizes native evidence only; whole-app/runtime certification is separate.
    """
    serial_build, serial_source = load_build(serial_build_record, source_spec, False)
    bindings = dict(runtime_sha256=serial_build['sha256'], patch_sha256=source_spec['patch_sha256'])
    tests = {'xc': validate_test(xc_record, serial_source, 'xc', [(1, 1), (1, 2)], **bindings),
             'core': validate_test(core_record, serial_source, 'core', [(1, 1), (1, 2)], **bindings)}
    result = {'schema': 'salmonpy.nlcc-native-validation/1', 'status': 'passed',
              'scope': 'native consistency only; not physical convergence or runtime capability certification',
              'build_executable_sha256': serial_build['sha256'],
              'source_patch_sha256': source_spec['patch_sha256'], 'tests': tests}
    if (mpi_build_record is None) != (mpi_record is None):
        raise ValueError('MPI build and test records must be provided together')
    if mpi_build_record is not None:
        mpi_build, mpi_source = load_build(mpi_build_record, source_spec, True)
        tests['mpi'] = validate_test(mpi_record, mpi_source, 'mpi',
            [(2, 1), (2, 2), (4, 1), (4, 2)], runtime_sha256=mpi_build['sha256'],
            patch_sha256=source_spec['patch_sha256'])
        tests['mpi']['build_executable_sha256'] = mpi_build['sha256']
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('serial_build', 'xc_test', 'core_test'):
        parser.add_argument(name, type=Path)
    parser.add_argument('--mpi-build-record', type=Path)
    parser.add_argument('--mpi-record', type=Path)
    parser.add_argument('--source-spec', type=Path, default=Path(__file__).resolve().parents[1]/'runtime/source.json')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = validate_native_records(args.serial_build, args.xc_test, args.core_test,
        mpi_build_record=args.mpi_build_record, mpi_record=args.mpi_record,
        source_spec=json.loads(args.source_spec.read_text()))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    tmp = args.output.with_name(args.output.name+'.tmp')
    tmp.write_text(json.dumps(result, indent=2)+'\n')
    tmp.replace(args.output)
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
