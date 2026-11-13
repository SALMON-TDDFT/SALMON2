"""Read-only setup diagnosis: python -m salmonpy.doctor [--json].

Loads the bundled kernels to check runtime dependencies and symbols, but never
runs a calculation, compiles code, installs packages or launches an executable.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import sys


def diagnose(*, executable=None, sha256=None):
    from . import ld, sbe
    native = {}
    for name, module in (('ld', ld), ('sbe', sbe)):
        try:
            module._native()
            native[name] = {'status': 'loadable', 'detail': '数値精度・計算規模は未評価'}
        except (OSError, RuntimeError, AttributeError) as error:
            native[name] = {'status': 'unavailable', 'detail': str(error)}
    optional = {name: 'discoverable' if importlib.util.find_spec(name) else 'missing'
                for name in ('matplotlib', 'ase', 'f90nml')}
    external = {'status': 'not_configured', 'capabilities': 'unverified'}
    if executable is not None:
        path = Path(executable).expanduser().resolve()
        external['path'] = str(path)
        if not path.is_file():
            external['status'] = 'missing'
        elif sha256 is None:
            external['status'] = 'hash_required'
        elif not isinstance(sha256, str) or not re.fullmatch('[0-9a-f]{64}', sha256):
            external['status'] = 'invalid_hash'
        else:
            try:
                digest = hashlib.sha256()
                with path.open('rb') as stream:
                    for block in iter(lambda: stream.read(1024*1024), b''):
                        digest.update(block)
                external['status'] = 'hash_matches' if digest.hexdigest() == sha256 else 'hash_mismatch'
            except OSError as error:
                external.update(status='unreadable', detail=str(error))
    next_steps = ['目的・使ってよい総メモリ・CPUコア数を確認してください。検出したCPU数は利用許可ではありません。']
    if any(v['status'] != 'loadable' for v in native.values()):
        next_steps.append('LD/SBE: 対応OS・Pythonの修復済みwheelを再導入してください。ソース構築はPACKAGING.mdを参照。')
    else:
        next_steps.append('小さな実行例: python -m salmonpy.examples.ld_response')
    if optional['matplotlib'] == 'missing':
        next_steps.append('描画を使う場合のみMatplotlibを導入してください。')
    if any(optional[n] == 'missing' for n in ('ase', 'f90nml')):
        next_steps.append('ASE連携を使う場合はASE.mdの追加依存を導入してください。')
    next_steps.append('DFT/TDDFT/バンド/FDTD: 対応する外部SALMONと監査済みSHA256、入力データを準備してください。ハッシュ一致だけでは機能対応・実行可能性を判定できません。')
    return {'schema_version': 'salmonpy.doctor/1',
            'environment': {'python': platform.python_version(), 'platform': platform.platform(),
                            'interpreter': sys.executable},
            'native': native, 'optional_packages': optional, 'external_salmon': external,
            'resources': {'detected_logical_cpus': os.cpu_count(), 'memory_budget': 'not_specified',
                          'cpu_budget': 'not_specified', 'maximum_problem_size': 'unknown',
                          'hard_memory_enforcement': False},
            'limitations': ['Optional packages are discovered, not imported or tested.',
                            'Native kernels are loaded; numerical calculations are not run.',
                            'External executable capabilities and GS compatibility are not tested.'],
            'next_steps': next_steps}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json', action='store_true', help='AI向けJSONを表示')
    parser.add_argument('--executable', help='外部SALMONのパス（起動しません）')
    parser.add_argument('--sha256', help='監査済み実行ファイルのSHA256')
    args = parser.parse_args()
    if args.sha256 and not args.executable:
        parser.error('--sha256 requires --executable')
    report = diagnose(executable=args.executable, sha256=args.sha256)
    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        print('SALMONpy 環境診断（計算・コンパイルは行いません）')
        print('Python:', report['environment']['python'])
        for name, value in report['native'].items():
            print(f'{name.upper()}: {value["status"]} — {value["detail"]}')
        print('追加依存:', ', '.join(f'{k}={v}' for k, v in report['optional_packages'].items()))
        print('外部SALMON:', report['external_salmon']['status'], '（機能対応は未検証）')
        print('検出CPU数:', report['resources']['detected_logical_cpus'], '／利用上限と計算規模は未設定・不明')
        for step in report['next_steps']:
            print('- ' + step)


if __name__ == '__main__':
    main()
