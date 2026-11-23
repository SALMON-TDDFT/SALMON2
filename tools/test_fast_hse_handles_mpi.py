"""Exercise real generation-checked C contexts: BUILD FULL_K_FIXTURE NEW_OUTPUT.

Launch with 2*--group-size host MPI ranks. Explicitly built current native library
required; no compilation/download. Compare mixed output with separate single-a/b
runs. Synthetic state ownership/restart verification, not material convergence.
"""
import argparse
import ctypes as c
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import numpy as np
from mpi4py import MPI


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('build', 'fixture', 'output'):
        parser.add_argument(key, type=Path)
    parser.add_argument('--group-size', type=int, default=2)
    parser.add_argument('--mode', choices=['single-a', 'single-b', 'mixed'], default='mixed')
    parser.add_argument('--legacy', action='store_true', help='isolated old-singleton comparison only')
    args = parser.parse_args()
    if args.legacy and args.mode == 'mixed':
        parser.error('--legacy only supports isolated controls')
    world = MPI.COMM_WORLD
    assert world.size == 2 * args.group_size
    group = world.rank // args.group_size
    comm = world.Split(group, world.rank)
    build = args.build if args.build.is_file() else args.build / 'build.json'
    receipt = json.loads(build.read_text())
    library = Path(receipt['library'])
    assert hashlib.sha256(library.read_bytes()).hexdigest() == receipt['sha256']
    lib = c.CDLL(str(library))
    vec = np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, shape=(3,), flags='C_CONTIGUOUS')
    ints = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, shape=(8,), flags='C_CONTIGUOUS')
    token_type = c.c_int64
    specs = {'protocol': [], 'init': [c.c_char_p, c.c_double, c.c_int, c.c_int, c.POINTER(token_type), vec],
             'new': [token_type, c.c_double, c.c_int, c.POINTER(token_type), vec],
             'step': [token_type, vec, vec, vec, c.POINTER(c.c_double)],
             'info': [token_type, c.POINTER(c.c_int), vec, vec],
             'history': [token_type, ints], 'transaction': [token_type, c.c_int],
             'checkpoint': [token_type, c.c_char_p, c.c_int],
             'label': [token_type, c.c_void_p, c.c_int], 'close': [token_type]}
    api = {}
    if not args.legacy:
        for name, spec in specs.items():
            fn = getattr(lib, 'sp_fast_hse_handle_' + name)
            fn.argtypes = spec
            fn.restype = c.c_int
            api[name] = fn
        assert api['protocol']() == 1
    else:
        # Explicit reference path, not a fallback for the handle implementation.
        old_specs = {'init': [c.c_char_p, c.c_double, c.c_int, vec],
                     'step': [vec, vec, vec, c.POINTER(c.c_double)],
                     'info': [c.POINTER(c.c_int), vec, vec],
                     'history_info': [ints], 'checkpoint': [c.c_char_p, c.c_int], 'close': []}
        old = {}
        for name, spec in old_specs.items():
            fn = getattr(lib, 'sp_fast_hse_' + name)
            fn.argtypes, fn.restype = spec, c.c_int
            old[name] = fn
        def old_init(path, dt, capacity, symmetry, handle, j):
            code = old['init'](path, dt, capacity, j)
            handle._obj.value = 1 if code == 0 else 0
            return code
        def old_label(token, buffer, capacity):
            c.memmove(buffer, b'native_fast_hse\0', len(b'native_fast_hse\0'))
            return 0
        api = {'init': old_init, 'label': old_label,
               'step': lambda token, *args: old['step'](*args),
               'info': lambda token, *args: old['info'](*args),
               'history': lambda token, *args: old['history_info'](*args),
               'checkpoint': lambda token, *args: old['checkpoint'](*args),
               'close': lambda token: old['close']()}
    lib.sp_attach.argtypes = [c.c_int64]
    lib.sp_attach.restype = c.c_int
    lib.sp_detach.argtypes = []
    lib.sp_detach.restype = c.c_int
    assert lib.sp_attach(comm.py2f()) == 0
    work = args.output.resolve() / f'group-{group}'
    if comm.rank == 0:
        shutil.copytree(args.fixture, work)
        text = (work / 'rt.inp').read_text()
        text = re.sub(r'\bnproc_k\s*=\s*\d+', 'nproc_k=' + str(comm.size), text)
        if group:
            text = text.replace("xc='hse06'", "xc='hse06',hse_omega=.12")
        (work / 'rt.inp').write_text(text)
    comm.Barrier()
    os.chdir(work)
    labels = {}
    handles = {}
    fields = {}
    rows = {'a': [], 'b': []}
    checks = []

    def label(handle):
        buffer = c.create_string_buffer(81)
        assert api['label'](handle, buffer, len(buffer)) == 0
        return buffer.value.decode()

    def opened(owner, model=None):
        dt, capacity = (.002, 34) if owner in ('b', 'copy-b') else (.001, 32)
        handle, j = token_type(), np.zeros(3)
        if model is None:
            code = api['init'](b'rt.inp', dt, capacity, 0, c.byref(handle), j)
        else:
            code = api['new'](model, dt, capacity, c.byref(handle), j)
        assert code == 0, ('open', owner, code)
        handles[owner] = handle.value
        fields[owner] = np.zeros(3)
        labels[owner] = label(handle.value)
        return handle.value

    def info(handle):
        index, a, j = c.c_int(), np.zeros(3), np.zeros(3)
        code = api['info'](handle, c.byref(index), a, j)
        return code, index.value, a, j

    def checkpoint(owner, tag, mode=0):
        path = work / f'{tag}-{comm.rank}.bin'
        code = api['checkpoint'](handles[owner], os.fsencode(path) if mode < 3 else b'', mode)
        assert code == 0, ('checkpoint', owner, mode, code)
        return path

    def advance(owner, step, append=True):
        a = np.zeros(3)
        if owner in ('b', 'copy-b'):
            a[:2] = .0001 * np.array([np.cos(.2 * step) - 1, np.sin(.2 * step)])
        else:
            a[2] = .0001 * np.sin(.3 * step)
        j, energy = np.zeros(3), c.c_double()
        code = api['step'](handles[owner], fields[owner], a, j, c.byref(energy))
        assert code == 0, ('step', owner, step, code)
        fields[owner] = a
        code, index, value, current = info(handles[owner])
        assert code == 0 and index == step
        np.testing.assert_array_equal(value, a)
        np.testing.assert_array_equal(current, j)
        row = np.r_[j, energy.value]
        if append:
            rows[owner].append(row)
        return row

    def files(owner):
        return {p.name: p.read_bytes() for p in work.glob(labels[owner] + '*.data')}

    def begin(owner):
        assert api['transaction'](handles[owner], 0) == 0
        return files(owner)

    def rollback(owner, saved):
        assert api['transaction'](handles[owner], 1) == 0
        if comm.rank == 0:
            for name, data in saved.items():
                (work / name).write_bytes(data)
        comm.Barrier()
        assert api['transaction'](handles[owner], 2) == 0
        fields[owner] = info(handles[owner])[2]

    if comm.size > 1 and not args.legacy:
        bad_handle, j = token_type(), np.zeros(3)
        assert api['init'](b'rt.inp', .001 + .001 * comm.rank, 32, 0,
                           c.byref(bad_handle), j) == 11
        assert bad_handle.value == 0
        assert api['init'](b'rt.inp', .001, 32 + comm.rank, 0,
                           c.byref(bad_handle), j) == 11
        assert bad_handle.value == 0
        checks.append('mixed-init-rejected-before-lifetime')
    mixed = args.mode == 'mixed'
    if mixed or args.mode == 'single-a':
        opened('a')
    if mixed or args.mode == 'single-b':
        opened('b', handles.get('a'))
    if mixed:
        assert handles['a'] != handles['b'] and labels['a'] != labels['b']
        # Same function, two distinct valid tokens: reject collectively before kernels.
        bad = handles['a'] if comm.rank == 0 else handles['b']
        if comm.size > 1:
            assert info(bad)[0] == 11
        for bad in (0, -1, 2**63 - 1):
            assert info(bad)[0] == 12
        tiny = c.create_string_buffer(1)
        assert api['label'](handles['b'], tiny, 1) == 13
        assert label(handles['b']) == labels['b']
        checks += ['invalid-and-mixed-rank-tokens']
    for step in range(1, 33):
        if mixed and step == 6:
            saved_a, saved_b = begin('a'), begin('b')
        if mixed and step == 13:
            saved_a = begin('a')
            assert api['transaction'](handles['a'], 2) == 0
            before = checkpoint('a', 'before-overflow').read_bytes()
            j, energy = np.zeros(3), c.c_double()
            bad = np.array([0., 0., 1e200])
            code = api['step'](handles['a'], fields['a'], bad, j, c.byref(energy))
            assert code in (*range(20, 26), 30, 32, 33, 34), ('overflow', code)
            assert checkpoint('a', 'after-overflow').read_bytes() == before
            assert files('a') == saved_a
            checks.append('numeric-failure-retains-selected-context')
        for owner in tuple(handles):
            advance(owner, step)
        if mixed and step == 6:
            b_before = info(handles['b'])
            rollback('a', saved_a)
            np.testing.assert_array_equal(advance('a', step, False), rows['a'][-1])
            after = info(handles['b'])
            assert b_before[:2] == after[:2]
            np.testing.assert_array_equal(b_before[2:], after[2:])
            rollback('b', saved_b)
            np.testing.assert_array_equal(advance('b', step, False), rows['b'][-1])
            checks.append('simultaneous-transactions')
        if mixed and step == 9:
            checkpoint('a', 'a-nine');checkpoint('b', 'b-nine')
        if mixed and step == 10:
            checkpoint('a', 'a-nine', 2);checkpoint('b', 'b-nine', 2)
            checkpoint('a', '', 3)
            fields['a'] = info(handles['a'])[2]
            np.testing.assert_array_equal(advance('a', step, False), rows['a'][-1])
            checkpoint('b', '', 4)
            checks.append('independent-staged-commit-cancel')
    if mixed:
        checkpoint('a', 'a-final')
        history = np.zeros(8, dtype=np.int32)
        assert api['history'](handles['a'], history) == 0
        np.save(work / f'a-history-{comm.rank}.npy', history)
        begin('a')
        checkpoint('b', 'b-thirtytwo');checkpoint('b', 'b-thirtytwo', 2)
        old = handles.pop('a')
        assert api['close'](old) == 0
        assert info(old)[0] == 12
        checkpoint('b', '', 3)
        fields['b'] = info(handles['b'])[2]
        checks.append('root-close-preserves-child-candidate')
    if 'b' in handles:
        advance('b', 33)
    if mixed:
        fresh = opened('c', handles['b'])
        assert fresh != old and fresh % 2**31 == old % 2**31
        assert info(old)[0] == 12 and api['close'](old) == 12
        for step in range(1, 4):
            np.testing.assert_array_equal(advance('c', step, False), rows['a'][step - 1])
        restored = opened('d', handles['b'])
        assert api['checkpoint'](restored, os.fsencode(work / f'a-nine-{comm.rank}.bin'), 1) == 0
        code, index, a, j = info(restored)
        assert code == 0 and index == 9
        fields['d'] = a
        for resumed in range(10, 33):
            np.testing.assert_array_equal(advance('d', resumed, False), rows['a'][resumed - 1])
        assert api['close'](handles.pop('d')) == 0
        checks.append('fresh-direct-load-without-replay')
        restored_b = opened('copy-b', handles['b'])
        assert api['checkpoint'](restored_b, os.fsencode(work / f'b-nine-{comm.rank}.bin'), 1) == 0
        fields['copy-b'] = info(restored_b)[2]
        for resumed in range(10, 34):
            np.testing.assert_array_equal(advance('copy-b', resumed, False), rows['b'][resumed - 1])
        assert api['close'](handles.pop('copy-b')) == 0
        begin('c')
        checkpoint('b', 'b-final-stage');checkpoint('b', 'b-final-stage', 2)
        # Grow the registry while numerical buffers and a restore candidate are live.
        temporary = []
        for _ in range(6):
            handle, j = token_type(), np.zeros(3)
            assert api['new'](handles['b'], .001, 32, c.byref(handle), j) == 0
            temporary.append(handle.value)
        assert len(set(temporary + list(handles.values()))) == 8
        for handle in temporary:
            assert api['close'](handle) == 0
            assert info(handle)[0] == 12
        assert api['transaction'](fresh, 2) == 0
        checkpoint('b', '', 3)
        assert api['close'](handles.pop('c')) == 0
        checks += ['generation-reuse', 'registry-growth-with-live-candidates']
    if 'b' in handles:
        advance('b', 34)
    for owner in list(handles):
        checkpoint(owner, owner + '-final')
        history = np.zeros(8, dtype=np.int32)
        assert api['history'](handles[owner], history) == 0
        np.save(work / f'{owner}-history-{comm.rank}.npy', history)
        assert api['close'](handles[owner]) == 0
        assert info(handles[owner])[0] != 0
    for owner, trace in rows.items():
        if trace:
            np.save(work / f'{owner}-trace-{comm.rank}.npy', trace)
    assert lib.sp_detach() == 0
    assert MPI.Is_initialized() and not MPI.Is_finalized()
    if comm.rank == 0:
        (work / 'validation.json').write_text(json.dumps({'status': 'pass', 'mode': args.mode,
            'group_size': comm.size, 'legacy': args.legacy, 'checks': checks, 'labels': labels,
            'library_sha256': receipt['sha256'], 'replayed_steps': 0}, indent=2) + '\n')
        print('Generation-checked C contexts PASS', args.mode, group, flush=True)
    comm.Free()


if __name__ == '__main__':
    main()
