import unittest
from pathlib import Path
import tempfile

class TDDFTLoaderTests(unittest.TestCase):
    def test_bundle_is_snapshotted_and_validated(self):
        from salmonpy.tddft import TDDFTInitialState
        from salmonpy.fdtd import SalmonExecutable
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); restart=root/'data_for_restart';restart.mkdir()
            for name in ('info.bin','occupation.bin','wfn.bin'):
                (restart/name).write_bytes(b'synthetic loader test, not physics')
            (root/'rt.inp').write_text('&calculation\n theory="tddft_pulse"\n/\n')
            exe=SalmonExecutable(root/'salmon','0'*64)
            initial=TDDFTInitialState.from_salmon(root,'rt.inp',executable=exe,max_steps=10)
            original=initial.manifest
            (restart/'wfn.bin').write_bytes(b'changed outside snapshot')
            self.assertEqual(initial.manifest,original)
            self.assertEqual(initial.backend,'SALMON-TDDFT-middlepoint')
            with self.assertRaises(AttributeError): initial.max_steps=12
            with self.assertRaises(ValueError): initial.new_state(dt=-1)
            with self.assertRaises(ValueError): initial.new_state(dt=.1,time0=1.)
            initial.close()
            with self.assertRaises(RuntimeError): initial.new_state(dt=.1)

    def test_reject_incomplete_bundle(self):
        from salmonpy.tddft import TDDFTInitialState
        from salmonpy.fdtd import SalmonExecutable
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);(root/'input').write_text('test')
            with self.assertRaises(ValueError):
                TDDFTInitialState.from_salmon(root,'input',executable=SalmonExecutable(root/'s','0'*64))


    def test_snapshot_cannot_recursively_copy_itself(self):
        from salmonpy.tddft import _copy_bundle
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaisesRegex(ValueError,'outside'):
                _copy_bundle(root,root/'nested')
            self.assertFalse((root/'nested').exists())

    def test_send_timeout_kills_nonreading_worker(self):
        import subprocess, sys, time
        from salmonpy.tddft import TDDFTState
        state = object.__new__(TDDFTState)
        state._closed = True
        state._timeout = .1
        state._proc = subprocess.Popen([sys.executable, '-c', 'import time; time.sleep(10)'],
            stdin=subprocess.PIPE, stdout=subprocess.DEVNULL, text=True)
        start = time.monotonic()
        try:
            with self.assertRaises(TimeoutError):
                state._send('x' * 1000000)
            self.assertLess(time.monotonic()-start, 3.)
            self.assertIsNotNone(state._proc.poll())
        finally:
            if state._proc.poll() is None: state._proc.kill()
            state._proc.wait()
            state._proc.stdin.close()
