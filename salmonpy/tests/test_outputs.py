import tempfile
import unittest
from pathlib import Path
import numpy as np

class OutputTests(unittest.TestCase):
    def test_roundtrip_and_no_overwrite(self):
        from salmonpy.outputs import save_arrays, OutputArchive
        with tempfile.TemporaryDirectory() as tmp:
            dest=Path(tmp)/'result'
            save_arrays(dest,{'time':np.arange(3.),'rho':np.eye(2,dtype=complex)},metadata={'units':{'time':'a.u.','rho':'dimensionless'}})
            archive=OutputArchive.open(dest)
            np.testing.assert_array_equal(archive.array('rho'),np.eye(2))
            self.assertEqual(archive.metadata['units']['time'],'a.u.')
            with self.assertRaises(FileExistsError): save_arrays(dest,{'x':[1]},metadata={})
            with (dest/'arrays.npz').open('ab') as f:f.write(b'changed')
            with self.assertRaisesRegex(ValueError,'checksum'): OutputArchive.open(dest)

    def test_batch_archive(self):
        from salmonpy.crystal import respond
        from salmonpy.fields import FieldSamples
        from salmonpy.outputs import OutputArchive
        class Initial:
            backend='archive-fixture'
            def new_state(self,dt,time0,initial_A):
                class State:
                    time=time0
                    current=np.zeros(3)
                    def step(self,A):self.time+=dt;return A.copy()
                return State()
        with tempfile.TemporaryDirectory() as tmp:
            drive=FieldSamples([0.,.1],[[0,0,0],[1,2,3]],'A')
            dest=Path(tmp)/'run';result=respond(Initial(),drive,output_directory=dest)
            archive=OutputArchive.open(dest)
            np.testing.assert_array_equal(archive.array('A'),drive.values)
            np.testing.assert_array_equal(archive.array('J'),result.current)
            self.assertEqual(result.output_directory,str(dest.resolve()))

    def test_generated_native_files_only(self):
        from salmonpy.outputs import _write_archive
        from salmonpy.tddft import _digest
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'worker';source.mkdir()
            (source/'data_for_restart').mkdir()
            (source/'data_for_restart/wfn.bin').write_bytes(b'initial state')
            (source/'input').write_text('old')
            excluded={str(p.relative_to(source)):_digest(p) for p in source.rglob('*') if p.is_file()}
            (source/'input').write_text('updated')
            (source/'current.data').write_text('# 1:Time[a.u.] 2:J[a.u.]\n0 1\n')
            archive=_write_archive(root/'saved',{'time':[0.]},{},native_source=source,initial_hashes=excluded)
            self.assertFalse((archive.directory/'native/data_for_restart/wfn.bin').exists())
            self.assertEqual((archive.directory/'native/input').read_text(),'updated')
            self.assertEqual(archive.table('native/current.data').units,('a.u.','a.u.'))
            with self.assertRaises(ValueError):archive.table('../elsewhere')

    def test_cube_axes_and_scalar_order(self):
        from salmonpy.outputs import read_cube
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'density.cube'
            p.write_text('Electron Density\nAll values here are in a.u.\n1 1 2 3\n2 .5 0 0\n2 .1 .6 0\n2 0 0 .7\n14 14 0 0 0\n0 1 2 3 4 5 6 7\n')
            grid=read_cube(p)
            self.assertEqual(grid.values[1,0,1],5)
            np.testing.assert_array_equal(grid.axes[1],[.1,.6,0])
            np.testing.assert_array_equal(grid.origin,[1,2,3])
            self.assertEqual(grid.unit_note,'All values here are in a.u.')
            p.write_text(p.read_text()+'8\n')
            with self.assertRaises(ValueError):read_cube(p)

    def test_identical_preexisting_output_is_not_discarded(self):
        from salmonpy.outputs import _write_archive, _digest
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);source=root/'worker';source.mkdir()
            p=source/'Si_rt.data';p.write_text('# 1:Time[a.u.]\n0\n')
            archive=_write_archive(root/'saved',{'time':[0.]},{},native_source=source,initial_hashes={p.name:_digest(p)})
            self.assertEqual((archive.directory/'native'/p.name).read_bytes(),p.read_bytes())

    def test_failed_finish_retains_native_outputs_after_close(self):
        from unittest.mock import Mock,patch
        from types import SimpleNamespace
        from salmonpy.tddft import TDDFTState
        import shutil
        state=TDDFTState.__new__(TDDFTState)
        state._closed=False;state._dt=.1;state._steps=0;state._history=[];state._curr_history=[]
        state._A0=np.zeros(3);state._J0=np.zeros(3);state._initial_hashes={};state._threads=1;state._input='input'
        state._executable=SimpleNamespace(sha256='0'*64);state._reader=None
        state._owner=tempfile.TemporaryDirectory()
        from salmonpy.tddft import _WorkingDirectory
        state._work=_WorkingDirectory();state._retain_outputs=False
        work=Path(state._work.name);(work/'run').mkdir();(work/'run/current.data').write_text('valuable result')
        state._proc=Mock();state._proc.poll.return_value=0;state._proc.returncode=0;state._proc.stdin=None
        try:
            with tempfile.TemporaryDirectory() as tmp,patch('salmonpy.outputs._write_archive',side_effect=OSError('disk full')):
                with self.assertRaisesRegex(RuntimeError,'retained') as caught: state.finish(Path(tmp)/'output')
            state.close()
            self.assertEqual((work/'run/current.data').read_text(),'valuable result')
            self.assertEqual(Path(caught.exception.recovery_directory),work)
        finally:shutil.rmtree(work,ignore_errors=True)

    def test_decomposed_current_blocks(self):
        from salmonpy.outputs import read_decomposed_current
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'bands.data'
            p.write_text('# 1:ik[none] 2:ib[none] 3:is[none] 4:curr_decomp_x[a.u.] 5:curr_decomp_y[a.u.] 6:curr_decomp_z[a.u.]\n1\n1 1 1 1 2 3\n1 2 1 4 5 6\n8\n1 1 1 7 8 9\n1 2 1 10 11 12\n')
            data=read_decomposed_current(p)
            np.testing.assert_array_equal(data.step,[1,8])
            self.assertEqual(data.step.dtype.kind,'i')
            self.assertEqual(data.indices.dtype.kind,'i')
            np.testing.assert_array_equal(data.indices,[[1,1,1],[1,2,1]])
            np.testing.assert_array_equal(data.current.sum(axis=1),[[5,7,9],[17,19,21]])
            self.assertEqual(data.units,('a.u.',)*3)

    def test_interrupted_finish_retains_native_outputs(self):
        from unittest.mock import Mock,patch
        from types import SimpleNamespace
        from salmonpy.tddft import TDDFTState
        import shutil
        state=TDDFTState.__new__(TDDFTState)
        state._closed=False;state._dt=.1;state._steps=0;state._history=[];state._curr_history=[]
        state._A0=np.zeros(3);state._J0=np.zeros(3);state._initial_hashes={};state._threads=1;state._input='input'
        state._executable=SimpleNamespace(sha256='0'*64);state._reader=None
        state._owner=tempfile.TemporaryDirectory()
        from salmonpy.tddft import _WorkingDirectory
        state._work=_WorkingDirectory();state._retain_outputs=False
        work=Path(state._work.name);(work/'run').mkdir();(work/'run/current.data').write_text('valuable result')
        state._proc=Mock();state._proc.poll.return_value=0;state._proc.returncode=0;state._proc.stdin=None
        try:
            with tempfile.TemporaryDirectory() as tmp,patch('salmonpy.outputs._write_archive',side_effect=KeyboardInterrupt()):
                with self.assertRaises(KeyboardInterrupt) as caught: state.finish(Path(tmp)/'output')
            state.close()
            self.assertEqual((work/'run/current.data').read_text(),'valuable result')
            self.assertEqual(Path(caught.exception.recovery_directory),work)
        finally:shutil.rmtree(work,ignore_errors=True)

