from pathlib import Path
import tempfile
import unittest
import numpy as np
from salmonpy.dft import DFT, Pseudopotential, SCFSettings, _read_status, GroundStateResult

class DFTTest(unittest.TestCase):
    def test_model_copies_and_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'Si_rps.dat';path.write_bytes(b'original')
            pseudo=Pseudopotential.from_file(path,atomic_number=14,local_channel=2)
            cell=np.eye(3)*10.26;positions=np.array([[0.,0.,0.]])
            model=DFT([14],cell,positions,[pseudo],electrons=4,states=8,grid=[12]*3,kmesh=[2]*3)
            cell[0,0]=99;positions[0,0]=.5;path.write_bytes(b'changed')
            self.assertEqual(model.cell[0,0],10.26)
            self.assertEqual(model.fractional_positions[0,0],0)
            self.assertEqual(pseudo.data,b'original')
            text=model.input_text(SCFSettings(max_iterations=100,tolerance=1e-8))
            self.assertIn("theory='dft'",text)
            self.assertIn('nscf=100',text)
            self.assertIn('_rps.dat',text)
            from dataclasses import replace
            hybrid=replace(model,xc='hse06')
            self.assertIn("xc='hse06'",hybrid.input_text())
            from salmonpy.fdtd import SalmonExecutable
            import hashlib
            binary=Path(tmp)/'solver';binary.write_bytes(b'not executable')
            engine=SalmonExecutable(binary,hashlib.sha256(binary.read_bytes()).hexdigest())
            with self.assertRaisesRegex(ValueError,'NativeGSBackend'):
                hybrid.run(engine,Path(tmp)/'must-not-run')
            self.assertFalse((Path(tmp)/'must-not-run').exists())
            with self.assertRaises(ValueError):model.cell.setflags(write=True)
            with self.assertRaises(ValueError):DFT([14],np.eye(3),[[0,0,0]],[pseudo],electrons=3,states=8,grid=[12]*3,kmesh=[2]*3)

    def test_native_status_and_incomplete_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'gs_status.dat'
            p.write_text('SALMONPY_GS_V1\n1 37\n1e-10 1e-9 -31.25\nrho_dne\n')
            self.assertEqual(_read_status(p)['iterations'],37)
            p.write_text(p.read_text().replace('1 37','0 37'))
            with self.assertWarnsRegex(RuntimeWarning,'SCF.*not converged'):
                self.assertFalse(_read_status(p)['converged'])
            for values in ('nan 1e-9 -31.25', '-1 1e-9 -31.25', '1e-10 0 -31.25'):
                p.write_text('SALMONPY_GS_V1\n0 37\n'+values+'\nrho_dne\n')
                with self.assertRaises(RuntimeError):_read_status(p)
            with self.assertRaises((ValueError,FileNotFoundError)):GroundStateResult.open(tmp)

    def test_extra_native_files_rejected_on_reopen(self):
        import json,hashlib
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);native=root/'native';(native/'data_for_restart').mkdir(parents=True)
            for name in ['info.bin','occupation.bin','wfn.bin']:(native/'data_for_restart'/name).write_bytes(b'x')
            (native/'gs.inp').write_text('input')
            (native/'gs_status.dat').write_text('SALMONPY_GS_V1\n1 37\n1e-10 1e-9 -31.25\nrho_dne\n')
            files={p.relative_to(root).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in native.rglob('*') if p.is_file()}
            (root/'ground-state.json').write_text(json.dumps({'schema':'salmonpy.gs/1','status':'complete','energy_hartree':-31.25,'iterations':37,'files':files}))
            GroundStateResult.open(root)
            (native/'hubbard_param.dat').write_text('unexpected physical input')
            with self.assertRaisesRegex(ValueError,'files'):GroundStateResult.open(root)

    def test_forces_shape_units_and_readonly(self):
        from salmonpy.dft import _read_forces
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'forces'
            p.write_text('SALMONPY_FORCES_V1\n2\n.1 -.2 0\n-.1 .2 0\n')
            f=_read_forces(p,2)
            np.testing.assert_array_equal(f,[[.1,-.2,0],[-.1,.2,0]])
            with self.assertRaises(ValueError):f.setflags(write=True)
            with self.assertRaises(ValueError):_read_forces(p,3)
            p.write_text('SALMONPY_FORCES_V1\n1\nnan 0 0\n')
            with self.assertRaises(ValueError):_read_forces(p,1)

class ValenceTest(unittest.TestCase):
    def test_ky_and_upf_charge(self):
        ky = Pseudopotential(14, 2, '_rps.dat', b'2000 .01 2 4.0D+0\n')
        self.assertEqual(ky.valence_electrons, 4)
        upf = Pseudopotential(14, 2, '.upf', b'<UPF><PP_HEADER pseudo_type="NC" z_valence="4.0"/></UPF>')
        self.assertEqual(upf.valence_electrons, 4)

    def test_no_guessing_or_rounding(self):
        for header in (b'fixture', b'20 .01 2 nan', b'20 .01 2 3.5', b'20 .01 2 0', b'20 .01 2 15'):
            with self.subTest(header=header), self.assertRaisesRegex(ValueError, 'valence'):
                Pseudopotential(14, 2, '_rps.dat', header).valence_electrons
        with self.assertRaisesRegex(ValueError, 'valence'):
            Pseudopotential(14, 2, '.upf', b'<UPF><PP_HEADER pseudo_type="NC"/></UPF>').valence_electrons

class WarmStartTest(unittest.TestCase):
    def test_wavefunction_only_copy_and_compatibility_before_launch(self):
        import hashlib, json
        from dataclasses import replace
        from types import SimpleNamespace
        from unittest.mock import patch
        from salmonpy.fdtd import SalmonExecutable
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            exe = root/'salmon'; exe.write_bytes(b'test executable')
            engine = SalmonExecutable(exe, hashlib.sha256(exe.read_bytes()).hexdigest())
            model = DFT([14], np.eye(3)*10, [[0,0,0]],
                        [Pseudopotential(14,2,'_rps.dat',b'20 .01 2 4\n')],
                        electrons=4, states=8, grid=(12,)*3, kmesh=(2,)*3)
            def native(*args, **kw):
                work = Path(kw['cwd'])
                kw['stdout'].write('end SALMON\n'); kw['stdout'].flush()
                (work/'gs_status.dat').write_text('SALMONPY_GS_V1\n0 5\n1e-5 1e-9 -1.0\nrho_dne\n')
                (work/'gs_forces.dat').write_text('SALMONPY_FORCES_V1\n1\n0 0 0\n')
                if "write_gs_restart_data='no'" not in kw['input']:
                    restart = work/'data_for_restart'; restart.mkdir()
                    for name in ('info.bin','occupation.bin','wfn.bin'):
                        (restart/name).write_bytes(name.encode())
                return SimpleNamespace(returncode=0)
            with patch('subprocess.run', side_effect=native), self.assertWarns(RuntimeWarning):
                source = model.run(engine, root/'source')
            self.assertFalse(source.converged)
            with patch('subprocess.run',side_effect=native), self.assertWarns(RuntimeWarning):
                scalar=model.run(engine,root/'scalar',save_restart=False)
            self.assertFalse(scalar.restart_saved)
            np.testing.assert_array_equal(scalar.forces_hartree_per_bohr,[[0,0,0]])
            self.assertFalse(list(scalar.directory.rglob('wfn.bin')))
            for operation in (lambda:scalar.band_initial_state(executable=engine),
                              lambda:scalar.tddft_initial_state(executable=engine),
                              lambda:model.run(engine,root/'scalar-guess',initial_guess=scalar)):
                with self.assertRaisesRegex(ValueError,'restart'):operation()
            self.assertFalse((root/'scalar-guess').exists())
            with patch('subprocess.run') as run:
                with self.assertRaisesRegex(TypeError,'save_restart'):model.run(engine,root/'bad-flag',save_restart='no')
                run.assert_not_called()
            self.assertFalse(json.loads((source.directory/'ground-state.json').read_text())['converged'])
            # Older runtimes may ignore the output request: warn, retain results, and record actual state.
            for names in [('info.bin',),('info.bin','occupation.bin','wfn.bin')]:
                def extra_output(*args,**kw):
                    result=native(*args,**kw);restart=Path(kw['cwd'])/'data_for_restart';restart.mkdir()
                    for name in names:(restart/name).write_bytes(b'synthetic')
                    return result
                with patch('subprocess.run',side_effect=extra_output),self.assertWarnsRegex(RuntimeWarning,'runtime wrote restart output'):
                    unexpected=model.run(engine,root/('extra-'+str(len(names))),save_restart=False)
                self.assertEqual(unexpected.restart_saved,len(names)==3)
            legacy_path=source.directory/'ground-state.json'
            original=legacy_path.read_bytes();legacy=json.loads(original);legacy.pop('restart_saved')
            legacy_path.write_text(json.dumps(legacy))
            self.assertTrue(GroundStateResult.open(source.directory).restart_saved)
            legacy_path.write_bytes(original)
            original = (source.directory/'ground-state.json').read_bytes()
            displaced = replace(model, fractional_positions=[[.01,0,0]])
            with patch('subprocess.run', side_effect=native) as run:
                result = displaced.run(engine, root/'warm', initial_guess=source)
                text = run.call_args.kwargs['input']
                for flag in ("yn_restart='y'", "read_gs_restart_data='wfn'", "yn_reset_step_restart='y'", "directory_read_data='./initial_guess/'"):
                    self.assertIn(flag, text)
            self.assertEqual({p.name for p in (result.directory/'native/initial_guess').iterdir()},
                             {'info.bin','occupation.bin','wfn.bin'})
            self.assertEqual(json.loads((result.directory/'ground-state.json').read_text())['initial_guess']['mode'], 'wavefunction')
            GroundStateResult.open(source.directory)
            self.assertEqual((source.directory/'ground-state.json').read_bytes(), original)
            def load_transverse(bundle, name, **options):
                text = (Path(bundle)/name).read_text()
                self.assertIn("trans_longi='tr'", text)
                self.assertIn("yn_restart='n'", text)
                return object()
            with patch('salmonpy.tddft.TDDFTInitialState.from_salmon', side_effect=load_transverse):
                result.tddft_initial_state(executable=engine)
            for change in ({'cell':np.eye(3)*11}, {'grid':(16,)*3}, {'kmesh':(1,)*3},
                           {'states':9}, {'electrons':6}, {'xc':'libxc_pz'},
                           {'pseudopotentials':[Pseudopotential(14,2,'_rps.dat',b'changed')]}):
                with self.subTest(change=change), patch('subprocess.run') as run:
                    with self.assertRaisesRegex(ValueError, 'initial guess'):
                        replace(model, **change).run(engine, root/'rejected', initial_guess=source)
                    run.assert_not_called()
                    self.assertFalse((root/'rejected').exists())
            (source.directory/'native/data_for_restart/wfn.bin').write_bytes(b'tamper')
            with patch('subprocess.run') as run, self.assertRaisesRegex(ValueError, 'checksum'):
                displaced.run(engine, root/'tampered', initial_guess=source)
            run.assert_not_called()
