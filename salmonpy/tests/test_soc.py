from dataclasses import replace
import unittest
from pathlib import Path
import numpy as np
from salmonpy.dft import DFT, Pseudopotential
from salmonpy.errors import InputError


def potential(z, charge, soc=True):
    attributes = 'relativistic="full" has_so="T"' if soc else 'relativistic="scalar" has_so="F"'
    data = (f'<UPF version="2.0.1"><PP_HEADER pseudo_type="NC" is_paw="F" '
            f'is_ultrasoft="F" {attributes} z_valence="{charge}" functional="PBE"/>'
            + ('<PP_SPIN_ORB/>' if soc else '') + '</UPF>').encode()
    return Pseudopotential(z, 2, '.upf', data)


class SOCInputTests(unittest.TestCase):
    def model(self):
        return DFT([74,34,34],np.diag([6.2,6.2,40.]),[[0,0,.5],[1/3,2/3,.42],[1/3,2/3,.58]],
                   [potential(74,14),potential(34,16)],electrons=46,states=92,
                   grid=[16,16,64],kmesh=[2,2,1],xc='libxc_pbe',spin_orbit=True)

    def test_material_validation_preparation_never_initializes_mpi(self):
        import importlib.util
        import io
        import json
        import tempfile
        from contextlib import redirect_stdout
        from unittest.mock import patch
        from salmonpy.pseudopotentials import PseudopotentialReference
        tool=Path(__file__).resolve().parents[2]/'tools/validate_soc_material.py'
        spec=importlib.util.spec_from_file_location('soc_material_validation_test',tool)
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp);refs={}
            entries=[]
            for symbol,z,charge in [('W',74,14),('Se',34,16)]:
                pp=potential(z,charge)
                (root/(symbol+'.upf')).write_bytes(pp.data)
                refs[symbol]=PseudopotentialReference(symbol+'.upf',z,'https://example.org/'+symbol+'.upf',
                    pp.sha256,'test','synthetic fixture')
                entries.append(dict(catalog_id=symbol,file=symbol+'.upf',local_channel=2))
            model=self.model()
            config=dict(schema='salmonpy.soc-validation/1',material_id='synthetic-wse2',
                pseudopotentials=entries,numbers=list(model.numbers),cell_bohr=model.cell.tolist(),
                fractional_positions=model.fractional_positions.tolist(),electrons=46,states=92,
                grid=[16,16,64],kmesh=[2,2,1],xc='libxc_pbe',temperature_k=300,
                scf=dict(max_iterations=400,tolerance=1e-8),response=dict(steps=16,dt=.001,A=[1e-5,0,0]))
            path=root/'input.json';path.write_text(json.dumps(config));output=root/'must-not-be-created'
            stdout=io.StringIO()
            import sys
            with patch('salmonpy.pseudopotentials.catalog_reference',side_effect=lambda name:refs[name]), \
                 patch.dict(sys.modules,{'mpi4py':None}), \
                 patch.object(sys,'argv',['validate_soc_material.py','missing-build.json',str(path),str(output)]), \
                 redirect_stdout(stdout):
                module.main()
            self.assertFalse(output.exists())
            self.assertIn('temperature_k=',stdout.getvalue())
            self.assertIn('nstate=92',stdout.getvalue())
            self.assertIn('no MPI initialization or calculation',stdout.getvalue())

    def test_shared_tddft_spin_pair_and_xc_contract(self):
        import tempfile
        from types import SimpleNamespace
        from salmonpy._shared_tddft import _profile
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            initial=SimpleNamespace(directory=root,input_file='rt.inp')
            text=self.model().input_text()
            (root/'rt.inp').write_text(text)
            _profile(initial)
            for changed in (text.replace("spin='noncollinear'","spin='unpolarized'"),
                            text.replace("yn_spinorbit='y'","yn_spinorbit='n'"),
                            text.replace("xc='libxc_pbe'","xc='tbmbj'")):
                (root/'rt.inp').write_text(changed)
                with self.assertRaises(InputError):_profile(initial)

    def test_soc_input_and_restart_identity(self):
        model=self.model()
        self.assertIn("spin='noncollinear'",model.input_text())
        self.assertIn("yn_spinorbit='y'",model.input_text())
        self.assertIn('num_kgrid=2,2,1',model.input_text())
        self.assertNotIn('yn_gamma_centered=',model.input_text())
        self.assertNotIn('dk_shift=',model.input_text())
        self.assertEqual(model._restart_signature()['spin_orbit'],True)
        self.assertTrue(model.pseudopotentials[0].has_spin_orbit)
        self.assertEqual(model.pseudopotentials[1].valence_electrons,16)
        with self.assertRaisesRegex(InputError,'states'):
            replace(model,states=23)
        with self.assertRaisesRegex(InputError,'spin_orbit'):
            replace(model,spin_orbit=False)
        with self.assertRaisesRegex(InputError,'spin_orbit'):
            replace(model,pseudopotentials=[potential(74,14,False),potential(34,16,False)])

    def test_scalar_input_stays_scalar(self):
        p=potential(14,4,False)
        model=DFT([14],np.eye(3)*10,[[0,0,0]],[p],4,8,[16]*3,[1]*3)
        self.assertIn("spin='unpolarized'",model.input_text())
        self.assertNotIn("yn_spinorbit='y'",model.input_text())
        self.assertFalse(p.has_spin_orbit)
        self.assertNotIn('spin_orbit',model._restart_signature())

    def test_catalog_never_substitutes_scalar_for_soc(self):
        from salmonpy.pseudopotentials import list_references, catalog_reference
        for element, charge in [('W',14),('Se',16)]:
            entries=list_references(element=element,xc='PBE',relativity='fully-relativistic')
            self.assertEqual(len(entries),1)
            reference=entries[0]
            self.assertTrue(reference.has_so)
            self.assertEqual(reference.valence_electrons,charge)
            self.assertEqual(catalog_reference(reference.catalog_id),reference)
            scalar=list_references(element=element,xc='PBE',relativity='scalar-relativistic')
            self.assertEqual(len(scalar),1)
            self.assertNotEqual(scalar[0].sha256,reference.sha256)
        self.assertEqual(list_references(element='W',xc='PZ',relativity='fully-relativistic'),())
        with self.assertRaises(InputError):catalog_reference('unknown')

    def test_pinned_reference_is_bound_to_potential_bytes(self):
        from salmonpy.pseudopotentials import PseudopotentialReference
        p=potential(74,14)
        ref=PseudopotentialReference('W.upf',74,'https://example.org/W.upf',p.sha256,'CC-BY-4.0','example')
        bound=replace(p,reference=ref)
        self.assertEqual(bound.reference.sha256,bound.sha256)
        with self.assertRaisesRegex(InputError,'identity'):
            replace(bound,data=p.data+b'changed')

    def test_soc_ase_preparation_never_runs_solver(self):
        import importlib.util, tempfile
        if importlib.util.find_spec('ase') is None:self.skipTest('optional ASE')
        from ase import Atoms
        from salmonpy.ase import dft_from_atoms
        from salmonpy.preparation import prepare_dft
        atoms=Atoms('WSe2',cell=[3.28,3.28,20.],scaled_positions=[[0,0,.5],[1/3,2/3,.58],[1/3,2/3,.42]],pbc=True)
        model=dft_from_atoms(atoms,pseudopotentials=[potential(74,14),potential(34,16)],
                             states=92,kmesh=(2,2,1),grid=(24,24,128),xc='libxc_pbe',spin_orbit=True)
        self.assertEqual(model.electrons,46)
        with tempfile.TemporaryDirectory() as tmp:
            from pathlib import Path
            prepared=prepare_dft(model,Path(tmp)/'prepared')
            self.assertIn("yn_spinorbit='y'",(prepared/'gs.inp').read_text())
            self.assertFalse((prepared/'stdout.log').exists())

    def test_nickel_metal_and_scf_conditions_are_explicit(self):
        from salmonpy.dft import SCFSettings
        from salmonpy.pseudopotentials import catalog_reference
        reference=catalog_reference('pseudodojo-pbe-0.4-fr-standard-ni')
        self.assertTrue(reference.has_so)
        self.assertEqual(reference.valence_electrons,18)
        self.assertEqual(reference.semicore_shells,('3s','3p'))
        self.assertEqual(reference.valence_configuration,'3s2 3p6 3d8 4s2')
        p=potential(28,18)
        model=DFT([28],np.array([[0,3.3,3.3],[3.3,0,3.3],[3.3,3.3,0]]),[[0,0,0]],[p],
                  18,36,[16]*3,[2]*3,xc='libxc_pbe',spin_orbit=True,electronic_temperature_k=300)
        scf=SCFSettings(max_iterations=3,mixing_method='simple_dm',mixing_rate=.2,
                        initial_wavefunction='random',random_seed=0)
        text=model.input_text(scf)
        for fragment in ("method_mixing='simple_dm'", "method_init_wf='random'",'iseed_number_change=0',
                         'temperature_k=3.00000000000000000e+02','num_kgrid=2,2,2'):
            self.assertIn(fragment,text)
        self.assertEqual(model._restart_signature()['electronic_temperature_k'],300)
        self.assertNotIn('dk_shift=',text)

    def test_soc_default_states_and_kelvin_temperature(self):
        model=replace(self.model(),states=None,electronic_temperature_k=300)
        self.assertEqual(model.states,2*model.electrons)
        self.assertIn('temperature_k=',model.input_text())
        with self.assertRaises(InputError):
            replace(model,electronic_temperature=.001)
