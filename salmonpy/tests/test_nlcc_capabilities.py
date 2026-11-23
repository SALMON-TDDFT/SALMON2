"""Only binary-bound completed validation can unlock NLCC Calculator forces."""
import copy,hashlib,json,tempfile,unittest,importlib.util
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from salmonpy.errors import InputError
from salmonpy.fdtd import SalmonExecutable
from salmonpy import runtime


def receipt(binary_sha):
    from salmonpy._runtime_capabilities import NLCC_PATCH_SHA256, NLCC_SOURCE_TREE_SHA256, reference_sha256
    return {'schema':'salmonpy.nlcc-forces/1','status':'passed','executable_sha256':binary_sha,
        'build_executable_sha256':binary_sha,'source_patch_sha256':NLCC_PATCH_SHA256,'source_tree_sha256':NLCC_SOURCE_TREE_SHA256,
        'scope':{'backend':'cpu','spin':'unpolarized','potential':'NC-UPF','xc':['pz','libxc_pz','libxc_pbe']},
        'native_validation':{'schema':'salmonpy.nlcc-native-validation/1','status':'passed',
            'build_executable_sha256':binary_sha,'source_patch_sha256':NLCC_PATCH_SHA256,'source_tree_sha256':NLCC_SOURCE_TREE_SHA256,
            'tests':{name:{'test':test,'cases':2,'maximum_reported_error':1e-10,'test_executable_sha256':'a'*64,'rank_thread_pairs':[[1,1],[1,2]],'record_sha256':'b'*64}
                for name,test in [('xc','test_nlcc_xc_variation'),('core','test_nlcc_core_force')]}},
        'application_validation':{'status':'passed','executable_sha256':binary_sha,
            'reference_sha256':reference_sha256(),'energy_error_hartree':1e-12,'maximum_force_error_hartree_per_bohr':1e-12}}

class NLCCCapabilities(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup);self.root=Path(self.tmp.name)
        self.exe=self.root/'build/salmon';self.exe.parent.mkdir();self.exe.write_bytes(b'synthetic runtime');self.exe.chmod(0o755)
        self.sha=hashlib.sha256(self.exe.read_bytes()).hexdigest()
    def records(self,token=None):
        from salmonpy._runtime_capabilities import NLCC_PATCH_SHA256, NLCC_SOURCE_TREE_SHA256
        token=receipt(self.sha) if token is None else token
        build={'status':'built-and-audited','executable':str(self.exe),'sha256':self.sha,
               'source':{'patch_sha256':NLCC_PATCH_SHA256,'source_tree_sha256':NLCC_SOURCE_TREE_SHA256},'build_configuration':{'mpi_enabled':False}}
        b=self.root/'build.json';v=self.root/'validation.json';b.write_text(json.dumps(build));v.write_text(json.dumps(token));return b,v
    def test_ordinary_constructor_has_no_nlcc_capability(self):
        engine=SalmonExecutable(self.exe,self.sha)
        self.assertEqual(engine.nlcc_force_xc,())
        with self.assertRaises((TypeError,AttributeError)):SalmonExecutable(self.exe,self.sha,nlcc_force_xc=('libxc_pbe',))
    def test_completed_matching_build_can_enable_and_scope_is_immutable(self):
        engine=runtime.load_validated_build(*self.records())
        self.assertIn('libxc_pbe',engine.nlcc_force_xc)
        before=engine.nlcc_force_xc
        with self.assertRaises((AttributeError,TypeError)):engine.nlcc_force_xc=()
        self.assertEqual(engine.nlcc_force_xc,before)
        from dataclasses import replace
        self.assertEqual(replace(engine).nlcc_force_xc,())
        token=receipt(self.sha)
        token['native_validation']['tests']['mpi']={'test':'test_nlcc_mpi_force','cases':4,'maximum_reported_error':1e-16,'test_executable_sha256':'a'*64,'rank_thread_pairs':[[2,1],[2,2],[4,1],[4,2]],'record_sha256':'b'*64,'build_executable_sha256':self.sha}
        build,proof=self.records(token);data=json.loads(build.read_text());data['build_configuration']['mpi_enabled']=True;build.write_text(json.dumps(data))
        self.assertIn('libxc_pbe',runtime.load_validated_build(build,proof).nlcc_force_xc)
        token['native_validation']['tests'].pop('mpi');proof.write_text(json.dumps(token))
        with self.assertRaises(InputError):runtime.load_validated_build(build,proof)

    def test_wrong_binary_source_failed_or_nonfinite_proofs_are_rejected(self):
        good=receipt(self.sha);bad=[]
        for key,val in [('executable_sha256','0'*64),('source_patch_sha256','0'*64),('status','failed')]:
            t=copy.deepcopy(good);t[key]=val;bad.append(t)
        t=copy.deepcopy(good);t['native_validation']['tests']['xc']['maximum_reported_error']=float('nan');bad.append(t)
        t=copy.deepcopy(good);t['native_validation']['tests'].pop('core');bad.append(t)
        t=copy.deepcopy(good);t['native_validation']['tests']['xc']['rank_thread_pairs']=[[1,1],[1,1]];bad.append(t)
        t=copy.deepcopy(good);t['application_validation']['reference_sha256']='0'*64;bad.append(t)
        t=copy.deepcopy(good);t['application_validation']['maximum_force_error_hartree_per_bohr']=.01;bad.append(t)
        t=copy.deepcopy(good);t['scope']['backend']='gpu';bad.append(t)
        for token in bad:
            with self.subTest(token=token),self.assertRaises(InputError):runtime.load_validated_build(*self.records(token))
        b,v=self.records();self.exe.write_bytes(b'changed')
        with self.assertRaises(InputError):runtime.load_validated_build(b,v)
    def test_unknown_source_tree_cannot_reuse_legacy_patch_identity(self):
        b,v=self.records();data=json.loads(b.read_text())
        missing=copy.deepcopy(receipt(self.sha));missing.pop('source_tree_sha256')
        with self.assertRaises(InputError):runtime.load_validated_build(*self.records(missing))
        b,v=self.records()
        data['source']['source_tree_sha256']='0'*64;b.write_text(json.dumps(data))
        with self.assertRaises(InputError):runtime.load_validated_build(b,v)
        from salmonpy._runtime_capabilities import NLCC_SOURCE_TREE_SHA256
        data['source']['source_tree_sha256']=NLCC_SOURCE_TREE_SHA256;b.write_text(json.dumps(data))
        self.assertIn('libxc_pbe',runtime.load_validated_build(b,v).nlcc_force_xc)
    def test_bundle_feature_without_valid_receipt_is_not_trusted(self):
        from salmonpy._runtime_capabilities import NLCC_PATCH_SHA256, NLCC_SOURCE_TREE_SHA256
        bundle=self.root/'bundle';(bundle/'bin').mkdir(parents=True);exe=bundle/'bin/salmon';exe.write_bytes(b'fixture');exe.chmod(0o755)
        sha=hashlib.sha256(exe.read_bytes()).hexdigest();host=runtime.host_platform()
        meta={'schema':'salmonpy.runtime/1','interface':runtime.INTERFACE,'platform':host,
              'features':['dft','libxc','nlcc-forces'],'executable':'bin/salmon',
              'provenance':{'patch_sha256':NLCC_PATCH_SHA256,'source_tree_sha256':NLCC_SOURCE_TREE_SHA256},'files':{'bin/salmon':sha}}
        (bundle/'runtime.json').write_text(json.dumps(meta))
        with self.assertRaises(InputError):runtime.load_runtime(bundle)
        token=receipt(sha);(bundle/'validation.json').write_text(json.dumps(token))
        meta['validated_capabilities']={'nlcc-forces':'validation.json'}
        meta['files']['validation.json']=hashlib.sha256((bundle/'validation.json').read_bytes()).hexdigest()
        (bundle/'runtime.json').write_text(json.dumps(meta))
        self.assertIn('libxc_pbe',runtime.load_runtime(bundle).nlcc_force_xc)
        proof=runtime.prepare_runtime_job(bundle,self.root/'job-proof.json',job_id='nlcc-test',job_token='n'*64)
        with patch.object(runtime,'_hash',wraps=runtime._hash) as hashes:
            loaded=runtime.load_runtime(bundle,validation=proof)
        self.assertIn('libxc_pbe',loaded.nlcc_force_xc)

    @unittest.skipUnless(importlib.util.find_spec('ase'),'optional ASE dependency')
    def test_validated_engine_can_supply_nlcc_forces_to_calculator(self):
        from ase import Atoms
        from salmonpy.ase_calculator import SalmonCalculator
        from salmonpy.dft import Pseudopotential,DFT
        engine=runtime.load_validated_build(*self.records())
        pp=Pseudopotential(14,2,'.upf',b'<UPF><PP_HEADER pseudo_type="NC" z_valence="4" core_correction="T"/></UPF>')
        calc=SalmonCalculator(executable=engine,directory=self.root/'out',pseudopotentials=[pp],states=8,grid=(12,)*3,kmesh=(2,)*3,xc='libxc_pbe')
        atoms=Atoms('Si',cell=np.eye(3)*5,pbc=True);atoms.calc=calc
        with patch.object(DFT,'run',return_value=SimpleNamespace(energy_hartree=-1.,forces_hartree_per_bohr=np.zeros((1,3)))) as run:
            atoms.get_potential_energy();np.testing.assert_array_equal(atoms.get_forces(),np.zeros((1,3)));self.assertEqual(run.call_count,1)
