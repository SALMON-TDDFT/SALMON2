"""Reject incomplete or falsely successful native CI results."""
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from validate_nlcc_ci_results import load_build, validate_test, validate_native_records


class NLCCCIValidation(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(); self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.source = self.root/'source'; self.source.mkdir()
        self.path = self.root/'test-record.json'
        self.base = {'source_root': str(self.source), 'executable_sha256': 'a'*64,
                     'runtime_executable_sha256':'d'*64, 'runtime_source_patch_sha256':'c'*64}

    def record(self, kind='mpi'):
        tests = {'mpi': 'test_nlcc_mpi_force', 'core': 'test_nlcc_core_force', 'xc': 'test_nlcc_xc_variation'}
        result = {**self.base, 'test': tests[kind], 'runs': []}
        for ranks in ((2, 4) if kind == 'mpi' else (1,)):
            for threads in (1, 2):
                if kind == 'mpi': stdout=''.join(f'MPI_NLCC_FORCE_COMPARE {fmt} {ranks} 6e-17\nMPI_NLCC_FORCE_PASS {ranks} 6e-17 {ranks}\n' for fmt in (1,2))
                elif kind == 'core': stdout=''.join(f'CORE_FORCE_VARIATION {g} 6e-10\n' for g in range(1,7))+'NLCC_CORE_FORCE_PASS\n'
                else: stdout=''.join(f'XC_VARIATION {f} {n} 1 1 0\n' for f in ('pz','libxc_pz','libxc_pbe') for n in (0,1))+'NLCC_XC_VARIATION_PASS\n'
                result['runs'].append(dict(ranks=ranks, threads=threads, returncode=0, stdout=stdout))
        return result

    def check(self, data, kind='mpi'):
        self.path.write_text(json.dumps(data))
        pairs=[(r,t) for r in ((2,4) if kind=='mpi' else (1,)) for t in (1,2)]
        return validate_test(self.path,self.source,kind,pairs,runtime_sha256='d'*64,patch_sha256='c'*64)

    def test_complete_records(self):
        for kind in ('xc','core','mpi'):
            self.assertGreater(self.check(self.record(kind),kind)['cases'],0)

    def test_partial_and_duplicate_matrices(self):
        data=self.record()
        for runs in (data['runs'][:-1],data['runs'][:-1]+[data['runs'][0]]):
            with self.assertRaisesRegex(ValueError,'rank/thread'):self.check({**data,'runs':runs})

    def test_mpi_nan_and_noncollective_rejection(self):
        for mode in ('nan','rejected','missing'):
            data=self.record();text=data['runs'][0]['stdout']
            if mode=='nan':text=text.replace('6e-17','nan',1)
            elif mode=='rejected':text=text.replace('MPI_NLCC_FORCE_PASS 2 6e-17 2','MPI_NLCC_FORCE_PASS 2 6e-17 1',1)
            else:text=''
            data['runs'][0]['stdout']=text
            with self.assertRaises(ValueError):self.check(data)

    def test_missing_xc_cases_or_fake_success(self):
        core=self.record('core')
        core['runs'][0]['stdout']=core['runs'][0]['stdout'].replace('CORE_FORCE_VARIATION 6 6e-10','CORE_FORCE_VARIATION 6 1.0')
        with self.assertRaisesRegex(ValueError,'tolerance'):self.check(core,'core')
        data=self.record('xc');data['runs'][0]['stdout']='NLCC_XC_VARIATION_PASS\n'
        with self.assertRaisesRegex(ValueError,'coverage'):self.check(data,'xc')
        data=self.record('xc');data['runs'][0]['stdout']=data['runs'][0]['stdout'].replace('1 1 0','1 1 0.1',1)
        with self.assertRaisesRegex(ValueError,'tolerance'):self.check(data,'xc')

    def test_test_records_are_bound_to_runtime_hash(self):
        data=self.record();data['runtime_executable_sha256']='e'*64
        with self.assertRaisesRegex(ValueError,'binding'):self.check(data)

    def test_normalized_serial_evidence_and_optional_mpi_pair(self):
        exe=self.root/'build/salmon';exe.parent.mkdir();(exe.parent/'CMakeCache.txt').write_text('CMAKE_HOME_DIRECTORY:INTERNAL='+str(self.source)+'\n');exe.write_bytes(b'not executed fixture')
        sha=hashlib.sha256(exe.read_bytes()).hexdigest()
        spec={'upstream_commit':'b'*40,'patch_sha256':'c'*64}
        build=self.root/'build.json'
        build.write_text(json.dumps({'status':'built-and-audited','build_configuration':{'mpi_enabled':False},
            'source':spec,'executable':str(exe),'sha256':sha}))
        for kind in ('xc','core'):
            data=self.record(kind);data['runtime_executable_sha256']=sha
            (self.root/(kind+'.json')).write_text(json.dumps(data))
        result=validate_native_records(build,self.root/'xc.json',self.root/'core.json',source_spec=spec)
        self.assertEqual(result['schema'],'salmonpy.nlcc-native-validation/1')
        self.assertEqual(result['status'],'passed')
        self.assertEqual(result['build_executable_sha256'],sha)
        self.assertEqual(result['tests']['xc']['rank_thread_pairs'],[[1,1],[1,2]])
        with self.assertRaisesRegex(ValueError,'together'):
            validate_native_records(build,self.root/'xc.json',self.root/'core.json',source_spec=spec,mpi_build_record=build)

    def test_build_hash_revision_and_capability(self):
        exe=self.root/'build/salmon';exe.parent.mkdir();(exe.parent/'CMakeCache.txt').write_text('CMAKE_HOME_DIRECTORY:INTERNAL='+str(self.source)+'\n');exe.write_bytes(b'not executed fixture')
        spec={'upstream_commit':'b'*40,'patch_sha256':'c'*64}
        data={'status':'built-and-audited','build_configuration':{'mpi_enabled':False},
              'source':spec,'executable':str(exe),'sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
        path=self.root/'build.json';path.write_text(json.dumps(data))
        load_build(path,spec,False)
        for changed in ({**data,'sha256':'0'*64},{**data,'source':{**spec,'patch_sha256':'0'*64}},
                        {**data,'build_configuration':{'mpi_enabled':True}}):
            path.write_text(json.dumps(changed))
            with self.assertRaises(ValueError):load_build(path,spec,False)

if __name__=='__main__':unittest.main()
