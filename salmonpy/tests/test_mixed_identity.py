import copy
import unittest
from types import SimpleNamespace
from salmonpy.distributed import ParallelPlan,SolverGroup
from salmonpy.material_response import _point_metadata
from salmonpy.examples._mixed_identity import sample_records,checkpoint_identity,check_identity


class MixedIdentityTests(unittest.TestCase):
    def fixture(self):
        points={10:(19,20,20),20:(20,20,20)}
        bindings={10:SimpleNamespace(key='sbe-point',material_record=None),
                  20:SimpleNamespace(key='si-tddft',material_record=None)}
        signature={'sbe_model':'a'*64,'library':'b'*64,'assets':{'rt.inp':'c'*64,'data_for_restart/wfn.bin':'d'*64}}
        record,mapping,ld=sample_records('run-a',signature,points,bindings,
            constitutive={'drude':{'plasma_frequency':20}},cell_ids={'1':'drude'})
        plan=ParallelPlan((SolverGroup('field',(0,),role='field'),SolverGroup('sbe',(1,),point_ids=(10,)),SolverGroup('tddft',(2,),point_ids=(20,))))
        metadata=_point_metadata(plan,record,mapping)
        coupled=SimpleNamespace(response=SimpleNamespace(context=SimpleNamespace(plan=plan),identity_record=lambda:copy.deepcopy(metadata)),constitutive_identity=ld)
        return coupled,record,mapping,ld

    def test_full_identity_is_detached_and_matching(self):
        coupled,record,mapping,ld=self.fixture()
        current=checkpoint_identity(coupled)
        check_identity(current,copy.deepcopy(current),'salmonpy.native-direct/4')
        self.assertEqual(current['point_metadata']['point_run']['run_id'],'run-a')
        self.assertEqual(current['constitutive']['cell_ids'],{'1':'drude'})
        current['constitutive']['materials']['drude']['material_id']='changed'
        self.assertEqual(ld['materials']['drude']['material_id'],'drude')

    def test_identity_mismatches_and_json_types_are_rejected(self):
        coupled,_,_,_=self.fixture();current=checkpoint_identity(coupled)
        for field in ('run_id','point_id','material_id','definition_digest','initial_state_id','initial_state_digest','ld','ownership','bool'):
            saved=copy.deepcopy(current)
            r=saved['point_metadata']['point_run']
            if field=='run_id':r['run_id']='other-run'
            elif field=='ld':saved['constitutive']['cell_ids']['1']='other-material'
            elif field=='ownership':saved['point_metadata']['mpi_point_map'][0]['ranks']=[2]
            elif field=='bool':saved['point_metadata']['mpi_point_map'][0]['ranks']=[True]
            else:r['points'][0][field]='e'*64
            with self.subTest(field=field),self.assertRaisesRegex(ValueError,'identity'):
                check_identity(current,saved,'salmonpy.native-direct/4')

    def test_current_schema_and_optional_identity(self):
        coupled,_,_,_=self.fixture();current=checkpoint_identity(coupled)
        check_identity(None,None,'salmonpy.native-direct/4')
        for schema in ('salmonpy.native-direct/1','salmonpy.native-direct/2','salmonpy.native-direct/3','salmonpy.native-replay/1'):
            with self.subTest(schema=schema),self.assertRaisesRegex(ValueError,'schema'):
                check_identity(None,None,schema)
        with self.assertRaises(ValueError):check_identity(current,None,'salmonpy.native-direct/4')
        with self.assertRaises(ValueError):check_identity(None,current,'salmonpy.native-direct/4')

    def test_no_implicit_identity_for_legacy_session(self):
        self.assertIsNone(checkpoint_identity(SimpleNamespace(response=SimpleNamespace())))

    def test_explicit_constitutive_identity_and_replay(self):
        explicit={'material_id':'custom','definition_digest':'e'*64,
                  'initial_state_id':'prepared','initial_state_digest':'f'*64}
        _,_,ld=sample_records('run',{}, {}, {},constitutive={'custom':{}},
            cell_ids={'1':'custom'},constitutive_records={'custom':explicit})
        self.assertEqual(ld['materials']['custom'],explicit)
        coupled,_,_,_=self.fixture();current=checkpoint_identity(coupled)
        check_identity(current,copy.deepcopy(current),'salmonpy.native-replay/2')
        with self.assertRaises(ValueError):check_identity(None,current,'salmonpy.native-replay/2')

    def test_separate_hse_bundle_changes_only_hse_identity(self):
        points={20:(20,20,20),30:(21,20,20)}
        bindings={20:SimpleNamespace(key='si-tddft',material_record=None),
                  30:SimpleNamespace(key='hse-point',material_record=None)}
        signature={'library':'b'*64,'assets':{'rt.inp':'c'*64},
                   'hse_assets':{'rt.inp':'d'*64},'fast_hse_backend':'hse-v1'}
        first=sample_records('run',signature,points,bindings,constitutive={},cell_ids={})[0]
        signature['hse_assets']['rt.inp']='e'*64
        second=sample_records('run',signature,points,bindings,constitutive={},cell_ids={})[0]
        self.assertEqual(first['points'][0],second['points'][0])
        self.assertNotEqual(first['points'][1]['definition_digest'],second['points'][1]['definition_digest'])
        self.assertNotEqual(first['points'][1]['initial_state_digest'],second['points'][1]['initial_state_digest'])

    def test_catalog_assets_identify_each_material_independently(self):
        from salmonpy.examples._mixed_identity import digest_record
        points={1:(19,20,20),2:(20,20,20)}
        bindings={k:SimpleNamespace(key=mid,material_record=None) for k,mid in ((1,'pbe'),(2,'tbmbj'))}
        signature={'library':'b'*64,'assets':{'rt.inp':'a'*64},
                   'standard_material_inputs':{'pbe':{'assets':{'rt.inp':'c'*64}},
                                               'tbmbj':{'assets':{'rt.inp':'d'*64}}}}
        first=sample_records('run',signature,points,bindings,constitutive={},cell_ids={})[0]
        for point,mid in zip(first['points'],('pbe','tbmbj')):
            assets=signature['standard_material_inputs'][mid]['assets']
            self.assertEqual(point['initial_state_digest'],digest_record({'assets':assets}))
        signature['standard_material_inputs']['tbmbj']['assets']['rt.inp']='e'*64
        second=sample_records('run',signature,points,bindings,constitutive={},cell_ids={})[0]
        self.assertEqual(first['points'][0],second['points'][0])
        self.assertNotEqual(first['points'][1]['definition_digest'],second['points'][1]['definition_digest'])
        self.assertNotEqual(first['points'][1]['initial_state_digest'],second['points'][1]['initial_state_digest'])
