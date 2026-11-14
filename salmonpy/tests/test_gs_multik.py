"""Opt-in multi-k contracts with synthetic protocol/native-file fixtures."""
from pathlib import Path
import hashlib,json,unittest
from unittest.mock import patch
import numpy as np
import test_gs_native as native_fixture
import test_gs_export as export_fixture
from salmonpy.gs_native import SerialRestartConverter,GammaRestartConverter
from salmonpy.gs_archive import save_ground_state
from salmonpy.errors import InputError


def k_text(points,weights,cell):
    lines=['# k-point distribution','# ik: k-point index','# kx,ky,kz: Reduced coordinate of k-points',
           '# wk: Weight of k-point','# 1:ik[none] 2:kx[none] 3:ky[none] 4:kz[none] 5:wk[none]']
    lines += [str(k+1)+' '+' '.join(format(float(v),'.15e') for v in (*point,weights[k])) for k,point in enumerate(points)]
    lines += ['# coefficients (2*pi/a [a.u.]) in kx, ky, kz',' '.join(format(v,'.15e') for v in 2*np.pi/np.diag(cell))]
    return '\n'.join(lines)+'\n'


class GSMultiKImport(unittest.TestCase):
    def setUp(self):
        native_fixture.GSNativeTests.setUp(self)
        record=json.loads(self.converter.build_record.read_text());record['multi_k_byte_parity']=True
        self.converter.build_record.write_text(json.dumps(record));self.serial=SerialRestartConverter.from_build(self.converter.build_record)
        self.waves=np.arange(1,2049).reshape(2,2,8,8,8)/4096.+1j*np.arange(2048).reshape(2,2,8,8,8)/8192.
        self.kw=dict(self.values,wavefunctions=self.waves,occupations=[[2.,2.],[2.,2.]],
                     eigenvalues=[[-.5,.2],[-.4,.3]],kpoints=[[-.25,0,0],[.25,0,0]],weights=[.25,.75],electrons=4.)
        self.multi=save_ground_state(Path(self.tmp.name)/'multi',orbital_block=1,**self.kw)

    def test_capability_gated_and_gamma_stays_restricted(self):
        with self.assertRaises(InputError):self.converter.write_bundle(self.multi,Path(self.tmp.name)/'gamma',executable=self.model.executable,pseudopotentials=[self.potential])
        record=json.loads(self.converter.build_record.read_text());record.pop('multi_k_byte_parity');self.converter.build_record.write_text(json.dumps(record))
        with self.assertRaises(InputError):SerialRestartConverter.from_build(self.converter.build_record)
        GammaRestartConverter.from_build(self.converter.build_record)

    def test_order_weights_file_and_receipt(self):
        from salmonpy.tddft import TDDFTInitialState
        dest=self.serial.write_bundle(self.multi,Path(self.tmp.name)/'native-multi',executable=self.model.executable,pseudopotentials=[self.potential])
        self.assertEqual((dest/'data_for_restart/wfn.bin').read_bytes(),np.asarray(self.waves,dtype='<c16').tobytes())
        self.assertIn("file_kw='kpoints.dat'",(dest/'rt.inp').read_text())
        values=[line.split() for line in (dest/'kpoints.dat').read_text().splitlines()]
        self.assertEqual(values[0],['2']);self.assertEqual([float(row[-1]) for row in values[1:]],[.25,.75])
        initial=TDDFTInitialState.from_salmon(dest,'rt.inp',executable=self.model.executable);initial.close()
        (dest/'kpoints.dat').write_text('changed')
        with self.assertRaises(InputError):TDDFTInitialState.from_salmon(dest,'rt.inp',executable=self.model.executable)

    def test_nonorthogonal_refused_before_conversion(self):
        kw=dict(self.kw,cell_bohr=[[4,.1,0],[0,4,0],[0,0,4]])
        other=save_ground_state(Path(self.tmp.name)/'nonorthogonal',**kw)
        with patch('salmonpy.gs_native.subprocess.run',side_effect=AssertionError('must not execute')):
            with self.assertRaises(InputError):self.serial.write_bundle(other,Path(self.tmp.name)/'bad-cell',executable=self.model.executable,pseudopotentials=[self.potential])


class GSMultiKExport(unittest.TestCase):
    def setUp(self):
        export_fixture.GSExportTests.setUp(self)
        from salmonpy.dft import DFT
        record=json.loads(self.converter.build_record.read_text());record['multi_k_byte_parity']=True
        self.converter.build_record.write_text(json.dumps(record));self.serial=SerialRestartConverter.from_build(self.converter.build_record)
        self.dft=DFT([2],np.eye(3)*4,[[0,0,0]],[self.potential],4,2,(8,8,8),(2,1,1),'pz')
        native=self.source/'native';(native/'gs.inp').write_text(self.dft.input_text())
        self.waves=np.arange(1,2049).reshape(2,2,8,8,8)/4096.+1j*np.arange(2048).reshape(2,2,8,8,8)/8192.
        (native/'data_for_restart/wfn.bin').write_bytes(np.asarray(self.waves,dtype='<c16').tobytes())
        (native/'data_for_restart/occupation.bin').write_bytes(np.full((2,2),2.,dtype='<f8').tobytes())
        (native/'gs_eigen.data').write_text('#esp: single-particle energies (eigen energies)\n#occ: occupation numbers, io: orbital index\n# 1:io, 2:esp[a.u.], 3:occ\nk= 1,  spin= 1\n1 -.5 2\n2 .2 2\nk= 2, spin= 1\n1 -.4 2\n2 .3 2\n')
        self.kt=k_text([[-.25,0,0],[.25,0,0]],[.5,.5],self.dft.cell)
        (native/'gs_k.data').write_text(self.kt)
        self.record.update(self.dft._restart_signature(),input_sha256=hashlib.sha256((native/'gs.inp').read_bytes()).hexdigest())
        self.rehash()
    rehash=export_fixture.GSExportTests.rehash

    def test_all_k_order_and_weighted_density(self):
        a=self.result.export_archive(Path(self.tmp.name)/'multi-export',converter=self.serial,orbital_block=1)
        for k in range(2):self.assertEqual(a.wavefunctions(k).tobytes(),self.waves[k].tobytes())
        np.testing.assert_array_equal(a.array('weights'),[.5,.5])
        np.testing.assert_allclose(a.array('density'),np.sum(abs(self.waves)**2,axis=(0,1)),rtol=1e-15)
        self.assertEqual(a.metadata['provenance']['kpoint_source'],'native/gs_k.data')

    def test_k_table_must_match_generated_source_mesh(self):
        p=self.source/'native/gs_k.data'
        wrong_point=self.kt.replace('-2.500000000000000e-01','-1.250000000000000e-01')
        wrong_weights=k_text([[-.25,0,0],[.25,0,0]],[.25,.75],self.dft.cell)
        for text in (wrong_point,wrong_weights):
            p.write_text(text);self.rehash()
            with self.assertRaises(InputError):self.result.export_archive(Path(self.tmp.name)/'wrong-mesh',converter=self.serial)

    def test_k_order_basis_and_eigen_order_rejected(self):
        p=self.source/'native/gs_k.data';e=self.source/'native/gs_eigen.data';original=e.read_text()
        for bad in (self.kt.replace('1 -2.500','2 -2.500'),self.kt.replace('coefficients','wrong-basis'),self.kt.replace('Reduced coordinate','Cartesian coordinate')):
            p.write_text(bad);self.rehash()
            with self.assertRaises(InputError):self.result.export_archive(Path(self.tmp.name)/'bad',converter=self.serial)
        p.write_text(self.kt);e.write_text(original.replace('k= 2,','k= 1,'));self.rehash()
        with self.assertRaises(InputError):self.result.export_archive(Path(self.tmp.name)/'bad-eigen',converter=self.serial)

if __name__=='__main__':unittest.main()
