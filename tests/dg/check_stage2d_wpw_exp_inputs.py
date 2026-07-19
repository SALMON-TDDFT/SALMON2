#!/usr/bin/env python3
from pathlib import Path
import hashlib

root=Path(__file__).resolve().parents[2]
case=root/'samples/exercise_dg_fragment_rt/diamond64_dc_flux_mac'
names=['inputfile_stage2d_dg_dc_wpw_smoke','inputfile_stage2d_wpw_fieldoff_smoke',
       'inputfile_stage2d_wpw_laser_smoke']
manifest=case/'stage2d_wpw_exp_manifest.tsv'
assert manifest.is_file(), 'missing manifest'
rows={}
for line in manifest.read_text().splitlines():
    if not line or line.startswith('#'): continue
    key,value=line.split('\t',1); rows[key]=value
for name in names:
    path=case/name
    assert path.is_file(), f'missing {name}'
    digest=hashlib.sha256(path.read_bytes()).hexdigest()
    assert rows[f'sha256:{name}']==digest, f'hash mismatch {name}'
gs=(case/names[0]).read_text().lower()
fieldoff=(case/names[1]).read_text().lower()
laser=(case/names[2]).read_text().lower()
for token in ("yn_dg_wpw_production = 'y'",'num_fragment(1:3) = 2,2,2',
              'n_plane_waves_dg = 128','k_cutoff_plane_wave = 2.0d0'):
    assert token in gs, token
for text in (fieldoff,laser):
    for token in ("yn_dg_wpw_checkpoint_rt = 'y'","yn_dg_fragment_rt = 'y'",
                  "yn_dg_length_gauge = 'y'",'dg_wpw_checkpoint_manifest'):
        assert token in text, token
assert "ae_shape1 = 'none'" in fieldoff and 'nt = 20' in fieldoff
assert "ae_shape1 = 'asin2cos'" in laser and 'nt = 20' in laser
assert rows['species']=='Si' and rows['natom']=='64' and rows['nelec']=='256'
assert rows['grid']=='24,24,24' and rows['fragments']=='2,2,2'
pseudo=root/rows['pseudopotential']
atom=root/rows['atom_file']
assert hashlib.sha256(pseudo.read_bytes()).hexdigest()==rows['sha256:pseudopotential']
assert hashlib.sha256(atom.read_bytes()).hexdigest()==rows['sha256:atom_file']
print('PASS stage2d Si64 WPW production input provenance')
