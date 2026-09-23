import json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from checkpoint import save_checkpoint,load_checkpoint,fingerprint

class CheckpointTests(unittest.TestCase):
 def test_roundtrip_and_guards(self):
  with tempfile.TemporaryDirectory() as folder:
   p=Path(folder)/'restart.npz';u=np.ones((2,1,3),complex);g=np.ones((2,1,1),complex)
   metadata=dict(format='SALMON_PT_HSE_RESTART_V1',step=1,dt=.32,amplitude=.0001,export_hash='abc',initial_hash='def',initial_energy=-1.)
   rows=[dict(step=1,time_au=.32)]
   save_checkpoint(p,u,g,metadata,rows)
   a,b,m,r,previous=load_checkpoint(p,dict(dt=.32,amplitude=.0001,export_hash='abc',initial_hash='def'))
   np.testing.assert_array_equal(a,u);self.assertEqual(m,metadata);self.assertEqual(r,rows);np.testing.assert_array_equal(previous,u.transpose(0,2,1))
   save_checkpoint(p,u,g,metadata,rows,previous=2*previous)
   np.testing.assert_array_equal(load_checkpoint(p,{})[-1],2*previous)
   with self.assertRaises(ValueError):load_checkpoint(p,dict(dt=.16))
   with self.assertRaises(ValueError):save_checkpoint(p,u,g,metadata,[])
   with self.assertRaises(ValueError):save_checkpoint(p,u*np.nan,g,metadata,rows)
 def test_failed_write_preserves_previous_checkpoint(self):
  with tempfile.TemporaryDirectory() as folder:
   p=Path(folder)/'restart.npz';p.write_bytes(b'previous')
   m=dict(format='SALMON_PT_HSE_RESTART_V1',step=0,dt=.32,amplitude=0.,export_hash='a',initial_hash='b',initial_energy=-1.)
   with patch('checkpoint.np.savez_compressed',side_effect=OSError('disk failure')):
    with self.assertRaises(OSError):save_checkpoint(p,np.ones((1,1,2)),np.ones((1,1,1)),m,[])
   self.assertEqual(p.read_bytes(),b'previous')
   self.assertEqual(list(Path(folder).iterdir()),[p])
 def test_hash_tracks_content(self):
  with tempfile.TemporaryDirectory() as folder:
   p=Path(folder)/'a';p.write_bytes(b'a');first=fingerprint([p]);p.write_bytes(b'b')
   self.assertNotEqual(first,fingerprint([p]))

if __name__=='__main__':unittest.main()
