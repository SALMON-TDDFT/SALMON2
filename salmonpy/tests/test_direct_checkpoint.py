"""Admission failures must precede native mutation."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from salmonpy.examples._direct_checkpoint import checkpoint, digest


class DirectCheckpointTests(unittest.TestCase):
    def world(self):
        return SimpleNamespace(rank=0, size=1, allgather=lambda x:[x], Barrier=lambda:None)

    def coupled(self):
        raw=SimpleNamespace(snapshot=lambda:None)
        response=SimpleNamespace(states={1:SimpleNamespace(_raw=raw)}, failed=False)
        return SimpleNamespace(response=response, field=None)

    def test_unsupported_backend_before_directory_creation(self):
        c=self.coupled();c.response.states={1:object()}
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'checkpoint'
            with self.assertRaisesRegex(ValueError,'capability'):checkpoint(self.world(),c,path,{})
            self.assertFalse(path.exists())

    def test_incomplete_rejected_before_state_change(self):
        with tempfile.TemporaryDirectory() as d:
            c=self.coupled()
            with self.assertRaises(ValueError):checkpoint(self.world(),c,d,{},restore=True)
            self.assertFalse(c.response.failed)

    def test_changed_signature_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'complete.json').write_text(json.dumps({'schema':'salmonpy.native-direct/1','signature':{'old':1},'ranks':1}))
            with self.assertRaisesRegex(ValueError,'mismatch'):checkpoint(self.world(),self.coupled(),p,{'steps':128},restore=True)

    def test_corrupted_rank_file_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);data=p/'rank-0.json';data.write_text('{}')
            m={'schema':'salmonpy.native-direct/1','signature':{'steps':128},'step':64,'ranks':1,'files':{data.name:digest(data), 'frames.json':'bad'}}
            (p/'complete.json').write_text(json.dumps(m));data.write_text('changed')
            with self.assertRaisesRegex(ValueError,'checksum'):checkpoint(self.world(),self.coupled(),p,{'steps':128},restore=True)
