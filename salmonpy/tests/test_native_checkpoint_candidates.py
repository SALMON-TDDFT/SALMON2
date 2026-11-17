import unittest
from types import SimpleNamespace
from salmonpy._native_checkpoint import NativeRestoreCandidate
from salmonpy.errors import InputError, RuntimeFailure


class Function:
    def __init__(self, call):
        self.call = call
    def __call__(self, *args):
        return self.call(*args)


class NativeCandidateTests(unittest.TestCase):
    def owner(self, status=0):
        calls=[]
        fn=Function(lambda path,mode: calls.append((path,mode)) or (status if mode==2 else 0))
        owner=SimpleNamespace(lib=SimpleNamespace(checkpoint=fn,checkpoint_protocol=Function(lambda:1)),
                              _call=lambda f,*args:f(*args))
        return owner,calls

    def test_stage_discard_and_no_destructive_fallback(self):
        owner,calls=self.owner()
        candidate=NativeRestoreCandidate(owner,'checkpoint','saved.bin')
        self.assertIs(candidate.stage(),candidate)
        candidate.discard();candidate.discard()
        self.assertEqual(calls,[(b'saved.bin',2),(b'',4)])
        with self.assertRaises(InputError):candidate.publish()

    def test_stage_failure_can_cleanup_without_publication(self):
        owner,calls=self.owner(7)
        candidate=NativeRestoreCandidate(owner,'checkpoint','broken.bin')
        with self.assertRaises(RuntimeFailure):candidate.stage()
        candidate.discard()
        self.assertEqual([m for p,m in calls],[2,4])

    def test_publish_does_not_reopen_archive(self):
        owner,calls=self.owner()
        candidate=NativeRestoreCandidate(owner,'checkpoint','saved.bin')
        candidate.stage();candidate.publish()
        self.assertEqual(calls,[(b'saved.bin',2),(b'',3)])
        with self.assertRaises(InputError):candidate.publish()

    def test_old_runtime_rejected_before_checkpoint(self):
        owner,calls=self.owner();del owner.lib.checkpoint_protocol
        with self.assertRaises(InputError):NativeRestoreCandidate(owner,'checkpoint','saved.bin')
        self.assertEqual(calls,[])

    def test_second_candidate_cannot_discard_first(self):
        owner,calls=self.owner()
        first=NativeRestoreCandidate(owner,'checkpoint','first.bin')
        second=NativeRestoreCandidate(owner,'checkpoint','second.bin')
        first.stage()
        with self.assertRaises(InputError):second.stage()
        second.discard()
        self.assertEqual(calls,[(b'first.bin',2)])
        first.discard()

    def test_stage_rejects_step_backup_without_reserving_or_native_call(self):
        owner,calls=self.owner();owner._step_guard=object()
        candidate=NativeRestoreCandidate(owner,'checkpoint','saved.bin')
        with self.assertRaises(InputError):candidate.stage()
        self.assertEqual(calls,[])
        self.assertIsNone(getattr(owner,'_restore_candidate',None))
