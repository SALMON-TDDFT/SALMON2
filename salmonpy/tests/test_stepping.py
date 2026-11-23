import unittest
from salmonpy.stepping import StepTarget,StepRecoveryError,guarded_step


class StepTests(unittest.TestCase):
    def test_reverse_rollback_preserves_opaque_tokens_and_retries(self):
        x=[0];calls=[];token=object()
        def save():calls.append('save');return token
        def undo(t):self.assertIs(t,token);x[0]=0;calls.append('undo')
        def release(t):self.assertIs(t,token);calls.append('release')
        target=StepTarget('x',x,save,undo,release)
        def bad():x[0]=8;raise ValueError('late')
        with self.assertRaises(StepRecoveryError) as e:guarded_step(bad,[target])
        self.assertTrue(e.exception.original_preserved);self.assertEqual(x,[0])
        self.assertEqual(calls,['save','undo','release'])
        self.assertEqual(guarded_step(lambda:4,[target]),4)

    def test_capture_failure_precedes_operation_and_cleans_backups(self):
        calls=[]
        def bad():raise ValueError('capture')
        a=StepTarget('a',object(),lambda:None,lambda t:None,lambda t:calls.append(t))
        b=StepTarget('b',object(),bad,lambda t:None)
        with self.assertRaises(StepRecoveryError) as e:guarded_step(lambda:calls.append('step'),[a,b])
        self.assertTrue(e.exception.original_preserved);self.assertEqual(calls,[None])

    def test_failed_rollback_never_claims_preservation(self):
        calls=[]
        def bad():raise ValueError('broken')
        a=StepTarget('a',object(),lambda:0,lambda t:calls.append('a'))
        b=StepTarget('b',object(),lambda:0,lambda t:bad())
        with self.assertRaises(StepRecoveryError) as e:guarded_step(bad,[a,b])
        self.assertFalse(e.exception.original_preserved)
        self.assertTrue(e.exception.rollback_errors);self.assertEqual(calls,['a'])

    def test_success_release_error_does_not_claim_rollback(self):
        def bad(t):raise ValueError('release')
        t=StepTarget('x',object(),lambda:None,lambda t:None,bad)
        with self.assertRaises(StepRecoveryError) as e:guarded_step(lambda:7,[t])
        self.assertFalse(e.exception.original_preserved)
        self.assertEqual(e.exception.phase,'release')

    def test_empty_targets_reject_before_operation(self):
        calls=[]
        with self.assertRaises(StepRecoveryError) as e:guarded_step(lambda:calls.append(1),[])
        self.assertEqual(e.exception.phase,'admission');self.assertEqual(calls,[])

    def test_partial_capture_release_gets_none(self):
        calls=[]
        def bad():raise ValueError('capture')
        t=StepTarget('x',object(),bad,lambda token:None,lambda token:calls.append(token))
        with self.assertRaises(StepRecoveryError):guarded_step(lambda:None,[t])
        self.assertEqual(calls,[None])
