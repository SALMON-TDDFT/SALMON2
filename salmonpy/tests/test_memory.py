import unittest
import warnings
from unittest.mock import patch
from salmonpy.memory import MemoryEstimate, MemoryBudgetWarning, memory_advice, process_memory


class MemoryAdviceTests(unittest.TestCase):
    def test_candidate_coexistence_and_reserve_counted_without_stop(self):
        with warnings.catch_warnings(record=True) as messages:
            warnings.simplefilter('always')
            report=memory_advice(MemoryEstimate(100,200,30),budget_bytes=300)
        self.assertEqual(report.estimated_peak_bytes,330)
        self.assertTrue(report.exceeds_budget)
        self.assertEqual(len(messages),1)
        self.assertIs(messages[0].category,MemoryBudgetWarning)
        self.assertEqual(report.to_dict()['budget_bytes'],300)

    def test_unknown_is_not_reported_as_zero_or_safe(self):
        report=memory_advice(MemoryEstimate(None,200,30),budget_bytes=300,warn=False)
        self.assertIsNone(report.estimated_peak_bytes)
        self.assertEqual(report.known_peak_bytes,230)
        self.assertTrue(report.unknown)
        self.assertIsNone(report.exceeds_budget)

    def test_warning_filter_error_does_not_become_computation_failure(self):
        with warnings.catch_warnings():
            warnings.simplefilter('error',MemoryBudgetWarning)
            report=memory_advice(MemoryEstimate(100,200),budget_bytes=1)
        self.assertTrue(report.exceeds_budget)
        self.assertTrue(report.warning_error)

    def test_probe_unavailable_is_reported_without_error(self):
        with patch('salmonpy.memory._peak_rss',side_effect=OSError('disabled')):
            value,basis=process_memory()
        self.assertIsNone(value)
        self.assertIn('disabled',basis)

    def test_restore_counts_native_and_snapshot_copies_without_loading(self):
        import tempfile
        from pathlib import Path
        from types import SimpleNamespace
        from salmonpy._memory_advice import restore_advice,RESERVE
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name,size in [('field-0.bin',100),('electronic-0.bin',200),('rank-0.json',4),('frames.json',8)]:
                (root/name).write_bytes(b'x'*size)
            raw=SimpleNamespace(_initial=SimpleNamespace(eigen=SimpleNamespace(shape=(2,1))))
            coupled=SimpleNamespace(field=object(),response=SimpleNamespace(states={1:SimpleNamespace(native=True),2:SimpleNamespace(_raw=raw)}),
                memory_reports=[],memory_budget_bytes=4*2**30)
            world=SimpleNamespace(rank=0,allgather=lambda x:[x])
            with patch('salmonpy._memory_advice.process_memory',return_value=(1024,'test RSS')):
                with patch.object(Path,'read_bytes',side_effect=AssertionError('estimate must not load payload')):
                    report=restore_advice(world,coupled,root,{})
        self.assertEqual(report.ranks[0].candidate_bytes,4*300+32*12+8*(16*4+48))
        self.assertEqual(report.known_peak_bytes,1024+RESERVE+report.ranks[0].candidate_bytes)

    def test_report_write_failure_is_only_diagnostic(self):
        from types import SimpleNamespace
        from pathlib import Path
        from salmonpy._memory_advice import _record
        records=[]
        world=SimpleNamespace(rank=0,allgather=lambda x:[x])
        with patch.object(Path,'write_text',side_effect=OSError('report disk unavailable')):
            report=_record(world,MemoryEstimate(10),100,'test',records,Path('unused'))
        self.assertFalse(report.exceeds_budget)
        self.assertIn('report disk unavailable',records[0]['report_write_errors'][0])

    def test_child_worker_unknown_survives_restore_estimate(self):
        import tempfile
        from pathlib import Path
        from types import SimpleNamespace
        from salmonpy._memory_advice import restore_advice
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            for name in ('worker-0.bin','rank-0.json','frames.json'):
                (root/name).write_bytes(b'x'*8)
            state=SimpleNamespace(_process=object(),checkpoint_protocol='salmonpy.direct-state/1')
            coupled=SimpleNamespace(field=None,response=SimpleNamespace(states={1:state}),memory_reports=[],memory_budget_bytes=4*2**30)
            world=SimpleNamespace(rank=0,allgather=lambda x:[x])
            with patch('salmonpy._memory_advice.process_memory',return_value=(100,'parent RSS')):
                report=restore_advice(world,coupled,root,{})
        self.assertIsNone(report.estimated_peak_bytes)
        self.assertTrue(any('child worker' in reason for reason in report.unknown))
