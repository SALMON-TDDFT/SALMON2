import unittest
from types import SimpleNamespace
from salmonpy.errors import InputError


class SymmetryHistoryProtocolTests(unittest.TestCase):
    def test_legacy_exact_taylor_does_not_require_new_capability(self):
        from salmonpy.fast_hse_mpi import _require_history_protocol
        _require_history_protocol(SimpleNamespace(),'off','hse_taylor4')
        _require_history_protocol(SimpleNamespace(),'off','hse_taylor4_full')

    def test_extended_modes_require_matching_core(self):
        from salmonpy.fast_hse_mpi import _require_history_protocol
        class Query:
            def __call__(self):return 1
        for factor,propagator in [('off','hse_ptcn'),('learned','hse_taylor4'),('strict','hse_taylor4'),('learned','hse_ptcn')]:
            with self.subTest(factor=factor,propagator=propagator):
                with self.assertRaisesRegex(InputError,'history'):
                    _require_history_protocol(SimpleNamespace(),factor,propagator)
                _require_history_protocol(SimpleNamespace(sp_fast_hse_symmetry_history_protocol=Query()),factor,propagator)

    def test_cn_uses_existing_native_interval_contract(self):
        from salmonpy.fast_hse_mpi import _symmetry_factor
        for body in ["exx_factor_mode='strict'", "exx_factor_mode='learned'", "exx_factor_mode='learned',exx_factor_exact_interval=2"]:
            with self.assertRaisesRegex(InputError,'PT-CN'):_symmetry_factor({'functional':body},'hse_ptcn')
        self.assertEqual(_symmetry_factor({'functional':"exx_factor_mode='learned',exx_factor_exact_interval=1"},'hse_ptcn'),'learned')
        self.assertEqual(_symmetry_factor({'functional':"exx_factor_mode='learned',exx_factor_exact_interval=2"},'hse_taylor4'),'learned')
