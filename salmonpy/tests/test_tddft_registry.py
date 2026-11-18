import unittest

class RegistryTests(unittest.TestCase):
    def test_explicit_public_registry_entry(self):
        from salmonpy.embedded_tddft import NativeTDDFTBackend,NativeTDDFTRegistry
        self.assertTrue(callable(NativeTDDFTBackend.open_registry))
        self.assertTrue(callable(NativeTDDFTRegistry.open_model))

    def test_close_attempts_all_models_in_token_order_and_can_retry(self):
        from salmonpy._tddft_registry import NativeTDDFTRegistry
        from salmonpy.errors import RuntimeFailure
        from types import SimpleNamespace
        events=[]
        native=SimpleNamespace(closed=False,_live=lambda:None)
        def detach():native.closed=True;events.append('detach')
        native.close=detach
        registry=NativeTDDFTRegistry(None,native,None)
        class Owner:
            def __init__(self,token,fail=False):self.model_token=token;self.fail=fail
            def close(self):
                events.append(self.model_token)
                if self.fail:self.fail=False;raise RuntimeFailure('injected close failure')
                registry._models.discard(self)
        first,second=Owner(-5,True),Owner(-2)
        registry._models={second,first}
        with self.assertRaisesRegex(RuntimeFailure,'injected close failure'):registry.close()
        self.assertEqual(events,[-5,-2])
        self.assertFalse(native.closed)
        registry.close()
        self.assertEqual(events,[-5,-2,-5,'detach'])
        self.assertTrue(native.closed)

    def test_mixed_python_model_selection_rejects_before_queued_state_is_consumed(self):
        from salmonpy._tddft_registry import _RegisteredModel
        from salmonpy.errors import InputError
        from types import SimpleNamespace
        from pathlib import Path
        from unittest.mock import Mock
        comm=SimpleNamespace(allgather=lambda value:[value,-99] if isinstance(value,int) else [value,value])
        native=SimpleNamespace(context=SimpleNamespace(solver_comm=comm),_live=lambda:None)
        model=object.__new__(_RegisteredModel)
        model.registry=SimpleNamespace(closed=False,_native=SimpleNamespace(closed=False))
        model.model_token=-1;model._native=native;model.closed=False;model.failed=False
        model._serial=0;model._initialized=True;model.work=Path('.')
        model.lib=SimpleNamespace(sp_tddft_model_new=Mock(side_effect=AssertionError('native spawn must not be called')))
        with self.assertRaisesRegex(InputError,'model selections differ'):model.new_state()
        model.lib.sp_tddft_model_new.assert_not_called()
