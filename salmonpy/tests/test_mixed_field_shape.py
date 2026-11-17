import unittest
from salmonpy.examples.mpi_mixed_native import _field_layout


class FieldShapeTests(unittest.TestCase):
    def test_default_and_expanded_centered_points(self):
        grid,points=_field_layout((40,40,40))
        self.assertEqual(grid.counts,(40,40,40))
        self.assertEqual(points,{10:(19,20,20),20:(20,20,20)})
        large,points=_field_layout((96,80,64))
        self.assertEqual(large.counts,(96,80,64))
        self.assertEqual(points,{10:(47,40,32),20:(48,40,32)})
        for p in points.values():
            self.assertEqual(tuple(large.coordinates[d][p[d]] for d in range(3)),(-1.,1.,1.) if p[0]==47 else (1.,1.,1.))

    def test_reject_invalid_fixture_shapes(self):
        for value in ((39,40,40),(41,40,40),(40,40),(40,40,True),(40.,40,40)):
            with self.subTest(value=value),self.assertRaises(ValueError):_field_layout(value)

    def test_clock_failure_is_diagnostic_only(self):
        from unittest.mock import patch
        from salmonpy.examples.mpi_mixed_native import _clock
        with patch('salmonpy.examples.mpi_mixed_native.time.perf_counter',side_effect=OSError('clock')):
            self.assertIsNone(_clock())
        with patch('salmonpy.examples.mpi_mixed_native.time.perf_counter',return_value=float('nan')):
            self.assertIsNone(_clock())
