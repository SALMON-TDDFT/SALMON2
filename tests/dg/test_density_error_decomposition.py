"""Tests for density differences; no SCF or Wannier regeneration needed."""
import importlib.util
from pathlib import Path
import unittest
import tempfile


class DensityDecompositionTest(unittest.TestCase):
    def test_snapshot_reader(self):
        path = Path(__file__).with_name("density_error_decomposition.py")
        spec = importlib.util.spec_from_file_location("decomposition", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertTrue(hasattr(module, "read_snapshot_set"), "production snapshot reader missing")
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / "density"
            files = [Path(str(prefix) + f".rank-{rank:08d}") for rank in range(2)]
            def shard(rank, point):
                return (f"SALMON_DG_DENSITY_DIAGNOSTIC_V1\n2 {rank} {rank+1} 1 2\n"
                        f"11 22 33\n4 0.001\n{point} 1 2 2 2 2 0.5\nEND_DENSITY_DIAGNOSTIC\n")
            for rank, file in enumerate(files):
                file.write_text(shard(rank, rank+1))
            result = module.read_snapshot_set(prefix)
            self.assertEqual(result["metrics"]["dg_increment"], 0)
            self.assertEqual(result["electron_counts"]["conventional"], 4)
            self.assertEqual(len(result["sha256"]), 2)
            for malformed in (shard(1, 1), shard(1, 2).replace("11 22 33", "11 23 33"),
                              shard(1, 2).replace("END_DENSITY_DIAGNOSTIC", ""),
                              shard(1, 2).replace("2 1 2 1 2", "2 1 1 1 2"),
                              shard(1, 2).replace("0.5", "nan")):
                files[1].write_text(malformed)
                with self.assertRaises(ValueError):
                    module.read_snapshot_set(prefix)

    def test_decomposition(self):
        path = Path(__file__).with_name("density_error_decomposition.py")
        self.assertTrue(path.exists(), "density decomposition implementation missing")
        spec = importlib.util.spec_from_file_location("decomposition", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        compare = module.compare_densities
        # Equal baseline and total norms must not hide opposite spatial shifts.
        result = compare([2, 2], [3, 1], [1, 3], [1, 1], relaxed=[2, 2])
        self.assertAlmostEqual(result["dc_lcfo_correction"], 0.5)
        self.assertAlmostEqual(result["dg_increment"], 1.0)
        self.assertAlmostEqual(result["frozen_total_displacement"], 0.5)
        self.assertAlmostEqual(result["relaxation_displacement"], 0.5)
        self.assertAlmostEqual(result["cross_term"], -1.0)
        self.assertFalse(result["absolute_dc_error_available"])
        # Weighted grid refinement leaves every diagnostic unchanged.
        split = compare([2]*4, [3, 3, 1, 1], [1, 1, 3, 3], [0.5]*4)
        for key in ("dc_lcfo_correction", "dg_increment", "cross_term"):
            self.assertAlmostEqual(result[key], split[key])
        same = compare([2, 2], [2, 2], [2, 2], [1, 1])
        self.assertEqual(same["dg_increment"], 0)
        self.assertIsNone(same["relaxation_displacement"])
        for args in (([], [], [], []), ([0], [1], [1], [1]),
                     ([1], [1, 2], [1], [1]), ([1], [1], [1], [-1]),
                     ([1], [float("nan")], [1], [1])):
            with self.assertRaises(ValueError):
                compare(*args)


if __name__ == "__main__":
    unittest.main()
