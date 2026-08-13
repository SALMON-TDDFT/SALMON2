#!/usr/bin/env python3
"""Focused tests for the Si8 MPI memory monitor."""

import importlib.util
from pathlib import Path
import unittest


RUNNER = Path(__file__).with_name("run_si8_overlapping_wannier_memory.py")


def load_runner():
    spec = importlib.util.spec_from_file_location("si8_memory_runner", RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load Si8 memory runner")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MonitorTests(unittest.TestCase):
    def test_vm_stat_available_bytes(self):
        runner = load_runner()
        sample = """Mach Virtual Memory Statistics: (page size of 16384 bytes)\nPages free: 10.\nPages inactive: 20.\nPages speculative: 3.\nPages purgeable: 2.\nPages active: 99.\n"""
        self.assertEqual(runner.parse_available_bytes(sample), 35 * 16384)

    def test_salmon_rss_by_descendant(self):
        runner = load_runner()
        sample = """100 1 1000 mpirun -np 8 salmon\n101 100 2048 salmon\n102 100 4096 salmon\n200 1 999 other\n"""
        self.assertEqual(runner.parse_rank_rss_kib(sample, 100), {101: 2048, 102: 4096})

    def test_latest_phase(self):
        runner = load_runner()
        text = "start\n[OW-GS-DIAGNOSTIC] phase one\nnoise\n[OW-GS-DIAGNOSTIC] phase two\n"
        self.assertEqual(runner.latest_phase(text), "[OW-GS-DIAGNOSTIC] phase two")

    def test_safety_gate(self):
        runner = load_runner()
        self.assertEqual(runner.safety_reason(7 << 30, {1: 100}, 8 << 30, 3 << 20), "available_memory")
        self.assertEqual(runner.safety_reason(20 << 30, {1: 4 << 20}, 8 << 30, 3 << 20), "rank_rss")
        self.assertIsNone(runner.safety_reason(20 << 30, {1: 2 << 20}, 8 << 30, 3 << 20))

    def test_runner_never_cycles_stop_continue(self):
        text = RUNNER.read_text()
        self.assertNotIn("SIGSTOP", text)
        self.assertNotIn("SIGCONT", text)


if __name__ == "__main__":
    unittest.main()
