"""Defaults and explicit history choices in newly generated WSe2 inputs."""
import unittest
import preflight


class PreparationTest(unittest.TestCase):
    def test_omitted_and_explicit_off(self):
        for text in ("&functional xc='hse06' /", "&functional exx_factor_mode='OFF' /",
                     '&functional exx_factor_mode="off" /',
                     "&functional xc='hse06' ! exx_factor_mode='learned'\n/"):
            preflight.check_factor_mode(text)

    def test_explicit_non_off_and_malformed_are_rejected(self):
        for text in ("&functional exx_factor_mode='learned' /",
                     "&functional exx_factor_mode='strict' /",
                     "&functional exx_factor_mode=learned /",
                     "&functional exx_factor_mode='off',exx_factor_mode='learned' /"):
            with self.assertRaises(AssertionError):
                preflight.check_factor_mode(text)


if __name__ == '__main__':
    unittest.main()
