"""Run installed-wheel tests; reject imports from the source checkout.

--suite full runs all tests. --suite core checks public imports and the native
LD/SBE contracts on every supported wheel/interpreter. No optional install,
compiler invocation or solver job is performed by this runner itself.
"""
import argparse
from pathlib import Path
import unittest

CORE_PATTERNS=('test_public_api.py','test_ld.py','test_sbe.py','test_manual.py','test_quickstart.py')


def build_suite(root,profile):
    tests=Path(root)/'salmonpy/tests'
    loader=unittest.TestLoader()
    if profile=='full':
        suite=loader.discover(str(tests))
        if not suite.countTestCases():raise RuntimeError('source distribution has no installed-wheel tests')
        return suite
    suite=unittest.TestSuite()
    for pattern in CORE_PATTERNS:
        selected=loader.discover(str(tests),pattern=pattern)
        if not selected.countTestCases():raise RuntimeError('source distribution is missing '+pattern)
        suite.addTests(selected)
    return suite


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root',type=Path)
    parser.add_argument('--suite',choices=('full','core'),default='full')
    args=parser.parse_args()
    import salmonpy
    root=args.root.resolve()
    loaded=Path(salmonpy.__file__).resolve()
    if loaded.is_relative_to(root):
        raise RuntimeError(f'Test imported checkout instead of wheel: {loaded}')
    print(f'Testing installed package: {loaded}; suite={args.suite}',flush=True)
    result=unittest.TextTestRunner(verbosity=2).run(build_suite(root,args.suite))
    raise SystemExit(not result.wasSuccessful())


if __name__=='__main__':main()
