"""Build an audited native library and run independent original-H regressions."""
from pathlib import Path
import os
import platform
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='salmon-original-h-') as temporary:
    build = Path(temporary)/'build'
    subprocess.run([sys.executable, str(root/'tools/hsex/build_q_openmp.py'), str(build)], check=True)
    library = build/('libhsex_q.dylib' if platform.system() == 'Darwin' else 'libhsex_q.so')
    env = dict(os.environ, HSEX_Q_KERNEL=str(library), OPENBLAS_NUM_THREADS='1',
               OMP_NUM_THREADS='12', OMP_THREAD_LIMIT='12', OMP_DYNAMIC='FALSE')
    for case in ('original_h_refinement', 'original_h_diagonal', 'refinement_large'):
        subprocess.run([sys.executable, str(root/'testsuites/hsex'/f'test_{case}.py')],
                       env=env, check=True)
