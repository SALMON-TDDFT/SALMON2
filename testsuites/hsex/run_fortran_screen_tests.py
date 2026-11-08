"""Build audited production screening modules and compare with Python references."""
from pathlib import Path
import os
import platform
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory(prefix='salmon-fortran-screen-') as temporary:
    build = Path(temporary) / 'build'
    env = dict(os.environ)
    if platform.system() != 'Darwin':
        env.setdefault('SCREEN_NATIVE_FLAGS', '')
    subprocess.run(['bash', str(root/'testsuites/hsex/build_fortran_screen.sh'), str(build)],
                   env=env, check=True)
    subprocess.run([str(build/'static_sternheimer')], check=True)
    for script, executable in [('test_fortran_screen.py', 'assembly'),
                               ('test_fortran_screen_head.py', 'head')]:
        subprocess.run([sys.executable, str(root/'testsuites/hsex'/script), str(build/executable)], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_distributed_spinor.py'),
                    str(build/'pipeline')], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_fortran_screen_pipeline.py'),
                    str(build/'pipeline'), '1', '1'], check=True)
