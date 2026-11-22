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
    subprocess.run([str(build/'update_schedule')], check=True)
    subprocess.run([str(build/'frozen_kinetic')], check=True)
    subprocess.run([str(build/'frozen_window')], check=True)
    subprocess.run([str(build/'frozen_response')], check=True)
    subprocess.run([str(build/'frozen_time')], check=True)
    subprocess.run([str(build/'frozen_limit')], check=True)
    subprocess.run([str(build/'coherent_optical')], check=True)
    subprocess.run([str(build/'projected_derivative')], check=True)
    subprocess.run([str(build/'grid_translation')], check=True)
    subprocess.run([str(build/'q2_head')], check=True)
    subprocess.run(['mpiexec','-n','2',str(build/'impulse_resolvent')],check=True)
    subprocess.run(['mpiexec','-n','4',str(build/'equilibrium_q')],check=True)
    for mode in ['window','window-exp','window-empty','window-missing','window-failed','window-key','window-one-reset','hermitian','hermitian-coherent','hermitian-q','hermitian-eta','hermitian-invalid','pairs', 'pairs-q', 'pairs-bad', 'pairs-coherent', 'pairs-prepare', 'pairs-invalid', 'valid', 'empty-pairs', 'cap', 'bad-action', 'coherent-reject', 'nonfinite', 'infinite',
                 'optical', 'optical-empty', 'optical-failure', 'optical-coherent', 'optical-nan', 'optical-eta-reject', 'q2-short', 'q2-conflict']:
        full_run = Path(temporary)/f'frozen-full-{mode}'
        full_run.mkdir()
        subprocess.run(['mpiexec', '-n', '4', str(build/'frozen_full'), mode], cwd=full_run, check=True)
    for mode in ['valid','metal','empty','short','negative','nan','infinite','layout','qfit','geometry']:
        subprocess.run(['mpiexec','-n','4',str(build/'density_gamma'),mode],check=True)
    subprocess.run(['mpiexec','-n','4',str(build/'density_full')],check=True)
    subprocess.run(['mpiexec','-n','4',str(build/'density_full'),'pairs'],check=True)
    subprocess.run(['mpiexec','-n','4',str(build/'density_full'),'hermitian'],check=True)
    frozen_run = Path(temporary)/'frozen-mpi'
    frozen_run.mkdir()
    subprocess.run(['mpiexec', '-n', '4', str(build/'frozen_q')], cwd=frozen_run, check=True)
    for mode in ['frequency','frequency-q0']:
        run=Path(temporary)/('frozen-'+mode)
        run.mkdir()
        subprocess.run(['mpiexec','-n','4',str(build/'frozen_q'),mode],cwd=run,check=True)
    frozen_failure = Path(temporary)/'frozen-remote-failure'
    frozen_failure.mkdir()
    subprocess.run(['mpiexec', '-n', '4', str(build/'frozen_q'), 'bad-action'], cwd=frozen_failure, check=True)
    frozen_output_failure = Path(temporary)/'frozen-output-failure'
    frozen_output_failure.mkdir()
    subprocess.run(['mpiexec', '-n', '4', str(build/'frozen_q'), 'existing-output'],
                   cwd=frozen_output_failure, check=True)
    subprocess.run([str(build/'forms')], check=True)
    subprocess.run([str(build/'q_batch')], check=True)
    subprocess.run([str(build/'projection_capture')], check=True)
    subprocess.run([str(build/'pdep')], check=True)
    subprocess.run([str(build/'response_space')], check=True)
    subprocess.run([str(build/'response_symmetry')], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_response_star.py'),
                    str(build/'response_star'), '4', '12'], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_reduced_q.py'),
                    str(build/'reduced_q'), '4', '12'], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_pdep_q.py'),
                    str(build/'pdep_q'), '4', '12'], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_reduced_w_layout.py'),
                    str(build/'pipeline')], check=True)
    subprocess.run([str(build/'static_sternheimer')], check=True)
    subprocess.run([str(build/'longwave')], check=True)
    subprocess.run([str(build/'gamma_optical')], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_gamma_pdep.py'),
                    str(build/'gamma_pdep'), '4', '12'], check=True)
    subprocess.run([str(build/'dielectric_inverse')], check=True)
    for script, executable in [('test_fortran_screen.py', 'assembly'),
                               ('test_fortran_screen_head.py', 'head')]:
        subprocess.run([sys.executable, str(root/'testsuites/hsex'/script), str(build/executable)], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_distributed_spinor.py'),
                    str(build/'pipeline')], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_fortran_screen_pipeline.py'),
                    str(build/'pipeline'), '5', '1'], check=True)
    subprocess.run([sys.executable, str(root/'testsuites/hsex/test_fortran_scalar_screen.py'),
                    str(build/'pipeline'), '5', '2'], check=True)
