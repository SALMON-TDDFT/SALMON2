"""Bounded native HSE laser checks using an existing compatible Si GS."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'hse_native'))
from check_input_smoke import make_input

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('exe', type=Path)
p.add_argument('--gs', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
p.add_argument('--case', default='ace', choices=['ace', 'full', 'half', 'quarter', 'zero', 'probe'])
a = p.parse_args()
d = a.work.resolve() / a.case
d.mkdir(parents=True, exist_ok=False)
(d / 'restart').symlink_to(a.gs.resolve(), target_is_directory=True)
pseudo = Path(__file__).resolve().parents[2] / 'testsuites/pseudo/Si_rps.dat'
nt, dt = {'ace': (16, .08), 'full': (16, .08), 'half': (32, .04),
          'quarter': (64, .02), 'zero': (16, .08), 'probe': (16, .08)}[a.case]
s = make_input(Path('restart'), pseudo, 4)
s = s.replace("theory='tddft_response'", "theory='tddft_pulse'")
s = s.replace('nt=1', f'nt={nt}').replace('dt=0.16', f'dt={dt}')
s = s.replace('checkpoint_interval=1', f'checkpoint_interval={nt}')
# All pulses use atomic units, z polarization and the same endpoint.
s = s.replace("ae_shape1='impulse'", "ae_shape1='Acos2'\n tw1=1.28\n omega1=.2\n I_wcm2_1=" + ('0' if a.case == 'zero' else '1e13'))
s = s.replace('e_impulse=0.0001', 'e_impulse=0')
if a.case == 'probe':
    s = s.replace('e_impulse=0', 'e_impulse=1e-4\n ae_shape2="impulse"\n epdir_re2=0,0,1\n T1_T2=.16')
if a.case == 'full':
    s += "\n&propagation\n propagator='hse_taylor4_full'\n/\n"
(d / 'inputfile').write_text(s)
with (d / 'run.log').open('w') as out:
    r = subprocess.run(['mpiexec', '-n', '4', str(a.exe.resolve())], input=s, text=True,
                       cwd=d, stdout=out, stderr=subprocess.STDOUT,
                       env=dict(os.environ, OMP_NUM_THREADS='1'), timeout=300)
text = (d / 'run.log').read_text()
result = dict(case=a.case, exit_code=r.returncode, completed=r.returncode == 0 and 'end SALMON' in text)
(d / 'status.json').write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps(result))
assert result['completed'], text[-2500:]
