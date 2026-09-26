"""One-step actual-solver diagnostic export smoke test; arguments match check.py."""
import os
import runpy
import sys
from pathlib import Path
os.environ['SALMON_HSE_SOLVER_DIAGNOSTIC'] = '1'
runpy.run_path(str(Path(__file__).with_name('check.py')), run_name='__main__')
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'samples/dc_hse'))
from locality import read_eigen_pair, eigen_diagnostics
for frag in range(1, 9):
    root = Path(sys.argv[2]) / f'data_dcdft/fragments/{frag:06d}'
    for tag in ('before_subspace','after_subspace','after_cg','after_orthogonalization'):
        path = root / f'hse_eigen_{tag}.bin'
        assert path.exists(), f'missing {path}'
        data = read_eigen_pair(path)
        assert data[0].shape == (8192, 64), data[0].shape
        report = eigen_diagnostics(*data)
        print(frag, tag, report['occupied_rms_residual_Ha'], report['max_orthogonality_error'])
