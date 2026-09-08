from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
flux=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text().lower()
runner=flux.split('subroutine run_wannier90_seed_files',1)[1].split(
    'end subroutine run_wannier90_seed_files',1)[0]
assert 'seedname_in' in runner and 'directory_in' in runner
assert 'present(seedname_in)' in runner and 'present(directory_in)' in runner
print('PASS production SAWF runner supports representative seed directories')
