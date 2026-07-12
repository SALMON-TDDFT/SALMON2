from pathlib import Path
import os
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
source = (ROOT / "src/gs/dc/lcfo_wannier_sawf_orchestrator.f90").read_text().lower()
flux = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text().lower()

required = [
    "type, public :: t_sawf_environment_receipt",
    "build_sawf_environment_execution_plan",
    "validate_sawf_environment_receipts",
    "representative_fragment",
    "operation_index",
    "generated_independently",
    "requires_execution",
    "same_supercell_fingerprint",
]
for token in required:
    assert token in source, f"missing hierarchical orchestrator contract: {token}"

assert "call build_sawf_environment_execution_plan" in flux

fc = shutil.which("gfortran")
if not fc:
    raise SystemExit("gfortran is required")

driver = r'''program check_orchestrator
  use lcfo_wannier_sawf_orchestrator
  implicit none
  type(t_sawf_environment_receipt),allocatable :: receipt(:)
  logical :: ok
  character(256) :: msg
  call build_sawf_environment_execution_plan([1,1,3],[1,2,1], &
    [.true.,.false.,.true.],'cell-A',receipt,ok,msg)
  call require(ok,'valid plan')
  call require(all(receipt%requires_execution.eqv.[.true.,.false.,.true.]),'execution schedule')
  receipt(1)%completed=.true.;receipt(3)%completed=.true.
  call validate_sawf_environment_receipts(receipt,'cell-A',ok,msg)
  call require(ok,'representative receipts accepted: '//trim(msg))
  receipt(1)%completed=.false.
  call validate_sawf_environment_receipts(receipt,'cell-A',ok,msg)
  call require(.not.ok,'missing representative receipt rejected')
  receipt(1)%completed=.true.
  call validate_sawf_environment_receipts(receipt,'cell-B',ok,msg)
  call require(.not.ok,'cross-supercell receipt rejected')
  write(*,'(a)') 'PASS hierarchical SAWF execution schedule and receipts'
contains
  subroutine require(condition,text)
    logical,intent(in)::condition;character(*),intent(in)::text
    if(.not.condition)then;write(*,'(a)')trim(text);error stop 1;end if
  end subroutine
end program'''

with tempfile.TemporaryDirectory(prefix="sawf-orchestrator-") as td:
    td = Path(td)
    (td / "driver.f90").write_text(driver)
    result = subprocess.run(
        [fc, str(ROOT / "src/gs/dc/lcfo_wannier_sawf_orchestrator.f90"),
         str(td / "driver.f90"), "-o", str(td / "check")],
        cwd=td, env=dict(os.environ), text=True, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if result.returncode:
        raise SystemExit(result.stdout)
    result = subprocess.run([str(td / "check")], cwd=td, text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(result.stdout)
    assert "PASS hierarchical SAWF" in result.stdout

print("PASS production hierarchical SAWF execution-plan contract")
