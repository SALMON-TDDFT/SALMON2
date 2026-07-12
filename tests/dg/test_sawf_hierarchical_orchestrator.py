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
    "complete_sawf_seed_bundle",
]
for token in required:
    assert token in source, f"missing hierarchical orchestrator contract: {token}"

assert "call build_sawf_environment_execution_plan" in flux
assert "call build_sawf_seed_bundles" in flux
assert "call atomic_create_directory(sawf_seed_bundles" in flux
assert "call select_sawf_local_complete_shells" in flux
assert "call solve_sawf_local_generalized_eigensystem" in flux
assert "call write_sawf_representative_local_seed" in flux
assert "call run_sawf_local_preprocessing" in flux
assert "call read_sawf_nnkp_neighbors" in flux
assert "call restrict_sawf_stabilizer_representation" in flux
assert "call build_sawf_local_band_representation" in flux
assert flux.index("call write_sawf_representative_local_seed") < flux.index(
    "call run_sawf_local_preprocessing"
) < flux.index("call read_sawf_nnkp_neighbors")
assert "preprocess_only" in flux
assert "esp_tot(1:nband_wann" not in flux.split(
    "subroutine write_sawf_representative_local_seed", 1
)[1].split("end subroutine write_sawf_representative_local_seed", 1)[0]

fc = shutil.which("gfortran")
if not fc:
    raise SystemExit("gfortran is required")

driver = r'''program check_orchestrator
  use lcfo_wannier_sawf_orchestrator
  implicit none
  type(t_sawf_environment_receipt),allocatable :: receipt(:)
  type(t_sawf_seed_bundle),allocatable :: bundle(:)
  integer,allocatable :: stabilizer(:)
  logical :: ok
  character(256) :: msg
  call build_sawf_environment_execution_plan([1,1,3],[1,2,1], &
    [.true.,.false.,.true.],'cell-A',receipt,ok,msg)
  call require(ok,'valid plan')
  call require(all(receipt%requires_execution.eqv.[.true.,.false.,.true.]),'execution schedule')
  call build_sawf_seed_bundles(receipt,'.','seed',bundle,ok,msg)
  call require(ok.and.size(bundle)==2,'representative seed bundles')
  call require(bundle(1)%environment==1.and.bundle(2)%environment==3,'bundle environments')
  call require(index(bundle(1)%directory,'environment-000001')>0,'isolated seed directory')
  call require(trim(bundle(1)%seedname)=='seed-env-000001','unique seed name')
  call select_sawf_environment_stabilizer(1,reshape([1,2,3,1,3,2],[3,2]), &
    [1,1,2,2],[1,2,1,2],[1,2,2,1],stabilizer,ok,msg)
  call require(ok.and.all(stabilizer==[1,2]),'actual-group local stabilizer')
  call complete_sawf_seed_bundle(bundle(1),receipt(1),ok,msg)
  call require(ok.and.receipt(1)%completed,'first artifact bundle completed: '//trim(msg))
  call require(receipt(1)%num_bands==2.and.receipt(1)%num_wann==2,'ragged local dimensions')
  bundle(2)%directory='missing-environment'
  call complete_sawf_seed_bundle(bundle(2),receipt(3),ok,msg)
  call require(.not.ok.and..not.receipt(3)%completed,'missing artifact rejected')
  bundle(2)%directory='./environment-000003'
  call complete_sawf_seed_bundle(bundle(2),receipt(3),ok,msg)
  call require(ok.and.receipt(3)%completed,'second artifact bundle completed: '//trim(msg))
  call propagate_sawf_representative_receipts(receipt,ok,msg)
  call require(ok.and.receipt(2)%completed,'replica receipt propagated: '//trim(msg))
  call require(receipt(2)%num_bands==receipt(1)%num_bands.and. &
    receipt(2)%num_wann==receipt(1)%num_wann,'replica ragged dimensions')
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
    for environment in (1, 3):
        directory = td / f"environment-{environment:06d}"
        directory.mkdir()
        seed = f"seed-env-{environment:06d}"
        for suffix in ("win", "dmn", "chk"):
            (directory / f"{seed}.{suffix}").write_text(f"{suffix}\n")
        (directory / f"{seed}.eig").write_text("1 1 -1.0\n2 1 2.0\n")
        (directory / f"{seed}.mmn").write_text("local overlaps\n2 1 3\n")
        (directory / f"{seed}.amn").write_text("local projections\n 2 1 2\n")
        (directory / f"{seed}.sawf-fingerprint").write_text("cell-A\n")
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
