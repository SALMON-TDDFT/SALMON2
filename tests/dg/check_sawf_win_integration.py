#!/usr/bin/env python3
import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def run(command, *, cwd=None, expect_success=True):
    result = subprocess.run(
        [str(item) for item in command],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    if expect_success and result.returncode:
        raise SystemExit("command failed:\n" + " ".join(map(str, command)) + "\n" + result.stdout)
    return result


parser = argparse.ArgumentParser()
parser.add_argument("--build-dir", default="build-mpi-eigenexa-wannier-lib")
args = parser.parse_args()
build_dir = (ROOT / args.build_dir).resolve()

with tempfile.TemporaryDirectory(prefix="sawf-win-") as temporary:
    temporary = Path(temporary)
    source_dir = temporary / "src"
    binary_dir = temporary / "build"
    run_dir = temporary / "run"
    source_dir.mkdir()
    run_dir.mkdir()

    edge_original = (
        b"\tSiTe_SyMmEtRy \t= false\r\n"
        + b"X" * 5003
        + b"   \r\n"
        + b"! site_symmetry = false  \r\n"
        + b"# symmetrize_eps = 9d-9\r\n"
        + b"site_symmetry_extra = keep\r\n"
        + b"symmetrize_epsilon = keep-too\r\n"
        + b"\tSYMMETRIZE_EPS\t= 4d-4\r\n"
        + b"tail-with-spaces   "
    )
    (run_dir / "edge.win").write_bytes(edge_original)
    (run_dir / "stale.win").write_bytes(
        b"base\nsite_symmetry = true\nsymmetrize_eps = 1d-6\n"
    )
    (run_dir / "stale.chk").write_bytes(b"same-dimension stale checkpoint\n3 3\n")

    (source_dir / "driver.f90").write_text(
        r'''program check_sawf_win
  use lcfo_wannier_sawf_win, only: activate_sawf_win,deactivate_sawf_win, &
    t_atomic_win_writer,begin_atomic_win,finish_atomic_win,abort_atomic_win
  use lcfo_wannier_command, only: select_wannier90_command,is_wannier90_reuse_command
  implicit none
  logical :: ok
  logical :: stale_exists
  type(t_atomic_win_writer) :: writer
  integer :: io_status,unit
  character(512) :: message
  character(1024) :: selected

  call select_wannier90_command('namelist','environment','compiled',selected)
  call require(trim(selected)=='namelist','namelist command precedence')
  call select_wannier90_command('','environment','compiled',selected)
  call require(trim(selected)=='environment','environment command precedence')
  call select_wannier90_command('','','compiled',selected)
  call require(trim(selected)=='compiled','compiled command fallback')
  call require(is_wannier90_reuse_command('REUSE'),'reuse classification')
  call require(.not.is_wannier90_reuse_command('export_only'),'export is not reuse')
  inquire(file='stale.chk',exist=stale_exists)
  call require(stale_exists .and. is_wannier90_reuse_command('reuse'), &
    'stale same-dimension checkpoint reuse gate')

  call write_base('seed.win')
  call activate_sawf_win('seed.win',1.25d-7,ok,message,temp_nonce=7001)
  call require(ok,'first activation: '//trim(message))
  call check_file('seed.win',1.25d-7)
  call activate_sawf_win('seed.win',1.25d-7,ok,message,temp_nonce=7002)
  call require(ok,'repeat activation: '//trim(message))
  call check_file('seed.win',1.25d-7)

  open(20,file='collision.win.tmp.7003.0',status='replace')
  write(20,'(a)') 'SENTINEL'
  close(20)
  call write_base('collision.win')
  call activate_sawf_win('collision.win',2.5d-8,ok,message,temp_nonce=7003)
  call require(ok,'collision retry: '//trim(message))
  call check_file('collision.win',2.5d-8)
  call require(first_line('collision.win.tmp.7003.0')=='SENTINEL','collision file changed')

  call write_base('preserved.win')
  call activate_sawf_win('missing/seed.win',1d-8,ok,message,temp_nonce=7004)
  call require(.not.ok,'missing source must fail')
  call require(count_keyword('preserved.win','site_symmetry')==0,'unrelated base file changed')

  open(21,file='directory.win',status='replace'); write(21,'(a)') 'BASE'; close(21)
  call execute_command_line('rm -f directory.win && mkdir directory.win')
  call activate_sawf_win('directory.win',1d-8,ok,message,temp_nonce=7005)
  call require(.not.ok,'directory source must fail without partial activation')

  call activate_sawf_win('edge.win',3.75d-9,ok,message,temp_nonce=7006)
  call require(ok,'byte-preserving edge activation: '//trim(message))
  call execute_command_line('cp edge.win edge.active.win')
  call deactivate_sawf_win('edge.win',ok,message,temp_nonce=7007)
  call require(ok,'byte-preserving edge deactivation: '//trim(message))
  call deactivate_sawf_win('stale.win',ok,message,temp_nonce=7008)
  call require(ok,'stale enabled win deactivation: '//trim(message))
  call deactivate_sawf_win('absent.win',ok,message,temp_nonce=7009)
  call require(ok,'absent win deactivation must succeed')

  call begin_atomic_win(writer,'atomic.win',unit,ok,message,temp_nonce=7010)
  call require(ok,'atomic base begin: '//trim(message))
  write(unit,'(a)',iostat=io_status) 'BASE WITH TRAILING   '
  call require(io_status==0,'atomic base write')
  call finish_atomic_win(writer,ok,message)
  call require(ok,'atomic base finish: '//trim(message))

  call execute_command_line('mkdir publish.win')
  call begin_atomic_win(writer,'publish.win',unit,ok,message,temp_nonce=7011)
  call require(ok,'publish failure begin')
  write(unit,'(a)',iostat=io_status) 'MUST NOT REPLACE DIRECTORY'
  call finish_atomic_win(writer,ok,message)
  call require(.not.ok,'atomic publish onto directory must fail')
  call abort_atomic_win(writer)
  call begin_atomic_win(writer,'missing/base.win',unit,ok,message,temp_nonce=7012)
  call require(.not.ok,'atomic base open in missing directory must fail')

  write(*,'(a)') 'PASS SAWF win byte-preserving activation and atomic base publication'
contains
  subroutine write_base(path)
    character(*), intent(in) :: path
    open(10,file=path,status='replace')
    write(10,'(a)') 'num_bands = 3'
    write(10,'(a)') 'num_wann = 3'
    write(10,'(a)') 'gamma_only = true'
    write(10,'(a)') 'mp_grid = 1 1 1'
    write(10,'(a)') 'begin projections'
    write(10,'(a)') 'random'
    write(10,'(a)') 'end projections'
    close(10)
  end subroutine

  subroutine check_file(path,tolerance)
    character(*), intent(in) :: path
    real(8), intent(in) :: tolerance
    character(128) :: expected
    write(expected,'(a,1x,es23.15)') 'symmetrize_eps =',tolerance
    call require(count_keyword(path,'site_symmetry')==1,'site_symmetry count')
    call require(count_keyword(path,'symmetrize_eps')==1,'symmetrize_eps count')
    call require(has_exact_line(path,'site_symmetry = true'),'site_symmetry value')
    call require(has_exact_line(path,trim(expected)),'symmetrize_eps value')
  end subroutine

  integer function count_keyword(path,keyword) result(count)
    character(*), intent(in) :: path,keyword
    character(512) :: line
    integer :: io
    count=0
    open(11,file=path,status='old')
    do
      read(11,'(a)',iostat=io) line
      if(io/=0) exit
      if(index(adjustl(line),trim(keyword))==1) count=count+1
    end do
    close(11)
  end function

  logical function has_exact_line(path,expected) result(found)
    character(*), intent(in) :: path,expected
    character(512) :: line
    integer :: io
    found=.false.
    open(12,file=path,status='old')
    do
      read(12,'(a)',iostat=io) line
      if(io/=0) exit
      if(trim(line)==trim(expected)) found=.true.
    end do
    close(12)
  end function

  function first_line(path) result(line)
    character(*), intent(in) :: path
    character(512) :: line
    open(13,file=path,status='old'); read(13,'(a)') line; close(13)
    line=trim(line)
  end function

  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if(.not.condition) then
      write(*,'(a)') 'FAIL '//trim(label)
      error stop 1
    end if
  end subroutine
end program
'''
    )
    (source_dir / "CMakeLists.txt").write_text(
        f'''cmake_minimum_required(VERSION 3.18)
project(sawf_win LANGUAGES Fortran)
add_library(sawf_win
  "{ROOT / 'src/gs/dc/lcfo_wannier_sawf_win.f90'}"
  "{ROOT / 'src/gs/dc/lcfo_wannier_command.f90'}")
add_executable(check_sawf_win driver.f90)
target_link_libraries(check_sawf_win PRIVATE sawf_win)
'''
    )
    run(["cmake", "-S", source_dir, "-B", binary_dir])
    run(["cmake", "--build", binary_dir, "-j", "2"])
    result = run([binary_dir / "check_sawf_win"], cwd=run_dir)
    print(result.stdout.strip())

    preserved = (
        b"X" * 5003
        + b"   \r\n"
        + b"! site_symmetry = false  \r\n"
        + b"# symmetrize_eps = 9d-9\r\n"
        + b"site_symmetry_extra = keep\r\n"
        + b"symmetrize_epsilon = keep-too\r\n"
        + b"tail-with-spaces   "
    )
    active = (run_dir / "edge.active.win").read_bytes()
    canonical = b"\r\nsite_symmetry = true\r\nsymmetrize_eps =   3.750000000000000E-09\r\n"
    if active != preserved + canonical:
        raise SystemExit("SAWF activation did not preserve CRLF/long/trailing/comment bytes")
    deactivated = (run_dir / "edge.win").read_bytes()
    if deactivated != preserved + b"\r\n":
        raise SystemExit(
            f"SAWF deactivation did not preserve all non-key bytes: {deactivated[-120:]!r}"
        )
    if (run_dir / "stale.win").read_bytes() != b"base\n":
        raise SystemExit("stale enabled win was not deactivated exactly")
    if (run_dir / "atomic.win").read_bytes() != b"BASE WITH TRAILING   \n":
        raise SystemExit("atomic base writer changed formatted output")
    if not (run_dir / "publish.win").is_dir():
        raise SystemExit("failed atomic publish replaced its target directory")
    if list(run_dir.glob("*.tmp.7011.*")):
        raise SystemExit("failed atomic publish left owned temporary files")

    mpifort = shutil.which("mpifort")
    mpiexec = shutil.which("mpiexec")
    if not mpifort or not mpiexec:
        raise SystemExit("Task7 collective spawn test requires mpifort and mpiexec")
    mpi_driver = source_dir / "mpi_spawn_failure.f90"
    mpi_driver.write_text(
        r'''program mpi_spawn_failure
  use mpi
  use lcfo_wannier_command, only: execute_wannier90_command
  implicit none
  integer :: ierr,rank,failure
  logical :: ok
  character(512) :: message,filename
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  failure=0; message=''
  if(rank==0) then
    call execute_wannier90_command("sh -c 'exit 7'",ok,message)
    failure=merge(0,1,ok)
  end if
  call MPI_Bcast(failure,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(message,len(message),MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  if(failure/=1 .or. index(message,'exit status=7')==0) call MPI_Abort(MPI_COMM_WORLD,2,ierr)
  write(filename,'(a,i0,a)') 'spawn-rank-',rank,'.done'
  open(20,file=trim(filename),status='replace'); write(20,'(a)') trim(message); close(20)
  call MPI_Barrier(MPI_COMM_WORLD,ierr)
  call MPI_Finalize(ierr)
end program
'''
    )
    mpi_exe = run_dir / "mpi_spawn_failure"
    run(
        [
            mpifort,
            ROOT / "src/gs/dc/lcfo_wannier_command.f90",
            mpi_driver,
            "-o",
            mpi_exe,
        ],
        cwd=run_dir,
    )
    mpi_result = run([mpiexec, "-n", "2", mpi_exe], cwd=run_dir)
    markers = sorted(run_dir.glob("spawn-rank-*.done"))
    if len(markers) != 2:
        raise SystemExit("collective spawn failure did not release both ranks:\n" + mpi_result.stdout)
    print("PASS 2-rank root command failure broadcast without hang")

source = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text()
seed = source.split("subroutine write_wannier_seed_files()", 1)[1].split(
    "end subroutine write_wannier_seed_files", 1
)[0]
validation = source.split("subroutine validate_sawf_projection_configuration()", 1)[1].split(
    "end subroutine validate_sawf_projection_configuration", 1
)[0]

assert seed.index("call deactivate_sawf_win_collective") < seed.index("if(nspin /= 1)")
assert seed.index("call deactivate_sawf_win_collective") < seed.index(
    "call validate_sawf_projection_configuration"
)
assert seed.index("call validate_sawf_projection_configuration") < seed.index(
    "call write_wannier_base_win_atomic"
)
base_writer = source.split("subroutine write_wannier_base_win_atomic", 1)[1].split(
    "end subroutine write_wannier_base_win_atomic", 1
)[0]
assert "call begin_atomic_win" in base_writer
assert "call finish_atomic_win" in base_writer
assert "status='replace'" not in base_writer
assert seed.index("call write_wannier_amn_file") < seed.index("call write_wannier_mmn_file")
assert seed.index("call write_wannier_mmn_file") < seed.index("call generate_sawf_dmn")
assert seed.index("call generate_sawf_dmn") < seed.index("call activate_sawf_win_collective")
assert seed.index("call activate_sawf_win_collective") < seed.index("call run_wannier90_seed_files")
nonoff = seed.split("if(trim(wannier_site_symmetry) /= 'off') then", 1)[1].split("end if", 1)[0]
assert nonoff.index("call generate_sawf_dmn") < nonoff.index("call activate_sawf_win_collective")
assert "call run_wannier90_seed_files" not in nonoff
assert seed.count("call generate_sawf_dmn") == 1
assert seed.count("call activate_sawf_win_collective") == 1
assert "wannier_dis_froz_max" in validation
assert "frozen window" in validation.lower()
assert "wannier_dis_win_max" not in validation
assert "call comm_bcast" in source.split("subroutine activate_sawf_win_collective", 1)[1].split(
    "end subroutine activate_sawf_win_collective", 1
)[0]
assert "trim(wannier_site_symmetry) /= 'off'" in seed
assert seed.index("call get_cached_wannier90_command") < seed.index("call reject_sawf_reuse_collective")
assert seed.index("call reject_sawf_reuse_collective") < seed.index("call write_wannier_base_win_atomic")
assert "call run_wannier90_seed_files(resolved_wannier_command)" in seed
assert "is_wannier90_export_only_requested()" not in seed
export_block = seed.split(
    "if(is_wannier90_export_only_command(resolved_wannier_command)) then", 1
)[1].split("return", 1)[0]
assert "trim(wannier_site_symmetry) /= 'off'" in export_block
assert "wannier90_command='import_only'" in export_block
assert "SALMON_WANNIER90_COMMAND=reuse" in export_block
assert export_block.index("wannier90_command='import_only'") < export_block.index(
    "SALMON_WANNIER90_COMMAND=reuse"
)

spawn = source.split("subroutine run_wannier90_seed_files", 1)[1].split(
    "end subroutine run_wannier90_seed_files", 1
)[0]
assert "get_environment_variable" not in spawn
assert "stop " not in spawn.lower()
assert "call execute_wannier90_command" in spawn
assert spawn.count("call comm_bcast") >= 2
assert "call lcfo_sawf_fatal" in spawn

reuse_gate = source.split("subroutine reject_sawf_reuse_collective", 1)[1].split(
    "end subroutine reject_sawf_reuse_collective", 1
)[0]
assert "is_wannier90_reuse_command" in reuse_gate
assert "call lcfo_sawf_fatal" in reuse_gate
assert seed.index("call reject_sawf_reuse_collective") < seed.index(
    "call run_wannier90_seed_files"
)
assert seed.index("call reject_sawf_reuse_collective") < seed.index(
    "call write_wannier90_global_basis_file"
)
checkpoint_import = source.split("subroutine write_wannier90_global_basis_file", 1)[1].split(
    "end subroutine write_wannier90_global_basis_file", 1
)[0]
assert "call read_wannier90_checkpoint_transform" in checkpoint_import
assert source.count("get_environment_variable('SALMON_WANNIER90_COMMAND'") == 1

sawf = (ROOT / "src/gs/dc/lcfo_wannier_sawf.f90").read_text().lower()
off_projection = sawf.split("if (trim(site_symmetry) == 'off') then", 1)[1].split("end if", 1)[0]
assert "'random'" in off_projection

print("PASS Task7 deactivation/reuse/base/dmn/activation/collective-spawn ordering")
