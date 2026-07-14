#!/usr/bin/env python3
"""Compile and exercise the real SAWF fragment-symmetry validator."""

from pathlib import Path
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[2]
BUILD = ROOT / "build-mpi-eigenexa-wannier-lib"

DRIVER = r"""
program check_sawf_fragment_symmetry_map
  use lcfo_wannier_sawf, only: t_sawf_symop
  use lcfo_wannier_sawf_band, only: validate_sawf_fragment_symmetry_map, &
    build_sawf_fragment_buffer_point_map
  implicit none
  type(t_sawf_symop) :: op
  integer :: mesh(3), origin(3,8), shape(3,8), buffer(3)
  integer, allocatable :: source_to_target(:),buffer_map(:)
  integer :: max_targets, ix, iy, iz, ifrag
  real(8) :: residual, center(3)
  logical :: grid_ok, fragment_ok, center_available
  character(256) :: message

  mesh=[32,32,32]
  shape=16
  ifrag=0
  do iz=0,1
    do iy=0,1
      do ix=0,1
        ifrag=ifrag+1
        origin(:,ifrag)=[16*ix,16*iy,16*iz]
      end do
    end do
  end do
  buffer=[6,6,6]

  call set_identity(op)
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  call require(grid_ok .and. fragment_ok,'identity must be accepted')
  call require(max_targets==1 .and. all(source_to_target==[(ifrag,ifrag=1,8)]), &
    'identity fragment map must be the identity')
  call require(.not.center_available,'identity must not report an inversion center')

  call set_inversion(op,15d0/32d0)
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  call require(grid_ok .and. fragment_ok,'q15 inversion must be accepted')
  call require(max_targets==1 .and. all(source_to_target==[(ifrag,ifrag=1,8)]), &
    'q15 inversion must map every fragment to itself')
  call require(center_available .and. maxval(abs(center-7.5d0))<1d-12, &
    'q15 inversion center must be the half-grid point 7.5')
  call build_sawf_fragment_buffer_point_map(op,mesh,origin(:,1),shape(:,1), &
    origin(:,1),shape(:,1),buffer,1d-10,buffer_map,grid_ok,message)
  call require(grid_ok .and. size(buffer_map)==28**3,'q15 buffer map')
  call require(buffer_map(1)==28**3 .and. buffer_map(28**3)==1, &
    'q15 extended buffer must be a continuous reversal')

  call set_inversion(op,1d0/8d0)
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  call require(grid_ok .and. .not.fragment_ok,'q4 inversion must fail whole-fragment mapping')
  call require(max_targets==8,'q4 inversion must touch all eight target fragments')

  call set_inversion(op,31d0/32d0)
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  call require(grid_ok .and. fragment_ok,'q31 inversion must be a valid fragment permutation')
  call require(max_targets==1 .and. all(source_to_target==[8,7,6,5,4,3,2,1]), &
    'q31 inversion must swap opposite fragments')
  call require(center_available .and. maxval(abs(center-15.5d0))<1d-12, &
    'q31 inversion center must remain half-grid valued')
  call build_sawf_fragment_buffer_point_map(op,mesh,origin(:,1),shape(:,1), &
    origin(:,8),shape(:,8),buffer,1d-10,buffer_map,grid_ok,message)
  call require(grid_ok .and. buffer_map(1)==28**3 .and. buffer_map(28**3)==1, &
    'q31 opposite-fragment buffer reversal')

  call set_inversion(op,1d0/10d0)
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  call require(.not.grid_ok .and. .not.fragment_ok,'nonintegral grid translation must be rejected')

  op%W=0
  op%W(1,2)=1; op%W(2,1)=1; op%W(3,3)=1; op%tau=0d0
  buffer=[6,7,6]
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  call require(grid_ok .and. .not.fragment_ok,'axis-exchange buffer mismatch must be rejected')

  call set_identity(op)
  op%W(1,1)=huge(0)
  buffer=[6,6,6]
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  call require(.not.grid_ok .and. .not.fragment_ok, &
    'pathological huge W must be rejected without unsafe integer rounding')

  call set_identity(op)
  buffer=[6,6,6]
  origin(:,2)=origin(:,1)
  call validate_sawf_fragment_symmetry_map(op,mesh,origin,shape,buffer,1d-10, &
    grid_ok,fragment_ok,max_targets,source_to_target,residual,center_available,center,message)
  call require(grid_ok .and. .not.fragment_ok,'overlapping fragment tiling must be rejected')

  write(*,'(a)') 'PASS SAWF fragment symmetry map validation'

contains
  subroutine set_identity(operation)
    type(t_sawf_symop), intent(out) :: operation
    operation%W=0
    operation%W(1,1)=1; operation%W(2,2)=1; operation%W(3,3)=1
    operation%tau=0d0
  end subroutine set_identity

  subroutine set_inversion(operation,tau)
    type(t_sawf_symop), intent(out) :: operation
    real(8), intent(in) :: tau
    operation%W=0
    operation%W(1,1)=-1; operation%W(2,2)=-1; operation%W(3,3)=-1
    operation%tau=tau
  end subroutine set_inversion

  subroutine require(condition,text)
    logical, intent(in) :: condition
    character(*), intent(in) :: text
    if(.not.condition) then
      write(*,'(2a)') 'FAIL: ',trim(text)
      error stop 1
    end if
  end subroutine require
end program check_sawf_fragment_symmetry_map
"""


def run(command, **kwargs):
    return subprocess.run(
        [str(item) for item in command],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        **kwargs,
    )


def main():
    cache_text = (BUILD / "CMakeCache.txt").read_text()
    compiler = next(
        line.split("=", 1)[1]
        for line in cache_text.splitlines()
        if line.startswith("CMAKE_Fortran_COMPILER:")
    )
    with tempfile.TemporaryDirectory(prefix="sawf-fragment-map-") as temporary:
        temporary = Path(temporary)
        source = temporary / "src"
        build = temporary / "build"
        source.mkdir()
        (source / "config.h").write_text("")
        (source / "sym_stub.f90").write_text(
            """module sym_sub
contains
subroutine read_symmetry_file(f,m,o,s)
character(*),intent(in)::f
real(8),allocatable,intent(out)::m(:,:,:)
logical,intent(out)::o
character(*),intent(out)::s
allocate(m(3,4,0)); o=.false.; s='stub'
end subroutine read_symmetry_file
end module sym_sub
"""
        )
        (source / "driver.f90").write_text(DRIVER)
        (source / "CMakeLists.txt").write_text(
            f"""cmake_minimum_required(VERSION 3.18)
project(sawf_fragment_map LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_library(sawf
  {source / 'sym_stub.f90'}
  {ROOT / 'src/gs/dc/lcfo_wannier_sawf.f90'}
  {ROOT / 'src/gs/dc/lcfo_wannier_sawf_band.f90'})
target_include_directories(sawf PRIVATE {source})
if(TARGET LAPACK::LAPACK)
  target_link_libraries(sawf PUBLIC LAPACK::LAPACK)
else()
  target_link_libraries(sawf PUBLIC ${{LAPACK_LIBRARIES}})
endif()
add_executable(check_fragment_map {source / 'driver.f90'})
target_link_libraries(check_fragment_map PRIVATE sawf)
"""
        )
        result = run(["cmake", "-S", source, "-B", build, f"-DCMAKE_Fortran_COMPILER={compiler}"])
        if result.returncode:
            raise SystemExit(result.stdout)
        result = run(["cmake", "--build", build, "-j", "2"])
        if result.returncode:
            raise SystemExit(result.stdout)
        result = run([build / "check_fragment_map"])
        if result.returncode:
            raise SystemExit(result.stdout)
        print(result.stdout.strip())

    flux_source = (ROOT / "src/gs/dc/lcfo_flux.f90").read_text()
    flux = flux_source.split("subroutine generate_sawf_dmn", 1)[1].split(
        "end subroutine generate_sawf_dmn", 1
    )[0]
    load = flux.index("call put_sawf_identity_first")
    validate = flux.index("call validate_sawf_fragment_symmetry_map", load)
    project = flux.index("call build_sawf_spd_projection_map", validate)
    publish = flux.index("call begin_sawf_dmn", validate)
    cache = flux.index("call prepare_sawf_fragment_state_cache", validate)
    assert load < validate < project < publish < cache
    assert "[DC-LCFO-SAWF-ALIGN]" in flux
    assert "call reduce_sawf_fragment_alignment_failure" in flux[validate:project]
    print("PASS SAWF runtime fragment-map validation ordering")

    band = (ROOT / "src/gs/dc/lcfo_wannier_sawf_band.f90").read_text()
    helper = band.index("subroutine validate_sawf_fragment_symmetry_map")
    prepare = band.index("call prepare_grid_operation", helper)
    residual = band.index("! Residual is measured", helper)
    assert prepare < residual
    assert "nint(scaled_value)" not in band[helper:band.index("end subroutine validate_sawf_fragment_symmetry_map", helper)]
    print("PASS SAWF residual diagnostics follow validated integer preparation")


if __name__ == "__main__":
    main()
