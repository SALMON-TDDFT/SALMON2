#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("--build-dir", type=Path, default=ROOT / "build-mpi-eigenexa-wannier-lib")
args = parser.parse_args()
BUILD = args.build_dir.resolve()

DRIVER = r"""
program check_sawf_band
  use lcfo_wannier_sawf, only: t_sawf_symop, t_sawf_projection_channel, &
    build_sawf_spd_projection_map, build_sawf_wannier_representation
  use lcfo_wannier_sawf_band
  implicit none
  type(t_sawf_symop) :: id, shift, inv, g1, g2, gp, swap
  type(t_sawf_symop),allocatable :: atom_op(:)
  type(t_sawf_projection_channel),allocatable :: channels(:)
  integer, allocatable :: map(:), m1(:), m2(:), mp(:)
  integer :: mesh(3), p, owner, local(3), mapped(3), origin(3,3), shape_frag(3,3)
  integer :: nmat_meta(1),nbasis_meta(2,1)
  integer :: global_index_meta(3,2,1)
  complex(8) :: psi3(3,3), psi4(4,4), d3(3,3), expected3(3,3)
  complex(8) :: d1(4,4), d2(4,4), dp(4,4), dleft(4,4), dright(4,4)
  complex(8) :: dset(4,4,3)
  complex(8) :: done(1,1), mixed(4,1)
  real(8),allocatable :: d_wann(:,:,:)
  real(8) :: smin,smax,closure,residual
  real(8) :: source_basis(2,2),target_basis(2,2),target_partial(2,1)
  real(8) :: leakage,residual_norm2,transformed_norm2
  logical :: ok
  character(512) :: msg

  call set_identity(id)
  mesh=[4,1,1]
  call build_sawf_periodic_grid_map(id,mesh,1d-8,map,ok,msg)
  call require(ok .and. all(map==[1,2,3,4]),'identity grid map')
  psi4=(0d0,0d0); do p=1,4; psi4(p,p)=1d0; enddo
  call accumulate_sawf_dmn_band(psi4,map,1d0,1,4,d1,ok,msg)
  call require(ok .and. maxval(abs(d1-identity4()))<1d-13,'identity D_band')

  inv=id; inv%W=-id%W
  call build_sawf_periodic_grid_map(inv,mesh,1d-8,map,ok,msg)
  call require(ok,'inversion grid map')
  call accumulate_sawf_dmn_band(psi4,map,1d0,1,4,d1,ok,msg)
  call require(ok .and. maxval(abs(d1-permutation_from_map(map)))<1d-13,'inversion D_band')

  ! Decisive forward-orientation oracle: a 3-cycle localized basis must give
  ! D(target,source)=1, identical to the active Wannier permutation.
  mesh=[3,1,1]; shift=id; shift%tau=[1d0/3d0,0d0,0d0]
  call build_sawf_periodic_grid_map(shift,mesh,1d-8,map,ok,msg)
  call require(ok .and. all(map==[2,3,1]),'three-cycle map')
  psi3=(0d0,0d0); do p=1,3; psi3(p,p)=1d0; enddo
  call accumulate_sawf_dmn_band(psi3,map,1d0,1,3,d3,ok,msg)
  expected3=(0d0,0d0); do p=1,3; expected3(map(p),p)=1d0; enddo
  call require(ok .and. maxval(abs(d3-expected3))<1d-13,'D_band is not forward P')
  call require(maxval(abs(d3-conjg(transpose(expected3))))>0.5d0,'3-cycle does not reject Mdagger')
  call build_sawf_spd_projection_map(3,channels,ok,msg); call require(ok,'three-atom spd map')
  allocate(atom_op(1)); allocate(atom_op(1)%atom_map(3))
  atom_op(1)%W=id%W; atom_op(1)%R=real(id%W,8); atom_op(1)%atom_map=[2,3,1]
  call build_sawf_wannier_representation(atom_op,channels,d_wann,ok,msg)
  call require(ok,'synthetic D_wann')
  call require(maxval(abs(d3-d_wann([1,10,19],[1,10,19],1)))<1d-13, &
    'D_band/D_wann forward intertwining mismatch when bases coincide')

  ! Noncommuting product: g1 swaps x/y, g2 translates x.  Active maps and
  ! representations must both satisfy g2*g1 and D(g2)D(g1).
  mesh=[2,2,1]
  g1=id; g1%W=reshape([0,1,0,1,0,0,0,0,1],[3,3])
  g2=id; g2%tau=[0.5d0,0d0,0d0]
  gp%W=matmul(g2%W,g1%W); gp%tau=matmul(real(g2%W,8),g1%tau)+g2%tau
  call build_sawf_periodic_grid_map(g1,mesh,1d-8,m1,ok,msg); call require(ok,'g1 map')
  call build_sawf_periodic_grid_map(g2,mesh,1d-8,m2,ok,msg); call require(ok,'g2 map')
  call build_sawf_periodic_grid_map(gp,mesh,1d-8,mp,ok,msg); call require(ok,'product map')
  do p=1,4; call require(mp(p)==m2(m1(p)),'active map product order'); enddo
  call accumulate_sawf_dmn_band(psi4,m1,1d0,1,4,d1,ok,msg); call require(ok,'Dg1')
  call accumulate_sawf_dmn_band(psi4,m2,1d0,1,4,d2,ok,msg); call require(ok,'Dg2')
  call accumulate_sawf_dmn_band(psi4,mp,1d0,1,4,dp,ok,msg); call require(ok,'Dproduct')
  call validate_sawf_representation_product(d2,d1,dp,1d-12,residual,ok,msg)
  call require(ok,'noncommuting representation product: '//trim(msg))
  dset(:,:,1)=d1; dset(:,:,2)=d2; dset(:,:,3)=dp
  call validate_sawf_operation_set_products(dset,[2],[1],[3],1d-12,residual,ok,msg)
  call require(ok,'operation-set representation validation: '//trim(msg))
  call accumulate_sawf_dmn_band(psi4,mp,1d0,1,2,dleft,ok,msg); call require(ok,'partition 1')
  call accumulate_sawf_dmn_band(psi4,mp,1d0,3,4,dright,ok,msg); call require(ok,'partition 2')
  call require(maxval(abs(dleft+dright-dp))<1d-13,'partition sum')

  call validate_sawf_dmn_band(dp,1d-10,smin,smax,closure,ok,msg)
  call require(ok .and. abs(smin-1d0)<1d-12 .and. abs(smax-1d0)<1d-12,'closed SVD')
  mixed(:,1)=(psi4(:,1)+psi4(:,2))/sqrt(2d0)
  call accumulate_sawf_dmn_band(mixed,mp,1d0,1,4,done,ok,msg); call require(ok,'nonclosed overlap')
  call validate_sawf_dmn_band(done,1d-10,smin,smax,closure,ok,msg)
  call require(.not.ok .and. index(msg,'not symmetry-closed')>0,'nonclosed rejection')

  ! Synthetic fragment tiling; swap maps one x-fragment across both targets.
  mesh=[4,4,1]
  origin(:,1)=[0,0,0]; origin(:,2)=[2,0,0]; origin(:,3)=[0,0,0]
  shape_frag(:,1)=[2,4,1]; shape_frag(:,2)=[2,4,1]; shape_frag(:,3)=[1,1,1]
  call validate_sawf_fragment_tiling(mesh,origin(:,1:2),shape_frag(:,1:2),ok,msg)
  call require(ok,'valid fragment tiling: '//trim(msg))
  swap=id; swap%W=reshape([0,1,0,1,0,0,0,0,1],[3,3])
  call map_sawf_periodic_grid_point(swap,mesh,1d-8,[1,1,1],mapped,ok,msg)
  call require(ok,'mapped fragment point 1')
  call locate_sawf_fragment_point(mapped,mesh,origin(:,1:2),shape_frag(:,1:2),owner,local,ok,msg)
  call require(ok .and. owner==1,'first target fragment')
  call map_sawf_periodic_grid_point(swap,mesh,1d-8,[1,4,1],mapped,ok,msg)
  call require(ok,'mapped fragment point 2')
  call locate_sawf_fragment_point(mapped,mesh,origin(:,1:2),shape_frag(:,1:2),owner,local,ok,msg)
  call require(ok .and. owner==2,'source fragment must span target fragments')
  call validate_sawf_fragment_tiling(mesh,origin(:,1:2),reshape([2,4,1,1,4,1],[3,2]),ok,msg)
  call require(.not.ok .and. index(msg,'missing')>0,'missing ownership rejection')
  call validate_sawf_fragment_tiling(mesh,reshape([0,0,0,1,0,0],[3,2]), &
    reshape([3,4,1,3,4,1],[3,2]),ok,msg)
  call require(.not.ok .and. index(msg,'duplicate')>0,'duplicate ownership rejection')

  shift%tau=[0.1d0,0d0,0d0]
  call build_sawf_periodic_grid_map(shift,[4,1,1],1d-8,map,ok,msg)
  call require(.not.ok .and. index(msg,'translation')>0,'incompatible tau rejection')
  call validate_sawf_seed_header(2,1,4,1,1,8,8,ok,msg)
  call require(.not.ok .and. index(msg,'fragment count')>0,'mixed fragment metadata rejection')
  call validate_sawf_seed_header(2,1,4,2,1,8,3,ok,msg)
  call require(.not.ok .and. index(msg,'insufficient')>0,'stale insufficient-band seed rejection')
  call validate_sawf_seed_header(2,1,4,2,1,8,8,ok,msg); call require(ok,'valid seed header')

  source_basis=reshape([1d0,0d0,0d0,1d0],[2,2])
  target_basis=source_basis
  call diagnose_sawf_fragment_basis_closure(source_basis,target_basis,[1,2],[1,2],1d0, &
    leakage,residual_norm2,transformed_norm2,ok,msg)
  call require(ok .and. leakage<1d-15 .and. residual_norm2<1d-15 .and. &
    abs(transformed_norm2-2d0)<1d-15,'identity basis closure diagnostic')
  target_partial(:,1)=[1d0,0d0]
  call diagnose_sawf_fragment_basis_closure(source_basis,target_partial,[1,2],[1,2],1d0, &
    leakage,residual_norm2,transformed_norm2,ok,msg)
  call require(ok .and. abs(leakage-0.5d0)<1d-15 .and. abs(residual_norm2-1d0)<1d-15 .and. &
    abs(transformed_norm2-2d0)<1d-15,'leaking basis closure diagnostic')
  nmat_meta=[8]; nbasis_meta(:,1)=[9,4]
  call validate_sawf_seed_basis_metadata(8,8,nmat_meta,nbasis_meta,ok,msg)
  call require(.not.ok .and. index(msg,'exceeds')>0,'basis metadata overflow rejection')
  nmat_meta=[6]; nbasis_meta(:,1)=[3,3]
  global_index_meta(:,:,1)=reshape([1,2,3,4,5,6],[3,2])
  call validate_sawf_seed_basis_metadata(3,4,nmat_meta,nbasis_meta,ok,msg,global_index_meta)
  call require(ok,'later-fragment global indices above nstate_frag must be accepted')
  global_index_meta(3,2,1)=7
  call validate_sawf_seed_basis_metadata(3,4,nmat_meta,nbasis_meta,ok,msg,global_index_meta)
  call require(.not.ok .and. index(msg,'exceeds n_mat')>0,'global index above n_mat rejection')
  write(*,'(a)') 'PASS SAWF band representation orientation, closure, and fragment ownership'
contains
  subroutine set_identity(op)
    type(t_sawf_symop),intent(out)::op
    op%W=0; op%W(1,1)=1; op%W(2,2)=1; op%W(3,3)=1; op%tau=0d0
  end subroutine
  function identity4() result(a)
    complex(8)::a(4,4); integer::k
    a=(0d0,0d0); do k=1,4; a(k,k)=1d0; enddo
  end function
  function permutation_from_map(pm) result(a)
    integer,intent(in)::pm(:); complex(8)::a(size(pm),size(pm)); integer::k
    a=(0d0,0d0); do k=1,size(pm); a(pm(k),k)=1d0; enddo
  end function
  subroutine require(cond,text)
    logical,intent(in)::cond; character(*),intent(in)::text
    if(.not.cond) then; write(*,'(a)') trim(text); error stop 1; endif
  end subroutine
end program
"""

MPI_DRIVER = r"""
program check_sawf_band_mpi
  use mpi
  use lcfo_wannier_sawf, only: t_sawf_symop
  use lcfo_wannier_sawf_band, only: build_sawf_periodic_grid_map, accumulate_sawf_dmn_band, &
    accumulate_sawf_dmn_band_blocks
  implicit none
  type(t_sawf_symop)::op
  integer,allocatable::map(:)
  complex(8)::psi(4,4),local(4,4),total(4,4),serial(4,4)
  complex(8)::source_block(2,2),target1(1,2),target2(1,2)
  complex(8)::block1(2,2),block2(2,2),block_local(2,2),block_total(2,2)
  integer,allocatable::empty(:)
  integer::i,rank,nproc,ierr,first,last
  logical::ok; character(256)::msg
  call MPI_Init(ierr); call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr); call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  if(nproc<2) then; if(rank==0) write(*,'(a)') 'MPI test requires at least two ranks'; call MPI_Abort(MPI_COMM_WORLD,2,ierr); endif
  op%W=0; op%W(1,1)=1; op%W(2,2)=1; op%W(3,3)=1; op%tau=[0.25d0,0d0,0d0]
  call build_sawf_periodic_grid_map(op,[4,1,1],1d-8,map,ok,msg)
  psi=(0d0,0d0); do i=1,4; psi(i,i)=1d0; enddo
  first=(4*rank)/nproc+1; last=(4*(rank+1))/nproc
  call accumulate_sawf_dmn_band(psi,map,1d0,first,last,local,ok,msg)
  call MPI_Allreduce(local,total,16,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  call accumulate_sawf_dmn_band(psi,map,1d0,1,4,serial,ok,msg)
  if(maxval(abs(total-serial))>1d-13) call MPI_Abort(MPI_COMM_WORLD,3,ierr)
  source_block=(0d0,0d0); source_block(1,1)=1d0; source_block(2,2)=1d0
  target1(1,:)=[cmplx(1d0,0d0,8),cmplx(0d0,0d0,8)]
  target2(1,:)=[cmplx(0d0,0d0,8),cmplx(1d0,0d0,8)]
  allocate(empty(0)); block_local=(0d0,0d0)
  if(rank==0) then
    call accumulate_sawf_dmn_band_blocks(source_block,target1,[1],[1],1d0,block1,ok,msg)
    call accumulate_sawf_dmn_band_blocks(source_block,target2,[2],[1],1d0,block2,ok,msg)
    block_local=block1+block2
  else
    call accumulate_sawf_dmn_band_blocks(source_block,target1,empty,empty,1d0,block_local,ok,msg)
  endif
  call MPI_Allreduce(block_local,block_total,4,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,ierr)
  if(maxval(abs(block_total-reshape([cmplx(1d0,0d0,8),cmplx(0d0,0d0,8), &
      cmplx(0d0,0d0,8),cmplx(1d0,0d0,8)],[2,2])))>1d-13) &
    call MPI_Abort(MPI_COMM_WORLD,4,ierr)
  if(rank==0) write(*,'(a)') 'PASS SAWF MPI partitioned D_band accumulation'
  call MPI_Finalize(ierr)
end program
"""

def run(cmd, **kwargs):
    return subprocess.run([str(x) for x in cmd], text=True, stdout=subprocess.PIPE,
                          stderr=subprocess.STDOUT, **kwargs)

cache=(BUILD/'CMakeCache.txt').read_text()
fc=next(x.split('=',1)[1] for x in cache.splitlines() if x.startswith('CMAKE_Fortran_COMPILER:'))
with tempfile.TemporaryDirectory(prefix='sawf-band-cmake-') as td:
    td=Path(td); src=td/'src'; bld=td/'build'; src.mkdir()
    (src/'config.h').write_text('')
    (src/'sym_stub.f90').write_text("""module sym_sub\ncontains\nsubroutine read_symmetry_file(f,m,o,s)\ncharacter(*),intent(in)::f; real(8),allocatable,intent(out)::m(:,:,:); logical,intent(out)::o; character(*),intent(out)::s\nallocate(m(3,4,0)); o=.false.; s='stub'\nend subroutine\nend module\n""")
    (src/'driver.f90').write_text(DRIVER); (src/'mpi_driver.f90').write_text(MPI_DRIVER)
    (src/'CMakeLists.txt').write_text(f"""cmake_minimum_required(VERSION 3.18)
project(sawf_band LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_library(sawf {src/'sym_stub.f90'} {ROOT/'src/gs/dc/lcfo_wannier_sawf.f90'} {ROOT/'src/gs/dc/lcfo_wannier_sawf_band.f90'})
target_include_directories(sawf PRIVATE {src})
if(TARGET LAPACK::LAPACK)
  target_link_libraries(sawf PUBLIC LAPACK::LAPACK)
else()
  target_link_libraries(sawf PUBLIC ${{LAPACK_LIBRARIES}})
endif()
add_executable(check_band {src/'driver.f90'})
target_link_libraries(check_band PRIVATE sawf)
find_package(MPI COMPONENTS Fortran)
if(MPI_Fortran_FOUND)
  add_executable(check_band_mpi {src/'mpi_driver.f90'})
  target_link_libraries(check_band_mpi PRIVATE sawf MPI::MPI_Fortran)
endif()
""")
    result=run(['cmake','-S',src,'-B',bld,f'-DCMAKE_Fortran_COMPILER={fc}'])
    if result.returncode: raise SystemExit(result.stdout)
    result=run(['cmake','--build',bld,'-j','2'])
    if result.returncode: raise SystemExit(result.stdout)
    result=run([bld/'check_band'])
    if result.returncode: raise SystemExit(result.stdout)
    print(result.stdout.strip())
    mpi_exe=bld/'check_band_mpi'
    mpiexec=None
    match=re.search(r'^MPIEXEC_EXECUTABLE:FILEPATH=(.+)$',(bld/'CMakeCache.txt').read_text(),re.M)
    if match and Path(match.group(1)).exists(): mpiexec=match.group(1)
    if mpi_exe.exists() and mpiexec:
      result=run([mpiexec,'-n','2',mpi_exe])
      if result.returncode: raise SystemExit(result.stdout)
      print(result.stdout.strip())
    else:
      print('SKIP SAWF MPI focused test: MPI Fortran or mpiexec unavailable')

base=(ROOT/'src/gs/dc/lcfo_wannier_sawf.f90').read_text().lower()
if 'zgesvd' in base or 'validate_sawf_dmn_band' in base:
    raise SystemExit('Task5 LAPACK dependency leaked into base SAWF module')
print('PASS base SAWF module remains LAPACK-free')

flux=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text()
output=flux.split('subroutine output',1)[1].split('end subroutine output',1)[0]
seed_write=output.find("filename = trim(base_directory)//binfile_wf_wannier_seed")
provenance_write=output.find("wavefunctions_wannier_seed.provenance")
sync=output.find("call comm_sync_all(dc%icomm_tot)")
global_seed=output.find("call write_wannier_seed_files()")
if min(seed_write,provenance_write,sync,global_seed)<0 or not (seed_write<provenance_write<sync<global_seed):
    raise SystemExit('current seed/provenance/sync/global seed call order regression')
seed_routine=flux.split('subroutine write_wannier_seed_files()',1)[1].split(
    'end subroutine write_wannier_seed_files',1)[0]
if 'call prepare_sawf_dmn_band' in seed_routine:
    raise SystemExit('Task5 must not invoke production D_band before Task6 group/AMN gate')
if 'write_wannier_amn_file' not in seed_routine or 'strict_current_run' not in flux:
    raise SystemExit('Task6 provenance/AMN prerequisites are missing')
print('PASS current-run seed provenance and deferred Task6 call order')
