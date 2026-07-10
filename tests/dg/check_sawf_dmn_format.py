#!/usr/bin/env python3
import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[2]
parser = argparse.ArgumentParser()
parser.add_argument("--build-dir", type=Path, default=ROOT / "build-mpi-eigenexa-wannier-lib")
args = parser.parse_args()
BUILD = args.build_dir.resolve()


def parse_dmn(path, nb, nw):
    text = Path(path).read_text()
    tokens = re.findall(r"\([^()]+,[^()]+\)|[-+]?\d+", text.splitlines()[1] + "\n" + "\n".join(text.splitlines()[2:]))
    dims = list(map(int, tokens[:4])); pos = 4
    if dims[0] != nb or dims[2:] != [1, 1]:
        raise AssertionError(f"bad dimensions {dims}")
    nsym = dims[1]
    ik2ir = [int(tokens[pos])]; pos += 1
    ir2ik = [int(tokens[pos])]; pos += 1
    kptsym = [int(tokens[pos+i]) for i in range(nsym)]; pos += nsym
    def zvals(count):
        nonlocal pos
        out=[]
        for tok in tokens[pos:pos+count]:
            a,b=tok.strip("()").split(",")
            out.append(complex(float(a.replace("D","E")),float(b.replace("D","E"))))
        pos += count
        return out
    dw=zvals(nw*nw*nsym); db=zvals(nb*nb*nsym)
    if pos != len(tokens): raise AssertionError("trailing or missing DMN fields")
    return dims,ik2ir,ir2ik,kptsym,dw,db


DRIVER = r'''
program check_sawf_dmn
  use lcfo_wannier_sawf, only: t_sawf_symop
  use lcfo_wannier_sawf_dmn
  implicit none
  type(t_sawf_dmn_writer) :: writer
  type(t_sawf_operation_index) :: operation_index
  type(t_sawf_symop) :: ops(3),ops6(6)
  type(t_sawf_symop) :: map_ops(4)
  complex(8) :: dw2(2,2),db2(2,2),amn2(2,2),bad(2,2)
  complex(8) :: dw3(3,3),db3(3,3),amn3(3,3)
  real(8) :: eig2(2),eig3(3),hres1,hres2,hres3
  logical :: ok
  character(512) :: msg
  integer :: i,j

  call identity_op(ops(1)); call identity_op(ops(2)); call identity_op(ops(3))
  ops(2)%tau=[1d0/3d0,0d0,0d0]; ops(3)%tau=[2d0/3d0,0d0,0d0]
  allocate(ops(1)%atom_map(1),ops(2)%atom_map(1),ops(3)%atom_map(1))
  do i=1,3; ops(i)%atom_map=1; enddo

  dw2=reshape([cmplx(0d0,0d0,8),cmplx(1d0,0d0,8), &
               cmplx(1d0,0d0,8),cmplx(0d0,0d0,8)],[2,2])
  db2=dw2; amn2=reshape([cmplx(1d0,1d0,8),cmplx(2d0,-1d0,8), &
                         cmplx(2d0,-1d0,8),cmplx(1d0,1d0,8)],[2,2])
  eig2=[1d0,1d0]
  call validate_sawf_dmn_covariances(dw2,db2,eig2,amn2,1d-10,ok,msg)
  call require(ok,'covariance pass: '//trim(msg))
  bad=db2; bad(1,2)=cmplx(0.5d0,0.2d0,8)
  call validate_sawf_dmn_covariances(dw2,bad,eig2,amn2,1d-10,ok,msg)
  call require(.not.ok,'AMN/closure failure not detected')
  call validate_sawf_dmn_covariances(dw2,db2,[1d0,2d0],amn2,1d-10,ok,msg)
  call require(.not.ok .and. index(msg,'Hamiltonian')>0,'Hamiltonian covariance failure')
  call begin_sawf_dmn(writer,'covariance_fail.dmn',2,2,1,1d-10,ok,msg)
  call require(ok,'covariance failure begin')
  call append_sawf_dmn_operation(writer,1,dw2,db2,[1d0,2d0],amn2,.true.,ok,msg)
  call require(.not.ok .and. index(msg,'Hamiltonian')>0, &
    'append must propagate covariance failure: '//trim(msg))
  call require(writer%appended==0,'failed covariance must not append an operation')
  call abort_sawf_dmn(writer)
  call validate_sawf_dmn_covariances(dw2,db2,[0d0,2d0],amn2,1d-10,ok,msg, &
    hamiltonian_residual=hres1)
  call require(.not.ok,'nondegenerate mixing baseline must fail')
  call validate_sawf_dmn_covariances(dw2,db2,[100d0,102d0],amn2,1d-10,ok,msg, &
    hamiltonian_residual=hres2)
  call validate_sawf_dmn_covariances(dw2,db2,27.211386245988d0*[0d0,2d0],amn2,1d-10,ok,msg, &
    hamiltonian_residual=hres3)
  call require(abs(hres1-hres2)<1d-12 .and. abs(hres1-hres3)<1d-12, &
    'Hamiltonian residual must be offset/unit invariant')
  amn2(1,1)=amn2(1,1)+cmplx(0.25d0,0.1d0,8)
  call validate_sawf_dmn_covariances(dw2,db2,eig2,amn2,1d-10,ok,msg)
  call require(.not.ok .and. index(msg,'AMN')>0,'orientation-sensitive AMN failure')

  call begin_sawf_dmn(writer,'identity.dmn',1,1,1,1d-10,ok,msg); call require(ok,'identity begin')
  call append_sawf_dmn_operation(writer,1,reshape([cmplx(1d0,0d0,8)],[1,1]), &
    reshape([cmplx(1d0,0d0,8)],[1,1]),[0d0],reshape([cmplx(1d0,0d0,8)],[1,1]),.true.,ok,msg)
  call require(ok,'identity append: '//trim(msg))
  call finish_sawf_dmn(writer,ops(1:1),ok,msg); call require(ok,'identity finish: '//trim(msg))
  call read_dmn_with_sitesym_order('identity.dmn',1,1,1)

  dw3=(0d0,0d0); db3=(0d0,0d0)
  do i=1,3
    dw3(mod(i,3)+1,i)=cmplx(1d0,0d0,8)
    db3(mod(i,3)+1,i)=cmplx(1d0,0d0,8)
  enddo
  amn3=(0d0,0d0); do i=1,3; amn3(i,i)=cmplx(1d0,0d0,8); enddo; eig3=0d0
  call validate_sawf_dmn_covariances(dw3,db3,eig3,amn3,1d-10,ok,msg)
  call require(ok,'forward three-cycle AMN covariance')
  call validate_sawf_dmn_covariances(conjg(transpose(dw3)),db3,eig3,amn3,1d-10,ok,msg)
  call require(.not.ok .and. index(msg,'AMN')>0,'three-cycle transpose orientation not rejected')
  call begin_sawf_dmn(writer,'cycle.dmn',3,3,3,1d-10,ok,msg); call require(ok,'cycle begin')
  call append_sawf_dmn_operation(writer,1,amn3,amn3,eig3,amn3,.true.,ok,msg); call require(ok,'cycle id')
  call append_sawf_dmn_operation(writer,2,dw3,db3,eig3,amn3,.false.,ok,msg); call require(ok,'cycle g')
  call append_sawf_dmn_operation(writer,3,matmul(dw3,dw3),matmul(db3,db3),eig3,amn3,.false.,ok,msg)
  call require(ok,'cycle g2')
  call finish_sawf_dmn(writer,ops,ok,msg); call require(ok,'cycle finish: '//trim(msg))
  call read_dmn_with_sitesym_order('cycle.dmn',3,3,3)

  ! Exact noncommuting S3 representation.  The six coordinate permutations
  ! force the generator traversal to validate both generator orders and cycles.
  call make_s3(ops6)
  call begin_sawf_dmn(writer,'s3.dmn',3,3,6,1d-10,ok,msg); call require(ok,'S3 begin')
  do i=1,6
    do j=1,3; db3(:,j)=cmplx(real(ops6(i)%W(:,j),8),0d0,8); enddo; dw3=db3
    call append_sawf_dmn_operation(writer,i,dw3,db3,eig3,amn3,i==1,ok,msg)
    call require(ok,'S3 append: '//trim(msg))
  enddo
  call finish_sawf_dmn(writer,ops6,ok,msg); call require(ok,'noncommuting S3 finish: '//trim(msg))

  do i=1,4
    call identity_op(map_ops(i)); allocate(map_ops(i)%atom_map(2)); map_ops(i)%atom_map=[1,2]
  enddo
  map_ops(2)%tau=[0.5d0,0d0,0d0]
  map_ops(3)%tau=map_ops(2)%tau; map_ops(3)%atom_map=[2,1]
  call build_sawf_operation_index(map_ops(1:3),1d-10,operation_index,ok,msg)
  call require(ok,'atom-map operation index build: '//trim(msg))
  call lookup_sawf_operation_product(operation_index,map_ops(1:3),3,1,j,ok,msg)
  call require(ok .and. j==3,'group lookup ignored composed atom map')
  map_ops(4)=map_ops(3)
  call build_sawf_operation_index(map_ops,1d-10,operation_index,ok,msg)
  call require(.not.ok .and. index(msg,'duplicate')>0,'duplicate full operation key not rejected')
  call begin_sawf_dmn(writer,'s3_bad.dmn',3,3,6,1d-10,ok,msg); call require(ok,'bad S3 begin')
  do i=1,6
    do j=1,3; db3(:,j)=cmplx(real(ops6(i)%W(:,j),8),0d0,8); enddo
    if(i==4) db3=amn3
    call append_sawf_dmn_operation(writer,i,db3,db3,eig3,amn3,i==1,ok,msg)
    call require(ok,'bad S3 operation must pass local gates')
  enddo
  call finish_sawf_dmn(writer,ops6,ok,msg)
  call require(.not.ok .and. index(msg,'group relation')>0, &
    'noncommuting group failure not detected: '//trim(msg))
  call abort_sawf_dmn(writer)

  open(77,file='stale.dmn',status='replace'); write(77,'(a)') 'OLD'; close(77)
  call begin_sawf_dmn(writer,'stale.dmn',2,2,1,1d-10,ok,msg); call require(ok,'failure begin')
  bad=(0d0,0d0)
  call append_sawf_dmn_operation(writer,1,bad,bad,eig2,reshape([cmplx(1d0,0d0,8), &
    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(1d0,0d0,8)],[2,2]),.true.,ok,msg)
  call require(.not.ok,'invalid operation must fail')
  call abort_sawf_dmn(writer)
  open(78,file='collision.dmn.tmp.4242.0',status='replace'); write(78,'(a)') 'SENTINEL'; close(78)
  call begin_sawf_dmn(writer,'collision.dmn',1,1,1,1d-10,ok,msg,temp_nonce=4242)
  call require(ok,'exclusive temp collision retry: '//trim(msg))
  call append_sawf_dmn_operation(writer,1,reshape([cmplx(1d0,0d0,8)],[1,1]), &
    reshape([cmplx(1d0,0d0,8)],[1,1]),[0d0],reshape([cmplx(1d0,0d0,8)],[1,1]),.true.,ok,msg)
  call finish_sawf_dmn(writer,ops(1:1),ok,msg); call require(ok,'collision finish')
  write(*,'(a)') 'PASS SAWF DMN writer and covariance gates'
contains
  subroutine identity_op(op)
    type(t_sawf_symop),intent(out)::op
    op%W=0; op%W(1,1)=1; op%W(2,2)=1; op%W(3,3)=1; op%R=real(op%W,8); op%tau=0d0
  end subroutine
  subroutine make_s3(group)
    type(t_sawf_symop),intent(out)::group(6)
    integer,parameter::p(3,6)=reshape([1,2,3, 2,1,3, 1,3,2, 2,3,1, 3,1,2, 3,2,1],[3,6])
    integer::g,k
    do g=1,6
      call identity_op(group(g)); group(g)%W=0
      do k=1,3; group(g)%W(p(k,g),k)=1; enddo
      group(g)%R=real(group(g)%W,8); allocate(group(g)%atom_map(1)); group(g)%atom_map=1
    enddo
  end subroutine
  subroutine read_dmn_with_sitesym_order(filename,nb_expected,nw_expected,nsym_expected)
    character(*),intent(in)::filename
    integer,intent(in)::nb_expected,nw_expected,nsym_expected
    integer::u,nb,nsym,nir,nk
    integer,allocatable::ik2ir(:),ir2ik(:),kptsym(:,:)
    complex(8),allocatable::dw(:,:,:,:),db(:,:,:,:)
    open(newunit=u,file=filename,status='old',action='read')
    read(u,*); read(u,*) nb,nsym,nir,nk
    call require(nb==nb_expected .and. nsym==nsym_expected .and. nir==1 .and. nk==1, &
      'sitesym header order')
    allocate(ik2ir(nk),ir2ik(nir),kptsym(nsym,nir))
    allocate(dw(nw_expected,nw_expected,nsym,nir),db(nb,nb,nsym,nir))
    read(u,*) ik2ir; read(u,*) ir2ik; read(u,*) kptsym; read(u,*) dw; read(u,*) db; close(u)
    call require(all(ik2ir==1) .and. all(ir2ik==1) .and. all(kptsym==1), &
      'sitesym Gamma index order')
    call require(maxval(abs(dw(:,:,:,1)-db(:,:,:,1)))<1d-12,'sitesym matrix stream order')
    deallocate(ik2ir,ir2ik,kptsym,dw,db)
  end subroutine
  subroutine require(condition,label)
    logical,intent(in)::condition; character(*),intent(in)::label
    if(.not.condition) then; write(*,'(a)') trim(label); error stop 1; endif
  end subroutine
end program
'''

cache=(BUILD/'CMakeCache.txt').read_text()
fc=next(x.split('=',1)[1] for x in cache.splitlines() if x.startswith('CMAKE_Fortran_COMPILER:'))
with tempfile.TemporaryDirectory(prefix='sawf-dmn-') as td:
    td=Path(td); src=td/'src'; bld=td/'build'; run=td/'run'; src.mkdir(); run.mkdir()
    (src/'config.h').write_text('')
    (src/'sym_stub.f90').write_text("""module sym_sub\ncontains\nsubroutine read_symmetry_file(f,m,o,s)\ncharacter(*),intent(in)::f; real(8),allocatable,intent(out)::m(:,:,:); logical,intent(out)::o; character(*),intent(out)::s\nallocate(m(3,4,0)); o=.false.; s='stub'\nend subroutine\nend module\n""")
    (src/'driver.f90').write_text(DRIVER)
    (src/'CMakeLists.txt').write_text(f"""cmake_minimum_required(VERSION 3.18)
project(sawf_dmn LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_library(sawf {src/'sym_stub.f90'} {ROOT/'src/gs/dc/lcfo_wannier_sawf.f90'} {ROOT/'src/gs/dc/lcfo_wannier_sawf_band.f90'} {ROOT/'src/gs/dc/lcfo_wannier_sawf_dmn.f90'})
target_include_directories(sawf PRIVATE {src})
target_link_libraries(sawf PUBLIC LAPACK::LAPACK)
add_executable(check_dmn {src/'driver.f90'})
target_link_libraries(check_dmn PRIVATE sawf)
""")
    p=subprocess.run(['cmake','-S',src,'-B',bld,f'-DCMAKE_Fortran_COMPILER={fc}'],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if p.returncode: raise SystemExit(p.stdout)
    p=subprocess.run(['cmake','--build',bld,'-j','2'],text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if p.returncode: raise SystemExit(p.stdout)
    p=subprocess.run([bld/'check_dmn'],cwd=run,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    if p.returncode: raise SystemExit(p.stdout)
    print(p.stdout.strip())
    fixture=parse_dmn(ROOT/'tests/dg/fixtures/sawf_gamma_identity.dmn',1,1)
    generated=parse_dmn(run/'identity.dmn',1,1)
    assert fixture[1:]==generated[1:], (fixture,generated)
    cycle=parse_dmn(run/'cycle.dmn',3,3)
    assert cycle[0]==[3,3,1,1] and cycle[1:4]==([1],[1],[1,1,1])
    expected=[0j,1+0j,0j,0j,0j,1+0j,1+0j,0j,0j]
    assert cycle[4][9:18]==expected, cycle[4][9:18]
    assert cycle[5][9:18]==expected, cycle[5][9:18]
    generated_s3=parse_dmn(run/'s3.dmn',3,3)
    fixture_s3=parse_dmn(ROOT/'tests/dg/fixtures/sawf_gamma_s3.dmn',3,3)
    assert generated_s3[1:]==fixture_s3[1:], 'writer-generated S3 differs from static roundtrip fixture'
    assert (run/'stale.dmn').read_text().strip()=='OLD'
    assert (run/'collision.dmn.tmp.4242.0').read_text().strip()=='SENTINEL'
    owned_leftovers=[p for p in run.glob('*.dmn.tmp*') if p.name!='collision.dmn.tmp.4242.0']
    assert not owned_leftovers, owned_leftovers

    w90_root=BUILD/'wannier90/src/wannier90-project'
    w90_binary=w90_root/'wannier90.x'
    w90_library=w90_root/'libwannier.a'
    w90_modules=w90_root/'src/obj'
    if not w90_binary.exists():
        print('SKIP bundled Wannier90 sitesym_read: wannier90.x unavailable')
    else:
        if not w90_library.exists() or not (w90_modules/'w90_sitesym.mod').exists():
            raise SystemExit('wannier90.x exists but bundled sitesym library/modules are unavailable')
        shutil.copyfile(run/'s3.dmn',run/'s3_fixture.dmn')
        reader=run/'w90_sitesym_reader.f90'
        reader.write_text(r'''
program w90_sitesym_fixture
  use w90_io, only: seedname
  use w90_parameters, only: num_bands,num_wann,num_kpts
  use w90_sitesym, only: sitesym_read,sitesym_dealloc,nsymmetry,nkptirr, &
    ik2ir,ir2ik,kptsym,d_matrix_wann,d_matrix_band
  implicit none
  complex(8) :: expected(3,3)
  seedname='s3_fixture'; num_bands=3; num_wann=3; num_kpts=1
  call sitesym_read()
  if(nsymmetry/=6 .or. nkptirr/=1) error stop 1
  if(any(ik2ir/=[1]) .or. any(ir2ik/=[1]) .or. any(kptsym(:,1)/=1)) error stop 2
  expected=reshape([cmplx(0d0,0d0,8),cmplx(1d0,0d0,8),cmplx(0d0,0d0,8), &
                    cmplx(0d0,0d0,8),cmplx(0d0,0d0,8),cmplx(1d0,0d0,8), &
                    cmplx(1d0,0d0,8),cmplx(0d0,0d0,8),cmplx(0d0,0d0,8)],[3,3])
  if(maxval(abs(d_matrix_wann(:,:,4,1)-expected))>1d-12) error stop 3
  if(maxval(abs(d_matrix_band(:,:,4,1)-expected))>1d-12) error stop 4
  if(maxval(abs(d_matrix_wann(:,:,4,1)-conjg(transpose(expected))))<0.5d0) error stop 5
  call sitesym_dealloc()
  write(*,'(a)') 'PASS bundled Wannier90 sitesym_read nontrivial S3 fixture'
end program
''')
        reader_exe=run/'w90_sitesym_reader'
        compile_cmd=[fc,'-I',w90_modules,reader,w90_library,'-o',reader_exe]
        if sys.platform=='darwin': compile_cmd += ['-framework','Accelerate']
        else: compile_cmd += ['-llapack','-lblas']
        p=subprocess.run([str(x) for x in compile_cmd],cwd=run,text=True,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if p.returncode: raise SystemExit('Bundled Wannier90 sitesym reader build failed:\n'+p.stdout)
        p=subprocess.run([reader_exe],cwd=run,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if p.returncode: raise SystemExit('Bundled Wannier90 sitesym reader failed:\n'+p.stdout)
        print(p.stdout.strip())

        actual=run/'w90_actual'; actual.mkdir()
        seed='s3_actual'
        (actual/f'{seed}.win').write_text('''num_bands = 3
num_wann = 3
num_iter = 2
site_symmetry = true
symmetrize_eps = 1d-8
gamma_only = true
mp_grid = 1 1 1
begin unit_cell_cart
ang
1.0 0.0 0.0
0.0 1.0 0.0
0.0 0.0 1.0
end unit_cell_cart
begin atoms_frac
H 0.0 0.0 0.0
end atoms_frac
begin projections
random
end projections
begin kpoints
0.0 0.0 0.0
end kpoints
''')
        p=subprocess.run([w90_binary,'-pp',seed],cwd=actual,text=True,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if p.returncode: raise SystemExit('Bundled wannier90.x -pp failed:\n'+p.stdout)
        nnkp=(actual/f'{seed}.nnkp').read_text().splitlines()
        begin=next(i for i,x in enumerate(nnkp) if x.strip().lower()=='begin nnkpts')
        end=next(i for i,x in enumerate(nnkp[begin+1:],begin+1) if x.strip().lower()=='end nnkpts')
        nnkpts_lines=[x.strip() for x in nnkp[begin+1:end] if x.strip()]
        if not nnkpts_lines:
            raise SystemExit('Bundled wannier90.x generated an empty nnkpts block')
        nntot=int(nnkpts_lines[0].split()[0])
        neighbors=[x.split() for x in nnkpts_lines[1:]]
        if len(neighbors)!=nntot or any(len(fields)<5 for fields in neighbors):
            raise SystemExit(f'Unexpected nnkp nnkpts block: nntot={nntot} neighbors={neighbors!r}')
        mmn=['synthetic identity overlaps',f'3 1 {nntot}']
        for fields in neighbors:
            mmn.append(' '.join(fields[:5]))
            for col in range(3):
                for row in range(3):
                    mmn.append(f'{1.0 if row==col else 0.0:.16e} 0.0')
        (actual/f'{seed}.mmn').write_text('\n'.join(mmn)+'\n')
        amn=['synthetic identity projections','3 1 3']
        for iw in range(1,4):
            for ib in range(1,4):
                amn.append(f'{ib} {iw} 1 {1.0 if ib==iw else 0.0:.16e} 0.0')
        (actual/f'{seed}.amn').write_text('\n'.join(amn)+'\n')
        (actual/f'{seed}.eig').write_text('\n'.join(f'{ib} 1 0.0' for ib in range(1,4))+'\n')
        shutil.copyfile(run/'s3.dmn',actual/f'{seed}.dmn')
        p=subprocess.run([w90_binary,seed],cwd=actual,text=True,
            stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        if p.returncode: raise SystemExit('Bundled wannier90.x site_symmetry run failed:\n'+p.stdout)
        print('PASS bundled wannier90.x site_symmetry nontrivial S3 run')

source=(ROOT/'src/gs/dc/lcfo_flux.f90').read_text()
seed=source.split('subroutine write_wannier_seed_files()',1)[1].split('end subroutine write_wannier_seed_files',1)[0]
amn_index=seed.index('call write_wannier_amn_file')
dmn_index=seed.index('call generate_sawf_dmn')
activation_index=seed.index('call activate_sawf_win_collective')
export_index=seed.index('is_wannier90_export_only_command(resolved_wannier_command)')
spawn_index=seed.index('call run_wannier90_seed_files(resolved_wannier_command)')
assert amn_index < dmn_index < activation_index < export_index < spawn_index
export_block=seed[export_index:spawn_index]
assert 'return' in export_block
assert "trim(wannier_site_symmetry) /= 'off'" in seed
base_writer=source.split('subroutine write_wannier_base_win_atomic',1)[1].split(
    'end subroutine write_wannier_base_win_atomic',1)[0]
assert 'site_symmetry' not in base_writer and 'symmetrize_eps' not in base_writer
generator=source.split('subroutine generate_sawf_dmn',1)[1].split('end subroutine generate_sawf_dmn',1)[0]
assert generator.index('prepare_sawf_fragment_state_cache') < generator.index('do iop=1,size(symmetry_operations)')
builder=source.split('subroutine build_sawf_dmn_operation_fragment_local',1)[1].split(
    'end subroutine build_sawf_dmn_operation_fragment_local',1)[0]
assert 'read_fragment_lcfo_seed_for_wannier_import' not in builder
assert 'source_seed_reads' in generator and 'source_reconstructions' in generator
assert 'sawf_target_cache_capacity=2' in source.replace(' ','').lower()
module=(ROOT/'src/gs/dc/lcfo_wannier_sawf_dmn.f90').read_text().lower()
assert 'd_band(:,:,:)' not in module and 'd_wann(:,:,:)' not in module
assert 'status=' in module and 'rename' in module and 'scratch' in module
assert 'file_storage_size' in module and 'checked_stream_position' in module
assert 'find_operation_product' not in module
lookup=module.split('subroutine find_index_matches',1)[1].split(
    'end subroutine find_index_matches',1)[0]
assert 'allocate(' not in lookup and 'size(operations)' not in lookup
print('PASS Task6 DMN/activation/export/spawn order and streaming source guard')
