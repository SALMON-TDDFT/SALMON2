#!/usr/bin/env python3
from pathlib import Path
import os, shutil, subprocess, tempfile

ROOT=Path(__file__).resolve().parents[2]
FC=shutil.which('gfortran')
if not FC: raise SystemExit('gfortran is required')
driver=r'''program check_gauge_apply
 use lcfo_wannier_sawf_templates, only: apply_sawf_gauge_connection, &
   stitch_and_apply_sawf_neighbor_pair
 implicit none
 complex(8)::q(2,2),basis(3,2),ww(2,2),wp(2,1),face_self(2,2),face_neighbor(2,2)
 complex(8)::overlap(2,2),neighbor_q(2,2),gauge(2,2),basis_before(3,2)
 complex(8)::b0(3,2),w0(2,2),p0(2,1),f0(2,2),n0(2,2)
 real(8)::c,s,residual
 logical::ok
 character(256)::msg
 integer::i
 c=cos(.31d0);s=sin(.31d0);q=reshape([cmplx(c,0d0,8),cmplx(s,0d0,8), &
   cmplx(-s,0d0,8),cmplx(c,0d0,8)],[2,2])
 b0=reshape([(cmplx(dble(i),0d0,8),i=1,6)],[3,2]); basis=b0
 w0=reshape([cmplx(2d0,0d0,8),cmplx(.2d0,-.1d0,8), &
   cmplx(.2d0,.1d0,8),cmplx(1d0,0d0,8)],[2,2]);ww=w0
 p0=reshape([cmplx(.4d0,.1d0,8),cmplx(-.2d0,.3d0,8)],[2,1]);wp=p0
 f0=w0*.3d0;face_self=f0;n0=w0*.1d0;face_neighbor=n0
 call apply_sawf_gauge_connection(q,q,basis,ww,wp,face_self,face_neighbor)
 call req(maxval(abs(basis-matmul(b0,q)))<1d-13,'basis')
 call req(maxval(abs(ww-matmul(conjg(transpose(q)),matmul(w0,q))))<1d-13,'WW')
 call req(maxval(abs(wp-matmul(conjg(transpose(q)),p0)))<1d-13,'WP')
 call req(maxval(abs(face_self-matmul(conjg(transpose(q)),matmul(f0,q))))<1d-13,'face self')
 call req(maxval(abs(face_neighbor-matmul(conjg(transpose(q)),matmul(n0,q))))<1d-13,'face neighbor')
 ! Atomic orchestration: a bad overlap must leave every payload unchanged.
 basis=b0;ww=w0;wp=p0;face_self=f0;face_neighbor=n0;basis_before=basis
 overlap=0;overlap(1,1)=1;neighbor_q=0;neighbor_q(1,1)=1;neighbor_q(2,2)=1
 call stitch_and_apply_sawf_neighbor_pair(overlap,1d-10,neighbor_q,basis,ww,wp, &
   face_self,face_neighbor,gauge,residual,ok,msg)
 call req(.not.ok.and.maxval(abs(basis-basis_before))<1d-15,'rank failure is atomic')
 ! A unitary overlap is accepted and applied.
 overlap=q;basis=b0;ww=w0;wp=p0;face_self=f0;face_neighbor=n0
 call stitch_and_apply_sawf_neighbor_pair(overlap,1d-10,neighbor_q,basis,ww,wp, &
   face_self,face_neighbor,gauge,residual,ok,msg)
 call req(ok.and.residual<1d-12,'valid stitched application')
 write(*,'(a)')'PASS gauge connection applied to basis, WW, WP, and DG face blocks'
contains
 subroutine req(x,t);logical,intent(in)::x;character(*),intent(in)::t
  if(.not.x)then;write(*,'(a)')trim(t);error stop 1;endif
 end subroutine
end program'''
with tempfile.TemporaryDirectory(prefix='sawf-gauge-apply-') as td:
 td=Path(td);(td/'driver.f90').write_text(driver)
 (td/'CMakeLists.txt').write_text(f'''cmake_minimum_required(VERSION 3.16)
project(gauge_apply LANGUAGES Fortran)
find_package(LAPACK REQUIRED)
add_executable(check "{ROOT/'src/gs/dc/lcfo_wannier_sawf_templates.f90'}" driver.f90)
target_link_libraries(check PRIVATE ${{LAPACK_LIBRARIES}})
''')
 env=dict(os.environ)
 for cmd in (["cmake","-S",str(td),"-B",str(td/'b'),f"-DCMAKE_Fortran_COMPILER={FC}"],
             ["cmake","--build",str(td/'b'),"-j","2"]):
  r=subprocess.run(cmd,env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
  if r.returncode:raise SystemExit(r.stdout)
 r=subprocess.run([str(td/'b/check')],env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
 if r.returncode:raise SystemExit(r.stdout)
 print(r.stdout.strip())
