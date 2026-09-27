"""Compiler-selection tests; fake executables are never used to compile code."""
from pathlib import Path
import os, subprocess, tempfile
repo=Path(__file__).resolve().parents[2]
with tempfile.TemporaryDirectory() as td:
 root=Path(td);binpath=root/'bin';binpath.mkdir()
 for name in ('mpifrtpx','mpifccpx'):
  p=binpath/name;p.write_text('#!/bin/sh\nexit 99\n');p.chmod(0o755)
 env={k:v for k,v in os.environ.items() if k not in ('CC','FC','CMAKE_TOOLCHAIN_FILE')}
 env['PATH']=str(binpath)+os.pathsep+env['PATH']
 def run(name,setup,expected,extra_env=None,fail=False):
  script=root/(name+'.cmake')
  script.write_text('cmake_minimum_required(VERSION 3.14)\nset(CMAKE_HOST_SYSTEM_NAME Linux)\n'+setup+'\ninclude("'+str(repo/'cmakefiles/select_platform.cmake')+'")\n'+expected)
  p=subprocess.run(['cmake','-P',str(script)],env={**env,**(extra_env or {})},text=True,capture_output=True)
  assert (p.returncode!=0) if fail else (p.returncode==0),(name,p.stdout,p.stderr)
  if fail:assert 'SALMON_PLATFORM' in p.stderr,(name,p.stderr)
 selected='if(NOT CMAKE_TOOLCHAIN_FILE MATCHES "/platforms/fugaku.cmake$")\nmessage(FATAL_ERROR "not selected")\nendif()\n'
 absent='if(DEFINED CMAKE_TOOLCHAIN_FILE)\nmessage(FATAL_ERROR "unexpected selection")\nendif()\n'
 run('auto','',selected)
 run('mac','set(CMAKE_HOST_SYSTEM_NAME Darwin)',absent)
 run('generic','set(SALMON_PLATFORM generic)',absent)
 run('explicit_fc','set(CMAKE_Fortran_COMPILER other-fc)',absent)
 run('explicit_cc','set(CMAKE_C_COMPILER other-cc)',absent)
 run('env_fc','',absent,{'FC':'other-fc'})
 run('env_cc','',absent,{'CC':'other-cc'})
 run('toolchain','set(CMAKE_TOOLCHAIN_FILE custom.cmake)','if(NOT CMAKE_TOOLCHAIN_FILE STREQUAL "custom.cmake")\nmessage(FATAL_ERROR "override")\nendif()')
 run('env_toolchain','',absent,{'CMAKE_TOOLCHAIN_FILE':'custom.cmake'})
 run('force','set(SALMON_PLATFORM fugaku)',selected)
 run('invalid','set(SALMON_PLATFORM bogus)','',fail=True)
 run('conflict','set(SALMON_PLATFORM fugaku)\nset(CMAKE_Fortran_COMPILER other-fc)','',fail=True)
 (binpath/'mpifccpx').unlink()
 run('missing','',absent)
 run('force_missing','set(SALMON_PLATFORM fugaku)','',fail=True)
 # Defaults and compatibility alias without invoking a compiler.
 for filename in ('fugaku.cmake','fujitsu-a64fx-ea.cmake'):
  for off in (False,True):
   script=root/'defaults.cmake'
   script.write_text('cmake_minimum_required(VERSION 3.14)\n'+('set(USE_MPI OFF CACHE BOOL "")\nset(USE_SCALAPACK OFF CACHE BOOL "")\n' if off else '')+'include("'+str(repo/'platforms'/filename)+'")\ninclude("'+str(repo/'cmakefiles/misc.cmake')+'")\noption_set(USE_MPI "MPI" OFF)\noption_set(USE_SCALAPACK "ScaLAPACK" OFF)\nif(NOT CMAKE_BUILD_TYPE STREQUAL "Release")\nmessage(FATAL_ERROR "no release")\nendif()\n'+('if(USE_MPI OR USE_SCALAPACK)' if off else 'if(NOT USE_MPI OR NOT USE_SCALAPACK)')+'\nmessage(FATAL_ERROR "bad defaults")\nendif()\n')
   subprocess.run(['cmake','-P',str(script)],env=env,check=True)
 # Re-inclusion after compiler detection must retain the chosen toolchain.
 (binpath/'mpifccpx').write_text('#!/bin/sh\nexit 99\n');(binpath/'mpifccpx').chmod(0o755)
 again='include("'+str(repo/'cmakefiles/select_platform.cmake')+'")\n'
 run('repeat','',selected+'set(CMAKE_Fortran_COMPILER mpifrtpx)\n'+again+selected)
 script=root/'debug.cmake'
 script.write_text('cmake_minimum_required(VERSION 3.14)\nset(CMAKE_BUILD_TYPE Debug CACHE STRING "")\ninclude("'+str(repo/'platforms/fugaku.cmake')+'")\nif(NOT CMAKE_BUILD_TYPE STREQUAL "Debug")\nmessage(FATAL_ERROR "overrode Debug")\nendif()\n')
 subprocess.run(['cmake','-P',str(script)],env=env,check=True)
 # An MPI-off override alone must not accidentally enable ScaLAPACK.
 script=root/'mpi-off.cmake'
 script.write_text('cmake_minimum_required(VERSION 3.14)\nset(USE_MPI OFF CACHE BOOL "")\ninclude("'+str(repo/'platforms/fugaku.cmake')+'")\nif(USE_SCALAPACK_DEFAULT)\nmessage(FATAL_ERROR "ScaLAPACK requires MPI")\nendif()\n')
 subprocess.run(['cmake','-P',str(script)],env=env,check=True)
 # Target toolchain must propagate to all HSE dependency projects.
 script=root/'external.cmake'
 toolchain=str(repo/'platforms/fugaku.cmake')
 script.write_text('cmake_minimum_required(VERSION 3.14)\nset(CMAKE_TOOLCHAIN_FILE "'+toolchain+'")\ninclude("'+toolchain+'")\ninclude("'+str(repo/'cmakefiles/Builder/hse_external_options.cmake')+'")\nlist(FIND SALMON_EXTERNAL_CMAKE_ARGS "-DCMAKE_TOOLCHAIN_FILE:STRING='+toolchain+'" found)\nif(found LESS 0)\nmessage(FATAL_ERROR "lost target toolchain")\nendif()\n')
 subprocess.run(['cmake','-P',str(script)],env=env,check=True)
 # Official configure.py entry point: capture its CMake invocation without
 # invoking unavailable Fujitsu compilers, then resolve the short toolchain name.
 capture=root/'configure-args.txt'
 fake=binpath/'cmake'
 fake.write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "'+str(capture)+'"\n')
 fake.chmod(0o755)
 import sys
 subprocess.run([sys.executable,str(repo/'configure.py'),'--arch=fujitsu-a64fx-ea',
                 '--enable-scalapack','--prefix='+str(root/'install')],cwd=root,env=env,check=True)
 args=capture.read_text().splitlines()
 assert 'CMAKE_TOOLCHAIN_FILE=fujitsu-a64fx-ea' in args,args
 assert 'USE_SCALAPACK=on' in args and 'CMAKE_BUILD_TYPE=Release' in args,args
 assert 'CMAKE_INSTALL_PREFIX='+str(root/'install') in args,args
 fake.unlink()
 script=root/'official-toolchain.cmake'
 script.write_text('cmake_minimum_required(VERSION 3.14)\nlist(APPEND CMAKE_MODULE_PATH "'+str(repo/'platforms')+'")\ninclude(fujitsu-a64fx-ea RESULT_VARIABLE resolved)\nif(NOT resolved MATCHES "fujitsu-a64fx-ea.cmake$")\nmessage(FATAL_ERROR "official arch not resolved")\nendif()\nif(NOT USE_MPI_DEFAULT OR NOT USE_SCALAPACK_DEFAULT)\nmessage(FATAL_ERROR "official arch defaults lost")\nendif()\n')
 subprocess.run(['cmake','-P',str(script)],env=env,check=True)
 # HSE cache assignments require Fortran allocatable assignment semantics.
 for filename in ('fugaku.cmake','fujitsu-a64fx-ea.cmake'):
  script=root/'alloc-assign.cmake'
  script.write_text('cmake_minimum_required(VERSION 3.14)\ninclude("'+str(repo/'platforms'/filename)+'")\n'+
   'foreach(mode DEBUG RELEASE)\n'+
   'if(NOT CMAKE_Fortran_FLAGS_${mode} MATCHES "(^| )-Nalloc_assign( |$)")\n'+
   'message(FATAL_ERROR "Missing allocatable assignment semantics")\nendif()\n'+
   'if(CMAKE_C_FLAGS_${mode} MATCHES "-Nalloc_assign")\n'+
   'message(FATAL_ERROR "Fortran flag leaked to C")\nendif()\nendforeach()\n')
  subprocess.run(['cmake','-P',str(script)],env=env,check=True)
print('Platform selection and toolchain default tests passed')
