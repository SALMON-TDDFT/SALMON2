#!/usr/bin/env python3
#
#   Copyright 2017-2025 SALMON developers
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#

#
# Autotools (configure) like script with Python.
#
from optparse import OptionParser, OptionGroup
import os
import os.path
import platform
import subprocess

SOURCE_DIR = os.path.dirname(__file__)

def on_or_off(v) :
  if v:
    return 'on'
  else:
    return 'off'

def debug_or_release(v) :
  if v:
    return 'Debug'
  else:
    return 'Release'

def add_option(dic, name, var) :
  if var is not None:
    dic[name] = on_or_off(var)

def add_env(dic, name, var) :
  if var is not None:
    dic[name] = var

def detect_brew_prefix(machine):
  try:
    prefix = subprocess.check_output(['brew', '--prefix'], stderr=subprocess.DEVNULL)
    prefix = prefix.decode().strip()
    if prefix:
      return prefix
  except Exception:
    pass
  if machine == 'arm64':
    return '/opt/homebrew'
  if machine == 'x86_64':
    return '/usr/local'
  return None

def add_macos_openmp_defaults(defines, env_vars, user_env):
  if platform.system() != 'Darwin':
    return

  machine = platform.machine()
  brew_prefix = detect_brew_prefix(machine)
  if brew_prefix is None:
    raise RuntimeError('Unsupported macOS CPU architecture: {0}'.format(machine))

  llvm_prefix = os.path.join(brew_prefix, 'opt', 'llvm')
  libomp_prefix = os.path.join(brew_prefix, 'opt', 'libomp')
  cc = os.path.join(llvm_prefix, 'bin', 'clang')
  cxx = os.path.join(llvm_prefix, 'bin', 'clang++')
  libomp = os.path.join(libomp_prefix, 'lib', 'libomp.dylib')

  if not os.path.exists(cc) or not os.path.exists(cxx):
    raise RuntimeError('Homebrew LLVM clang/clang++ not found under {0}'.format(llvm_prefix))
  if not os.path.exists(libomp):
    raise RuntimeError('Homebrew libomp not found under {0}'.format(libomp_prefix))

  if user_env.get('CC') == '/usr/bin/clang' or user_env.get('CXX') == '/usr/bin/clang++':
    raise RuntimeError('Do not use Apple /usr/bin/clang for macOS OpenMP builds')

  env_vars.setdefault('CC', cc)
  env_vars.setdefault('CXX', cxx)
  env_vars.setdefault('OMPI_CC', cc)
  env_vars.setdefault('OMPI_CXX', cxx)
  env_vars['CFLAGS'] = append_env_flags(user_env.get('CFLAGS'), '-fopenmp', env_vars.get('CFLAGS'))
  env_vars['CPPFLAGS'] = append_env_flags(user_env.get('CPPFLAGS'), '-I{0}/include'.format(libomp_prefix),
                                         env_vars.get('CPPFLAGS'))
  env_vars['LDFLAGS'] = append_env_flags(user_env.get('LDFLAGS'), '-L{0}/lib'.format(libomp_prefix),
                                        env_vars.get('LDFLAGS'))
  env_vars['LIBS'] = append_env_flags(user_env.get('LIBS'), '-lomp', env_vars.get('LIBS'))

  defines.setdefault('CMAKE_C_COMPILER', cc)
  defines.setdefault('CMAKE_CXX_COMPILER', cxx)
  defines['OpenMP_ROOT'] = libomp_prefix
  defines['OpenMP_C_FLAGS'] = '-fopenmp'
  defines['OpenMP_C_LIB_NAMES'] = 'omp'
  defines['OpenMP_omp_LIBRARY'] = libomp

def append_env_flags(user_value, extra, current_value=None):
  if current_value:
    base = current_value
  elif user_value:
    base = user_value
  else:
    base = ''
  if base:
    return '{0} {1}'.format(base, extra)
  return extra

def remove_failed_cmake_cache():
  cache_path = os.path.join(os.getcwd(), 'CMakeCache.txt')
  if os.path.exists(cache_path):
    os.remove(cache_path)

usage  = "usage: %prog [options]"
parser = OptionParser(usage)

parser.add_option('-n', '--dry-run', action='store_true', default=False, dest='dry_run', help='don\'t actually run.')
parser.add_option('-v', '--verbose', action='store_true', default=False, dest='verbose', help='show verbose messages.')

group = OptionGroup(parser, 'Build target')
group.add_option('-a', '--arch',   action='store', default=None, dest='arch',   help='cross compile mode. ARCH format should be <COMPILER>-<SYSTEM>')
group.add_option('--prefix',       action='store', default=None, dest='prefix', help='package install prefix.')
parser.add_option_group(group)

group = OptionGroup(parser, 'Library options')
group.add_option('--enable-mpi',        action='store_true',  dest='mpi')
group.add_option('--disable-mpi',       action='store_false', dest='mpi',       help='enable/disable MPI parallelization.')
group.add_option('--enable-scalapack',  action='store_true',  dest='scalapack')
group.add_option('--disable-scalapack', action='store_false', dest='scalapack', help='enable/disable computations with ScaLAPACK library.')
group.add_option('--enable-eigenexa',   action='store_true',  dest='eigenexa')
group.add_option('--disable-eigenexa',  action='store_false', dest='eigenexa', help='enable/disable computations with EigenExa library (SALMON will build it internally)')
group.add_option('--enable-libxc',      action='store_true',  dest='libxc')
group.add_option('--disable-libxc',     action='store_false', dest='libxc', help='enable/disable Libxc library.')
group.add_option('--enable-fftw',        action='store_true',  dest='fftw')
group.add_option('--disable-fftw',       action='store_false', dest='fftw', help='enable/disable computations FFTW library.')
group.add_option('--enable-wannier90',   action='store_true',  dest='wannier90')
group.add_option('--disable-wannier90',  action='store_false', dest='wannier90', help='enable/disable computations with Wannier90 library (SALMON will build it internally).')
group.add_option('--with-lapack',       action='store', type=str, default=None, dest='lapack_installdir', help='specify install path to LAPACK/ScaLAPACK')
group.add_option('--with-libxc',        action='store', type=str, default=None, dest='libxc_installdir', help='specify install path to Libxc package')

parser.add_option_group(group)

group = OptionGroup(parser, 'Optimization options for stencil computations')
group.add_option('--explicit-vec',          action='store_true',  dest='explicit_vec',      help='enable explicit vectorization. it requires --simd-set option to be set.')
group.add_option('--compiler-vec',          action='store_false', dest='explicit_vec',      help='entrust optimization to a compiler.')
group.add_option('--simd-set',              action='store',       dest='simd',              help='specifies SIMD instruction set. (e.g. AVX, AVX_512, HPC_ACE2...)')
group.add_option('--enable-array-padding',  action='store_true',  dest='padding')
group.add_option('--disable-array-padding', action='store_false', dest='padding',           help='enable/disable array padding for the cache utilization.')
parser.add_option_group(group)

group = OptionGroup(parser, 'Debug options')
group.add_option('-d', '--debug',   action='store_true',  default=False, dest='debug', help='enable debug build.')
group.add_option('-r', '--release', action='store_false',                dest='debug', help='enable release build.')
parser.add_option_group(group)


(options, args) = parser.parse_args()

### check options
if (options.libxc_installdir is not None):
  options.libxc = True

dict = {}
env_dict = {}
user_env_dict = {}
for var in args:
  (k, v) = var.split('=', 1)
  user_env_dict[k] = v
  env_dict[k] = v

if options.arch is not None:
  arch = options.arch.lower()
  toolchain = arch
  if not os.path.isabs(toolchain):
    if not toolchain.endswith('.cmake'):
      toolchain = toolchain + '.cmake'
    candidate = os.path.join(SOURCE_DIR, 'platforms', toolchain)
    if os.path.exists(candidate):
      toolchain = candidate
  dict['CMAKE_TOOLCHAIN_FILE']     = toolchain
if options.prefix is not None:
  dict['CMAKE_INSTALL_PREFIX']     = options.prefix
dict['CMAKE_BUILD_TYPE']           = debug_or_release(options.debug)
dict['CMAKE_VERBOSE_MAKEFILE']     = on_or_off(options.verbose)

add_env(dict, 'LAPACK_INSTALLDIR', options.lapack_installdir)
add_env(dict, 'LIBXC_INSTALLDIR',  options.libxc_installdir)

add_option(dict, 'USE_MPI',             options.mpi)
add_option(dict, 'USE_SCALAPACK',       options.scalapack)
add_option(dict, 'USE_EIGENEXA',        options.eigenexa)
add_option(dict, 'USE_LIBXC',           options.libxc)
add_option(dict, 'USE_FFTW',            options.fftw)
add_option(dict, 'USE_WANNIER90',       options.wannier90)

add_option(dict, 'USE_OPT_ARRAY_PADDING',          options.padding)
add_option(dict, 'USE_OPT_EXPLICIT_VECTORIZATION', options.explicit_vec)

if options.simd is not None:
  dict['SIMD_SET'] = options.simd.upper()

add_macos_openmp_defaults(dict, env_dict, user_env_dict)
if not options.dry_run:
  remove_failed_cmake_cache()

define = ''
for k,v in dict.items():
  define = '{0} -D {1}={2}'.format(define, k, v)

env = ''
for k, v in env_dict.items():
  env = '{0} {1}="{2}"'.format(env, k, v)

### configuration
comm = '{2} cmake {0} {1}'.format(define, SOURCE_DIR, env)
print('    $ %s' % comm)
if not options.dry_run:
  os.system(comm)
