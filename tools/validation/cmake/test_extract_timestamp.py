"""Exercise the project's policy preamble with an offline ExternalProject archive."""
from pathlib import Path
import io, subprocess, tarfile, tempfile
repo=Path(__file__).resolve().parents[3]
with tempfile.TemporaryDirectory() as td:
 root=Path(td);source=root/'source';source.mkdir();build=root/'build'
 archive=root/'payload.tar.gz'
 with tarfile.open(archive,'w:gz') as tf:
  info=tarfile.TarInfo('payload/check.txt');info.size=2;info.mtime=946684800
  tf.addfile(info,io.BytesIO(b'ok'))
 preamble=(repo/'CMakeLists.txt').read_text().split('if ("${CMAKE_CURRENT_SOURCE_DIR}"',1)[0]
 (source/'CMakeLists.txt').write_text(preamble+'''\nproject(timestamp_probe NONE)
include(ExternalProject)
ExternalProject_Add(payload URL "'''+str(archive)+'''"
 SOURCE_DIR "${CMAKE_BINARY_DIR}/unpacked"
 CONFIGURE_COMMAND "" BUILD_COMMAND "" INSTALL_COMMAND "")
''')
 subprocess.run(['cmake','-Werror=dev','-S',str(source),'-B',str(build)],check=True)
 subprocess.run(['cmake','--build',str(build)],check=True)
 assert (build/'unpacked/check.txt').read_bytes()==b'ok'
 version=subprocess.check_output(['cmake','--version'],text=True).split()[2]
 if tuple(map(int,version.split('.')[:2]))>=(3,24):
  assert (build/'unpacked/check.txt').stat().st_mtime>946684800,'archive time retained'
print('ExternalProject timestamp policy regression passed')
