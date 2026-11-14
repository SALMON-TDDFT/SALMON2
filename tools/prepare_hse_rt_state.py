"""Prepare a fresh pinned fast-PR HSE snapshot with optional RT-state overlay.

python tools/prepare_hse_rt_state.py LOCAL_PR_GIT NEW_OUTPUT [--reference-only]
No download/build/PR edits. Source audit precedes every preparation/build; local
CMake guards retain mandatory checks even though the PR excludes those guards.
Failed evidence is preserved and explicitly not approved for execution.
"""
import argparse,hashlib,json,os,re,shutil,subprocess,sys,tarfile,tempfile
from pathlib import Path,PurePosixPath


def prepare(source,directory,*,reference_only=False,explicit_routing=False):
    root=Path(__file__).resolve().parents[1];spec=json.loads((root/'runtime/hse-rt-state.json').read_text())
    patch=root/'runtime'/spec['patch'];checker=Path(__file__).with_name('check_forbidden_eigensolvers.py')
    if hashlib.sha256(patch.read_bytes()).hexdigest()!=spec['patch_sha256']:raise ValueError('HSE overlay hash changed')
    routing=None
    if explicit_routing:
        if reference_only:raise ValueError('explicit routing cannot use reference-only preparation')
        routing=json.loads((root/'runtime/hse-explicit-routing.json').read_text())
        if routing['parent_patch_sha256']!=spec['patch_sha256'] or routing['base_commit']!=spec['upstream_commit']:
            raise ValueError('explicit routing requires its pinned parent overlay')
        child=root/'runtime'/routing['patch']
        if hashlib.sha256(child.read_bytes()).hexdigest()!=routing['patch_sha256']:raise ValueError('routing overlay hash changed')
    source=Path(source).resolve();dest=Path(directory).resolve()
    if dest.is_relative_to(source) or dest.is_relative_to(root):raise ValueError('new external directory required')
    commit=spec['upstream_commit']
    if subprocess.check_output(['git','-C',str(source),'rev-parse',commit+'^{commit}'],text=True).strip()!=commit:raise ValueError('pinned PR commit unavailable')
    dest.mkdir(parents=True,exist_ok=False);marker=dest/'salmonpy-hse-rt-state.json'
    record={**spec,'status':'preparing','overlay_applied':not reference_only,'explicit_rt_routing':explicit_routing}
    marker.write_text(json.dumps(record,indent=2)+'\n')
    try:
        with tempfile.TemporaryFile() as raw:
            subprocess.run(['git','-C',str(source),'archive','--format=tar',commit],stdout=raw,check=True);raw.seek(0)
            with tarfile.open(fileobj=raw,mode='r:') as archive:
                for member in archive:
                    path=PurePosixPath(member.name)
                    if path.is_absolute() or '..' in path.parts or '\\' in member.name:raise ValueError('unsafe archive path')
                    target=dest.joinpath(*path.parts)
                    if member.isdir():target.mkdir(parents=True,exist_ok=True)
                    elif member.isfile():
                        target.parent.mkdir(parents=True,exist_ok=True)
                        with archive.extractfile(member) as inp,target.open('wb') as out:shutil.copyfileobj(inp,out)
                        target.chmod(member.mode & 0o777)
                    else:raise ValueError('unsupported archive member')
        subprocess.run([sys.executable,str(checker),'--source-root',str(dest/'src')],check=True)
        for name,expected in spec['reference_source_hashes'].items():
            if hashlib.sha256((dest/name).read_bytes()).hexdigest()!=expected:raise ValueError('pinned source changed: '+name)
        if not reference_only:
            env={k:v for k,v in os.environ.items() if k not in ('GIT_DIR','GIT_WORK_TREE','GIT_INDEX_FILE','GIT_COMMON_DIR')}
            env['GIT_CEILING_DIRECTORIES']=str(dest.parent)
            subprocess.run(['git','apply','--check',str(patch)],cwd=dest,env=env,check=True)
            subprocess.run(['git','apply','--whitespace=nowarn',str(patch)],cwd=dest,env=env,check=True)
            if routing is not None:
                subprocess.run(['git','apply','--check',str(child)],cwd=dest,env=env,check=True)
                subprocess.run(['git','apply','--whitespace=nowarn',str(child)],cwd=dest,env=env,check=True)
                record['routing_overlay']=routing
        # User-mandated build audits override the upstream PR's omission of them.
        (dest/'tools').mkdir(exist_ok=True);shutil.copy2(checker,dest/'tools/check_forbidden_eigensolvers.py')
        guard=dest/'cmakefiles/salmonpy_prohibition.cmake'
        guard.write_text('''find_package(Python3 REQUIRED COMPONENTS Interpreter)
set(SALMONPY_CHECKER "${CMAKE_SOURCE_DIR}/tools/check_forbidden_eigensolvers.py")
execute_process(COMMAND "${Python3_EXECUTABLE}" "${SALMONPY_CHECKER}" --source-root "${CMAKE_SOURCE_DIR}/src"
  RESULT_VARIABLE SALMONPY_AUDIT_STATUS)
if(NOT SALMONPY_AUDIT_STATUS EQUAL 0)
  message(FATAL_ERROR "Required eigensolver source audit failed")
endif()
add_custom_target(salmonpy_prohibition_sources ALL
  COMMAND "${Python3_EXECUTABLE}" "${SALMONPY_CHECKER}" --source-root "${CMAKE_SOURCE_DIR}/src" VERBATIM)
''')
        cmake=dest/'CMakeLists.txt';text=cmake.read_text()
        text,count=re.subn(r'(project\(SALMON\b[\s\S]*?\))',r'\1\ninclude(${CMAKE_SOURCE_DIR}/cmakefiles/salmonpy_prohibition.cmake)',text,count=1)
        if count!=1:raise ValueError('unsupported CMake project boundary')
        cmake.write_text(text)
        cmake=dest/'src/CMakeLists.txt';cmake.write_text(cmake.read_text()+'''
add_dependencies(${TARGET_NAME} salmonpy_prohibition_sources)
add_custom_command(TARGET ${TARGET_NAME} POST_BUILD
  COMMAND "${Python3_EXECUTABLE}" "${SALMONPY_CHECKER}" --object-tree "${CMAKE_BINARY_DIR}" VERBATIM)
''')
        subprocess.run([sys.executable,str(checker),'--source-root',str(dest/'src')],check=True)
        record.update(status='prepared-and-source-audited',checker_sha256=hashlib.sha256(checker.read_bytes()).hexdigest(),
            modified_source_hashes={name:hashlib.sha256((dest/name).read_bytes()).hexdigest() for name in spec['modified_sources'] if (dest/name).is_file()},
            cmake_guard_sha256=hashlib.sha256(guard.read_bytes()).hexdigest())
        marker.write_text(json.dumps(record,indent=2)+'\n');return dest
    except BaseException:
        record['status']='failed-not-approved';marker.write_text(json.dumps(record,indent=2)+'\n');raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('source');p.add_argument('output');p.add_argument('--reference-only',action='store_true');p.add_argument('--explicit-routing',action='store_true')
    a=p.parse_args();print(prepare(a.source,a.output,reference_only=a.reference_only,explicit_routing=a.explicit_routing))
