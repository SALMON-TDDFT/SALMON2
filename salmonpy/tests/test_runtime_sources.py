"""Bundled source inventories reject corruption before preparation mutations."""
import hashlib,importlib.util,json,tempfile,unittest,os
from pathlib import Path
import salmonpy
package=Path(salmonpy.__file__).resolve().parent
root=package/'_runtime_source' if (package/'_runtime_source').is_dir() else package.parent
spec=importlib.util.spec_from_file_location('source_store',root/'tools/runtime_sources.py')
store=importlib.util.module_from_spec(spec);spec.loader.exec_module(store)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class RuntimeSources(unittest.TestCase):
 def test_bundled_trees_are_complete_and_byte_exact(self):
  profile=json.loads((root/'runtime/source.json').read_text())
  self.assertEqual(profile['source_tree'],'sources/solver')
  self.assertEqual(profile['source_profile'],'fast-pr')
  from salmonpy._runtime_capabilities import NLCC_SOURCE_TREE_SHA256
  self.assertEqual(NLCC_SOURCE_TREE_SHA256,profile['source_tree_sha256'],
                   'maintained solver changed: refresh NLCC evidence before advancing the capability identity')
  for name in ('hse-rt-state.json','hse-explicit-routing.json'):
   other=json.loads((root/'runtime'/name).read_text())
   self.assertEqual(other['source_tree'],profile['source_tree'])
   self.assertEqual(other['source_tree_sha256'],profile['source_tree_sha256'])
  files=store.manifest(root/'runtime',profile,'source_tree_manifest','source_tree_sha256','salmonpy.source-tree/1')
  store.check_tree(store.path_in(root/'runtime',profile['source_tree']),files)
  # A populated local checkout can hide ignored, never-committed source assets.
  if (root/'.git').exists():
   import subprocess
   prefix='runtime/'+profile['source_tree']+'/'
   tracked=set(subprocess.check_output(['git','ls-files','-z','--',prefix],cwd=root).decode().split('\0'))
   self.assertEqual(sorted(prefix+n for n in files if prefix+n not in tracked),[],
                    'manifest files must be tracked, not merely present locally')

  effective=store.inventory(root/'runtime')['profiles']['fast_hse']
  self.assertEqual(set(effective),set(files))
  for kernel,external in [('sbe_types.f90','src/ssbe/sbe_types.f90'),
                          ('sbe_velocity_kernel.f90','src/ssbe/sbe_velocity_kernel.f90'),
                          ('lorentz_drude.f90','src/maxwell/lorentz_drude.f90')]:
   self.assertEqual(effective[external]['sha256'],sha(root/'native'/kernel))
   self.assertEqual(effective[external]['source'],'sources/solver/'+external)
  worker=(root/'runtime/sources/solver/src/rt/tddft_worker.f90').read_text()
  self.assertEqual([n for n,line in enumerate(worker.splitlines(),1) if len(line)>132],[],
                   'worker source must compile with the standard free-form line limit')
  parser=(root/'runtime/sources/solver/src/io/inputoutput.f90').read_text()
  self.assertIn('!SALMONPY_END_INPUT',parser)
  self.assertIn('if(.not.worker_embedded) call worker_configure()',parser)
  for name in ('src/maxwell/embedded_eh.f90','src/io/model_inputs.f90','src/io/step_output.f90',
               'src/gs/gs_seed.f90','src/rt/embedded_tddft.f90'):
   self.assertEqual(effective[name]['sha256'],sha(root/'runtime/sources/solver'/name))
   self.assertEqual(effective[name]['source'],'sources/solver/'+name)
 def setUp(self):
  self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
  self.root=Path(self.tmp.name);self.tree=self.root/'tree';self.tree.mkdir();(self.tree/'x.f90').write_text('old\n')
  self.manifest=self.root/'tree.json';self.manifest.write_text(json.dumps({'schema':'salmonpy.source-tree/1','files':{'x.f90':sha(self.tree/'x.f90')}}))
  self.spec={'source_tree':'tree','source_tree_manifest':'tree.json','source_tree_sha256':sha(self.manifest)}
 def test_copy_matches_and_source_unchanged(self):
  out=self.root/'out';store.copy_tree(self.root,self.spec,out)
  self.assertEqual((out/'x.f90').read_bytes(),(self.tree/'x.f90').read_bytes())
  import tarfile,io
  archive=self.root/'history.tar.gz'
  with tarfile.open(archive,'w:gz') as stream:stream.add(self.tree/'x.f90',arcname='x.f90')
  archived={**self.spec,'source_tree':archive.name,'source_archive_sha256':sha(archive)}
  archive_out=self.root/'historical';store.copy_tree(self.root,archived,archive_out)
  self.assertEqual((archive_out/'x.f90').read_bytes(),(self.tree/'x.f90').read_bytes())
  self.assertEqual(store.effective_files(self.root,archived)['x.f90']['member'],'x.f90')
  # Reference payloads may contain only known regular files; reject links before writes.
  with tarfile.open(archive,'w:gz') as stream:
   member=tarfile.TarInfo('x.f90');member.type=tarfile.SYMTYPE;member.linkname='../outside';stream.addfile(member)
  archived['source_archive_sha256']=sha(archive)
  with self.assertRaisesRegex(ValueError,'links/special'):store.copy_tree(self.root,archived,self.root/'bad-archive')
  self.assertFalse((self.root/'bad-archive').exists())
  archived['source_archive_sha256']='0'*64
  with self.assertRaisesRegex(ValueError,'archive hash'):store.copy_tree(self.root,archived,self.root/'wrong-archive')

  store.check_shared_kernels(self.root,out,{'tree/x.f90':'x.f90'})
  (out/'x.f90').write_text('changed')
  with self.assertRaisesRegex(ValueError,'shared kernel'):
   store.check_shared_kernels(self.root,out,{'tree/x.f90':'x.f90'})
 def test_missing_changed_extra_and_manifest_hash_rejected(self):
  (self.tree/'x.f90').write_text('changed')
  with self.assertRaises(ValueError):store.copy_tree(self.root,self.spec,self.root/'out')
  self.assertFalse((self.root/'out').exists())
  (self.tree/'x.f90').write_text('old\n');(self.tree/'extra').write_text('extra')
  with self.assertRaises(ValueError):store.copy_tree(self.root,self.spec,self.root/'out')
  (self.tree/'extra').unlink();self.manifest.write_text('{}')
  with self.assertRaises(ValueError):store.copy_tree(self.root,self.spec,self.root/'out')
 def test_replacements_preflight_and_new_files(self):
  replacement=self.root/'replacement';replacement.write_text('new\n')
  files={'x.f90':{'before_sha256':sha(self.tree/'x.f90'),'after_sha256':sha(replacement),'source':'replacement'},'added.f90':{'before_sha256':None,'after_sha256':sha(replacement),'source':'replacement'}}
  manifest=self.root/'replace.json';manifest.write_text(json.dumps({'schema':'salmonpy.source-replacements/1','files':files}))
  spec={'replacements':'replace.json','replacements_sha256':sha(manifest)}
  resolved=store.effective_files(self.root,self.spec,overlays=(spec,))
  self.assertEqual(set(resolved),{'x.f90','added.f90'})
  self.assertEqual(resolved['x.f90']['sha256'],sha(replacement))
  self.assertEqual(resolved['x.f90']['source'],'replacement')
  self.assertEqual((self.tree/'x.f90').read_text(),'old\n')
  original=store.effective_files(self.root,self.spec)
  self.assertEqual(original['x.f90']['source'],'tree/x.f90')
  with self.assertRaises(ValueError):store.effective_files(self.root,self.spec,overlays=(spec,spec))
  self.assertFalse((self.root/'out').exists())
  (self.tree/'added.f90').write_text('keep')
  with self.assertRaises(ValueError):store.replace_files(self.root,spec,self.tree)
  self.assertEqual((self.tree/'x.f90').read_text(),'old\n')
  (self.tree/'added.f90').unlink();store.replace_files(self.root,spec,self.tree)
  self.assertEqual((self.tree/'x.f90').read_text(),'new\n')
  self.assertEqual((self.tree/'added.f90').read_text(),'new\n')
 def test_traversal_rejected(self):
  self.spec['source_tree']='../outside'
  with self.assertRaises(ValueError):store.copy_tree(self.root,self.spec,self.root/'out')
 def test_symlink_rejected(self):
  target=self.root/'external';target.write_text('old\n');(self.tree/'x.f90').unlink()
  try:(self.tree/'x.f90').symlink_to(target)
  except OSError:self.skipTest('symlink permission unavailable')
  with self.assertRaises(ValueError):store.copy_tree(self.root,self.spec,self.root/'out')

 @unittest.skipUnless(hasattr(os,'mkfifo'),'POSIX FIFO test')
 def test_fifo_rejected_before_copy(self):
  os.mkfifo(self.tree/'pipe')
  with self.assertRaisesRegex(ValueError,'special source'):
   store.copy_tree(self.root,self.spec,self.root/'out')
  self.assertFalse((self.root/'out').exists())

 def test_pip_bytecode_is_not_copied_or_used(self):
  py=self.tree/'helper.py';py.write_text('pass\n')
  files={'x.f90':sha(self.tree/'x.f90'),'helper.py':sha(py)}
  self.manifest.write_text(json.dumps({'schema':'salmonpy.source-tree/1','files':files}))
  self.spec['source_tree_sha256']=sha(self.manifest)
  cache=self.tree/'__pycache__';cache.mkdir();(cache/'helper.cpython-314.pyc').write_bytes(b'untrusted bytecode')
  store.copy_tree(self.root,self.spec,self.root/'out')
  self.assertFalse((self.root/'out/__pycache__').exists())
  (cache/'other.cpython-314.pyc').write_bytes(b'unlisted')
  with self.assertRaises(ValueError):store.copy_tree(self.root,self.spec,self.root/'another')
