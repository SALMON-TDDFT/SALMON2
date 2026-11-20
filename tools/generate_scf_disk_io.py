"""Generate explicit Fortran SCF stream IO from owned numerical type declarations.

Developer-only, no compilation. Fail on unknown syntax or unapproved pointers.
The native format is versioned and binary/topology-pinned, not a portable GS archive.
"""
import argparse,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]/'runtime/sources/solver/src'
FILES=['common/structures.f90','io/model_inputs.f90','gs/main_dft.f90','xc/exx_rt_state.f90','xc/exx_ace.f90','xc/exx_symmetry.f90','xc/exx_factor_history.f90','xc/exx_spatial.f90','xc/exx_distributed_gauge.f90','xc/exx_sparse_orbitals.f90']
SKIP={('s_dft_system','salmonpy_hse_owner'),('s_exx_rt_state','native_kernel'),('s_exx_ace','metric_comm'),('spatial_exx_state','local_batch'),('spatial_exx_state','create_local_backend'),('s_exx_gauge','comm')}

def split(text):
    result=[];start=depth=0;quote=None
    for i,c in enumerate(text):
        if quote:
            if c==quote:quote=None
        elif c in "'\"":quote=c
        elif c in '([':depth+=1
        elif c in ')]':depth-=1
        elif c==',' and depth==0:result.append(text[start:i].strip());start=i+1
    result.append(text[start:].strip());return result

def types():
    result={}
    for file in FILES:
        text=re.sub(r'!.*','',(ROOT/file).read_text());text=re.sub(r'&\s*\n\s*&?','',text)
        active=None
        for line in text.splitlines():
            s=line.strip();m=re.match(r'type\s+(\w+)\s*$|type\s*,\s*public\s*::\s*(\w+)',s,re.I)
            if m:active=next(v for v in m.groups() if v);result[active]=[];continue
            if active:
                if re.match(r'end\s*type|contains',s,re.I):active=None;continue
                if s.startswith('#'):result[active].append(s)
                elif '::' in s:
                    left,right=s.split('::',1)
                    for item in split(right):
                        item=item.split('=',1)[0].strip();m=re.fullmatch(r'(\w+)(?:\((.*)\))?',item)
                        if not m:raise ValueError((active,item))
                        name,shape=m.groups()
                        if (active,name) in SKIP:continue
                        if 'pointer' in left:
                            result[active].append(('POINTER',name,False,0));continue
                        sub=re.match(r'type\((\w+)\)',left,re.I)
                        result[active].append((name,sub[1] if sub else None,bool('allocatable' in left),len(split(shape)) if shape else 0))
    return result

def generate():
    data=types();seen=set()
    def visit(t):
        if t in seen:return
        seen.add(t)
        for f in data[t]:
            if not isinstance(f,str):
                if f[0]=='POINTER':raise ValueError(('unapproved pointer',t,f[1]))
                if f[1]:visit(f[1])
    visit('s_scf_payload')
    def routine(t):
        body=[f'subroutine scf_io_{t}(u,reading,value,ios,limit)', ' integer,intent(in)::u', ' logical,intent(in)::reading',
              f' class({t}),intent(inout)::value' if t=='s_model_inputs' else f' type({t}),intent(inout)::value',
              ' integer,intent(inout)::ios',' integer(8),intent(in)::limit',
              ' integer::lo(7),hi(7),k,j1,j2,j3,j4,j5,j6,j7,istat',' integer(8)::n,e', ' logical::present_array',' if(ios/=0)return']
        def io(expr):return [' if(reading)then',f'  read(u,iostat=ios){expr}',' else',f'  write(u,iostat=ios){expr}',' endif',' if(ios/=0)return']
        for f in data[t]:
            if isinstance(f,str):body.append(f);continue
            name,sub,alloc,rank=f;v='value%'+name
            if alloc:
                body += [f' present_array=allocated({v})']+io('present_array')
                body += [' if(reading)then',f'  if(allocated({v}))deallocate({v})',' endif',' if(present_array)then',
                         '  if(.not.reading)then',f'   lo(1:{rank})=lbound({v});hi(1:{rank})=ubound({v})','  endif']
                body += io(f'lo(1:{rank}),hi(1:{rank})')
                body += ['  if(reading)then','   n=1_8',f'   do k=1,{rank}','    e=int(hi(k),8)-int(lo(k),8)+1_8',
                         '    if(e<0_8.or.e>limit)then','     ios=1;return','    endif',
                         '    if(e>0_8)then','     if(n>limit/e)then','      ios=1;return','     endif','    endif','    n=n*e','   enddo']
                dims=','.join(f'lo({i}):hi({i})' for i in range(1,rank+1))
                body += [f'   allocate({v}({dims}),stat=istat)','   if(istat/=0)then','    ios=1;return','   endif','  endif']
            if sub:
                for j in range(rank,0,-1):body += [f' do j{j}=lbound({v},{j}),ubound({v},{j})']
                expr=v+('('+','.join(f'j{j}' for j in range(1,rank+1))+')' if rank else '')
                body += [f' call scf_io_{sub}(u,reading,{expr},ios,limit)',' if(ios/=0)return']
                body += [' enddo']*rank
            else:body+=io(v)
            if alloc:body+=[' endif']
        body+=['end subroutine'];return '\n'.join(body)+'\n'
    outputs={ROOT/'io/scf_model_disk_io.inc':routine('s_model_inputs')}
    content='! Generated by tools/generate_scf_disk_io.py; numerical fields only.\n'
    hsetypes={'s_exx_rt_state','s_exx_ace','s_exx_symmetry_map','s_factor_history','spatial_exx_state','s_exx_gauge','s_sparse_orbitals'}
    for t in sorted(seen-{'s_model_inputs'}):
        if t in hsetypes:content+='#ifdef USE_HSE\n'
        if t=='s_scf_payload':
            # Private model components are serialized by their owning module.
            r=routine(t).replace('call scf_io_s_model_inputs(u,reading,value%continuation_inputs,ios,limit)',
                                  'call value%continuation_inputs%disk_io(u,reading,ios,limit)')
        else:r=routine(t)
        content+=r
        if t in hsetypes:content+='#endif\n'
    outputs[ROOT/'gs/scf_disk_io.inc']=content
    return outputs

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--check',action='store_true');args=p.parse_args()
    import json
    catalog=Path(__file__).resolve().parents[1]/'salmonpy/scf_compatibility.json'
    if catalog.exists():
        expected=['src/'+x for x in FILES]+['src/gs/scf_disk_io.inc','src/io/scf_model_disk_io.inc']
        data=json.loads(catalog.read_text())
        if args.check:
            if data.get('codec_sources')!=expected:raise ValueError('SCF codec inventory is stale')
        else:
            data['codec_sources']=expected;catalog.write_text(json.dumps(data,indent=2)+'\n')
    for path,text in generate().items():
        if args.check:
            if path.read_text()!=text:raise ValueError('regenerate '+str(path))
        else:path.write_text(text)
if __name__=='__main__':main()
