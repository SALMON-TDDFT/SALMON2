"""Synthetic retained GS export under actual MPI; never loads a SALMON runtime/GS.

mpiexec -n 4 python tools/test_retained_gs_archive_mpi.py NEW_OUTPUT
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from types import SimpleNamespace
from unittest.mock import patch,PropertyMock
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from salmonpy.dft import DFT,Pseudopotential
from salmonpy.native_gs import NativeGSState,NativeGSArray,NativeGSInfo
from salmonpy.gs_archive import GroundStateArchive
from salmonpy.errors import InputError


def fixture(comm,decomposition):
    nk,nb=4,5
    grid=(8,10,12)
    model=DFT((2,),np.diag([8.,10.,12.]),[[.1,.2,.3]],
        (Pseudopotential(2,0,'_rps.dat',b'synthetic-only'),),2,nb,grid,(2,2,1),xc='libxc_pbe')
    rng=np.random.default_rng(1119)
    waves=rng.normal(size=(nk,nb,12,10,8))+1j*rng.normal(size=(nk,nb,12,10,8))
    occupations=np.zeros((nk,nb));occupations[:,0]=2
    eigenvalues=np.arange(nk*nb,dtype=float).reshape(nk,nb)/7
    rank=comm.rank
    coords=[]
    for parts in decomposition:
        coords.append(rank%parts);rank//=parts
    limits=[]
    for n,c,parts in zip((nk,nb,*grid),coords,decomposition):
        limits.append((n*c//parts,n*(c+1)//parts))
    (ks,ke),(bs,be),(xs,xe),(ys,ye),(zs,ze)=limits
    tile=waves[ks:ke,bs:be,zs:ze,ys:ye,xs:xe].transpose(4,3,2,1,0)
    raw=np.full((xe-xs+2,ye-ys+4,ze-zs+6,1,be-bs,ke-ks,1),complex(999,999))
    raw[1:-1,2:-2,3:-3,0,:,:,0]=tile
    values={
        'density':NativeGSArray(np.zeros((xe-xs,ye-ys,ze-zs)),(xs+1,ys+1,zs+1),'density','bohr^-3'),
        'orbitals':NativeGSArray(raw,(xs,ys-1,zs-2,1,bs+1,ks+1,1),'orbitals','native'),
        'occupations':NativeGSArray(occupations.T[:,:,None].copy(),(1,1,1),'occupations','1'),
        'eigenvalues':NativeGSArray(eigenvalues.T[:,:,None].copy(),(1,1,1),'eigenvalues','hartree')}
    session=SimpleNamespace(_live=lambda:None,comm=comm,
        backend=SimpleNamespace(decomposition=decomposition,sha256='c'*64))
    state=NativeGSState(session,1,None,'synthetic',model);state._identity='synthetic-shared-owner'
    state.array=lambda name:values[name]
    return state,values,waves


def hashes(root):
    return {str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in root.rglob('*') if p.is_file()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output',type=Path);args=parser.parse_args()
    from mpi4py import MPI
    comm=MPI.COMM_WORLD
    if comm.size!=4:raise ValueError('use four MPI ranks')
    error=None
    if comm.rank==0:
        try:args.output.mkdir(parents=True,exist_ok=False)
        except Exception as exc:error=str(exc)
    error=comm.bcast(error,root=0)
    if error:raise RuntimeError(error)
    info=NativeGSInfo(-2.,1e-12,1e-9,10,True,1,5)
    cases={'k':(4,1,1,1,1),'orbital':(1,4,1,1,1),'spatial':(1,1,2,2,1),
           'k-orbital':(2,2,1,1,1),'k-spatial':(2,1,2,1,1)}
    reports=[]
    with patch.object(NativeGSState,'info',new_callable=PropertyMock,return_value=info):
        reference=None
        if comm.rank==0:
            single=SimpleNamespace(size=1,rank=0,allgather=lambda v:[v])
            state,_,_=fixture(single,(1,1,1,1,1))
            reference=hashes(state.export_archive(args.output/'serial').directory)
        for label,dec in cases.items():
            state,values,waves=fixture(comm,dec)
            before=hashlib.sha256(values['orbitals'].values.tobytes()).hexdigest()
            archive=state.export_archive(args.output/label,orbital_block=16)
            assert (archive is not None)==(comm.rank==0)
            assert hashlib.sha256(values['orbitals'].values.tobytes()).hexdigest()==before
            if comm.rank==0:
                assert hashes(archive.directory)==reference,label
                for k in range(4):np.testing.assert_array_equal(archive.wavefunctions(k),waves[k])
            reports.append({'case':label,'status':'pass','serial_files_exact':True})
        dec=cases['k-spatial']
        for failure in ('bounds','replicas','overlap','gap','metadata','options','allocation','writer'):
            state,values,waves=fixture(comm,dec)
            before=hashlib.sha256(values['orbitals'].values.tobytes()).hexdigest()
            if failure=='bounds' and comm.rank==3:
                old=values['orbitals'];lower=list(old.lower_indices);lower[5]=99
                values['orbitals']=NativeGSArray(old.values,tuple(lower),old.quantity,old.unit)
            elif failure=='replicas' and comm.rank==3:
                values['occupations'].values[0,0,0]=1.9
            elif failure=='overlap':
                # Every rank supplies rank0's tile: in-range but duplicated/incomplete.
                other=SimpleNamespace(size=4,rank=0,allgather=comm.allgather)
                _,duplicate,_=fixture(other,dec);values.update(duplicate)
                before=hashlib.sha256(values['orbitals'].values.tobytes()).hexdigest()
            if failure=='gap' and comm.rank==3:
                old=values['orbitals']
                values['orbitals']=NativeGSArray(old.values[:,:,:,:,:-1,:,:],old.lower_indices,old.quantity,old.unit)
            before_all={name:hashlib.sha256(value.values.tobytes()).hexdigest() for name,value in values.items()}
            destination=args.output/('failed-'+failure)
            def fail(*a,**kw):raise OSError('injected archive write failure')
            from contextlib import nullcontext
            context=nullcontext()
            if failure=='writer':context=patch('salmonpy.gs_archive._write',side_effect=fail)
            if failure=='allocation':context=patch('salmonpy._gs_archive_mpi.np.empty',side_effect=MemoryError('injected root allocation failure'))
            chosen_info=NativeGSInfo(-2.,1e-12,1e-9,10,True,2,5) if failure=='metadata' and comm.rank==3 else info
            caught=None
            try:
                with context,patch.object(NativeGSState,'info',new_callable=PropertyMock,return_value=chosen_info):
                    state.export_archive(destination,orbital_block=2 if failure=='options' and comm.rank==3 else 16)
            except (InputError,OSError,MemoryError) as exc:caught=(type(exc).__name__,str(exc))
            errors=comm.allgather(caught)
            assert all(v is not None for v in errors),(failure,errors)
            assert all(v==errors[0] for v in errors),(failure,errors)
            assert not state.closed
            assert {name:hashlib.sha256(value.values.tobytes()).hexdigest() for name,value in values.items()}==before_all
            if comm.rank==0:assert not destination.exists()
            clean,_,_=fixture(comm,dec)
            retried=clean.export_archive(destination)
            if comm.rank==0:assert hashes(retried.directory)==reference
            reports.append({'case':'reject-'+failure,'status':'pass','retry_exact':True})
        state,_,_=fixture(comm,dec)
        caught=None
        try:state.export_archive(args.output/'serial')
        except FileExistsError as exc:caught=type(exc).__name__
        assert comm.allgather(caught)==['FileExistsError']*comm.size
        if comm.rank==0:assert hashes(args.output/'serial')==reference
        reports.append({'case':'existing-output-preserved','status':'pass'})
    if comm.rank==0:
        result={'status':'pass','scope':'synthetic real-MPI export only; no SALMON/SCF/GS input',
                'ranks':comm.size,'cases':reports}
        (args.output/'result.json').write_text(json.dumps(result,indent=2)+'\n')
        print(json.dumps(result))


if __name__=='__main__':main()
