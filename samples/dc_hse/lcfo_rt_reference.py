"""Strict Gamma/unpolarized complex LCFO reader for RT validation.

Reads the SALMON V1 wire format. Reconstructs the *same symmetrized* Hamiltonian
as LCFO diagonalization and separately exposes the directed-block discrepancy.
This is not a self-consistent TDHSE driver or a dielectric-response calculator.
"""
import argparse
import json
import math
from pathlib import Path
import numpy as np


def require(ok, message):
    if not ok:
        raise ValueError(message)


class Wire:
    def __init__(self, path):
        self.data = Path(path).read_bytes()
        self.pos = 0

    def raw(self, n):
        require(0 <= n <= len(self.data)-self.pos, 'truncated/oversized LCFO payload')
        out = self.data[self.pos:self.pos+n]
        self.pos += n
        return out

    def array(self, dtype, count):
        require(count >= 0, 'negative LCFO element count')
        return np.frombuffer(self.raw(np.dtype(dtype).itemsize*count), dtype=dtype).copy()

    def integers(self, n):
        return self.array('<i4', n)

    def complex(self, rows, cols):
        out = self.array('<c16', int(rows)*int(cols)).reshape((rows, cols), order='F')
        require(np.isfinite(out).all(), 'nonfinite LCFO matrix')
        return out


def read_file(path, kind):
    w = Wire(path)
    require(w.raw(16).rstrip() == b'SLCFO_COMPLEX_V1', 'invalid LCFO magic')
    require(np.array_equal(w.integers(6), [1, 0x01020304, 32, 64, kind, 0]), 'unsupported LCFO format')
    header_size = int(w.array('<i8', 1)[0])
    run = w.raw(96)
    meta = w.integers(20)
    geom = w.array('<f8', 10)
    # Narrow explicit initial scope; do not accidentally treat multi-k as Gamma.
    require(meta[16] == 1 and meta[17] == 1, 'RT reference currently requires unpolarized one-k LCFO')
    require(np.all(meta > 0), 'invalid LCFO dimensions')
    require(np.array_equal(meta[:3], meta[6:9]*meta[9:12]), 'inconsistent domain tiling')
    require(np.all(meta[3:6] >= meta[6:9]), 'invalid fragment dimensions')
    require(np.all(meta[12:15] <= meta[:3]), 'invalid fragment origin')
    require(np.isfinite(geom).all() and geom[9] > 0, 'invalid geometry')
    cell = geom[:9].reshape((3, 3), order='F')
    require(np.isclose(abs(np.linalg.det(cell))/math.prod(map(int,meta[:3])),geom[9],rtol=1e-10),
            'cell/grid volume mismatch')
    k = w.array('<f8', 3); weight = w.array('<f8', 1)
    require(np.max(abs(k)) < 1e-14 and abs(weight[0]-1) < 1e-14, 'Gamma with unit k weight required')
    require(w.pos == header_size == 336, 'LCFO header size mismatch')
    require(w.integers(1)[0] == 1, 'unexpected k record')
    payload = int(w.array('<i8', 1)[0]); end = w.pos+payload
    require(0 <= payload <= len(w.data)-w.pos-124, 'invalid record size')
    spin, nb = map(int, w.integers(2))
    require(spin == 1 and 0 <= nb <= meta[18], 'invalid basis dimensions')
    result = dict(run=run, meta=meta, geom=geom, nb=nb)
    if kind == 1:
        result['basis'] = w.complex(math.prod(map(int, meta[6:9])), nb)
    else:
        nmat = int(w.integers(1)[0]); nf = math.prod(map(int, meta[9:12]))
        require(1 <= meta[15] <= nf, 'invalid fragment id')
        counts = w.integers(nf)
        require(np.all(counts >= 0) and np.all(counts <= meta[18]), 'invalid basis counts')
        require(sum(map(int,counts)) == nmat and nmat >= meta[19], 'invalid global dimension')
        require(counts[meta[15]-1] == nb, 'inconsistent local basis count')
        result.update(counts=counts, nmat=nmat)
        if kind == 2:
            rows = w.integers(nb)
            start = sum(map(int, counts[:meta[15]-1]))
            require(np.array_equal(rows, np.arange(start+1, start+nb+1)), 'invalid coefficient row indices')
            result['coefficients'] = w.complex(nb, int(meta[19]))
        elif kind == 3:
            result['diagonal'] = w.complex(nb, nb)
            nh = int(w.integers(1)[0]);require(0 <= nh <= 26, 'invalid halo count')
            halos = []; seen = set()
            for _ in range(nh):
                src, dx, dy, dz, ns = map(int, w.integers(5))
                require(1 <= src <= nf and ns == counts[src-1], 'invalid halo source')
                direction=(dx,dy,dz)
                require(max(map(abs,direction)) == 1, 'invalid halo direction')
                require((src,direction) not in seen, 'duplicate halo')
                seen.add((src,direction))
                halos.append((src, direction, w.complex(ns,nb)))
            result['halos'] = halos
    require(w.pos == end, 'LCFO record byte count mismatch')
    require(w.raw(16).rstrip() == b'SLCFO_DONE_V1', 'missing completion footer')
    require(w.raw(96) == run, 'header/footer run mismatch')
    require(w.integers(1)[0] == 1, 'incomplete k records')
    require(w.array('<i8',1)[0] == len(w.data) and w.pos == len(w.data), 'LCFO file size mismatch')
    return result


def load_lcfo(fragments):
    root=Path(fragments)
    first=read_file(root/'000001/hamiltonian_local.bin',3)
    nf=math.prod(map(int, first['meta'][9:12])); n=first['nmat']
    counts=first['counts']; offsets=np.r_[0,np.cumsum(counts)]
    h=np.zeros((n,n),complex); directed=np.zeros_like(h)
    coeff=np.zeros((n,int(first['meta'][19])),complex)
    bases=[]; origins=set(); origin_by_fragment={}; halo_geometry=[]; basis_error=0.
    for frag in range(1,nf+1):
        folder=root/f'{frag:06d}'
        records=[read_file(folder/name,kind) for name,kind in
                 [('basis_functions.bin',1),('wavefunctions.bin',2),('hamiltonian_local.bin',3)]]
        for rec in records:
            require(rec['run']==first['run'],'mixed LCFO runs')
            # Only fragment index and origin may differ from the first fragment.
            common=np.r_[0:12,16:20]
            require(np.array_equal(rec['meta'][common],first['meta'][common]),'inconsistent LCFO metadata')
            require(np.array_equal(rec['geom'],first['geom']),'inconsistent geometry')
            require(rec['meta'][15]==frag and rec['nb']==counts[frag-1],'inconsistent fragment index/count')
            if 'counts' in rec:require(np.array_equal(rec['counts'],counts),'inconsistent basis layout')
        b,c,hh=records
        require(all(np.array_equal(rec['meta'],b['meta']) for rec in records),'inconsistent local metadata')
        origin=tuple(map(int,b['meta'][12:15]-1));core=b['meta'][6:9]
        require(all(origin[d] % core[d] == 0 for d in range(3)) and origin not in origins,'invalid core tiling')
        origins.add(origin);origin_by_fragment[frag]=np.array(origin)
        basis=b['basis']; nb=basis.shape[1]
        if nb:
            basis_error=max(basis_error,float(np.max(abs(basis.conj().T@basis*first['geom'][9]-np.eye(nb)))))
        require(basis_error < 1e-10,'nonorthonormal LCFO basis')
        bases.append(basis)
        sl=slice(offsets[frag-1],offsets[frag])
        coeff[sl]=c['coefficients'];h[sl,sl]=hh['diagonal'];directed[sl,sl]=hh['diagonal']
        for src,direction,block in hh['halos']:
            require(all(direction[d]==0 for d in range(3) if b['meta'][9+d]==1),
                    'halo points along unpartitioned axis')
            expected=(np.array(origin)-np.array(direction)*core) % b['meta'][:3]
            halo_geometry.append((src,expected))
            sr=slice(offsets[src-1],offsets[src])
            h[sr,sl]+=0.5*block;h[sl,sr]+=0.5*block.conj().T
            directed[sr,sl]+=block
    for src,expected in halo_geometry:
        require(np.array_equal(origin_by_fragment[src],expected),'halo source geometry mismatch')
    herm=float(np.max(abs(h-h.conj().T))/max(1.,np.max(abs(h))))
    require(herm<1e-10,'nonhermitian LCFO Hamiltonian')
    # Match native diagonalization AFTER retaining its pre-symmetrization check.
    h=0.5*(h+h.conj().T)
    ortho=float(np.max(abs(coeff.conj().T@coeff-np.eye(coeff.shape[1]))))
    require(ortho<1e-10,'nonorthonormal LCFO eigenvectors')
    hc=h@coeff; eigenvalues=np.sum(coeff.conj()*hc,axis=0)
    residual=float(np.max(np.linalg.norm(hc-coeff*eigenvalues[None,:],axis=0)))
    require(residual<1e-8 and np.max(abs(eigenvalues.imag))<1e-10,'LCFO eigenpair residual failed')
    require(np.all(np.diff(eigenvalues.real)>=-1e-10),'unordered LCFO eigenstates')
    return dict(hamiltonian=h,coefficients=coeff,bases=bases,offsets=offsets,
                eigenvalues=eigenvalues.real,run_id=first['run'].decode().strip(),
                fragment_count=nf,dimension=n,origins=[origin_by_fragment[f] for f in range(1,nf+1)],
                grid=first['meta'][:3],core_grid=first['meta'][6:9],dv=float(first['geom'][9]),basis_orthogonality_error=basis_error,
                coefficient_orthogonality_error=ortho,eigen_residual_Ha=residual,
                hermiticity_relative=herm,
                directed_antihermitian_relative=float(np.linalg.norm(directed-directed.conj().T)/max(1.,np.linalg.norm(directed))),
                self_consistent_rt_validated=False)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('fragments',type=Path);parser.add_argument('output',type=Path)
    args=parser.parse_args();data=load_lcfo(args.fragments)
    summary={key:value for key,value in data.items() if not isinstance(value,(np.ndarray,list))}
    summary['lowest_eigenvalues_Ha']=data['eigenvalues'][:8].tolist()
    args.output.write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
