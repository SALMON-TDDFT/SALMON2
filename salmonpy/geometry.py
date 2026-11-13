"""Material-grid previews for isolated, unrotated SALMON FDTD geometries."""
from dataclasses import dataclass
import numpy as np
from .ld import _array


def vector(value,name,positive=False):
    a=_array(value,name,(3,))
    if positive and np.any(a<=0): raise ValueError(f'{name} must be positive')
    return tuple(map(float,a))


@dataclass(frozen=True,slots=True)
class Grid:
    counts: tuple
    spacing: tuple

    def __post_init__(self):
        n=np.asarray(self.counts)
        if n.shape!=(3,) or n.dtype.kind not in 'iu' or np.any(n<2):
            raise ValueError('counts must be three integers >= 2')
        object.__setattr__(self,'counts',tuple(map(int,n)))
        object.__setattr__(self,'spacing',vector(self.spacing,'spacing',True))

    @property
    def coordinates(self):
        return tuple((np.arange(n)-(n-1)/2)*d for n,d in zip(self.counts,self.spacing))

    @property
    def size(self): return tuple(n*d for n,d in zip(self.counts,self.spacing))


@dataclass(frozen=True,slots=True)
class Shape:
    kind: str
    size: tuple
    center: tuple=(0.,0.,0.)
    material: int=1

    def __post_init__(self):
        if self.kind not in ('box','ellipsoid'): raise ValueError('kind must be box or ellipsoid')
        object.__setattr__(self,'size',vector(self.size,'size',True))
        object.__setattr__(self,'center',vector(self.center,'center'))
        if isinstance(self.material,bool) or not isinstance(self.material,(int,np.integer)) or self.material<1:
            raise ValueError('material must be a positive integer')
        object.__setattr__(self,'material',int(self.material))


@dataclass(frozen=True,slots=True)
class Geometry:
    grid: Grid
    shapes: tuple=()

    def __post_init__(self):
        if not isinstance(self.grid,Grid): raise TypeError('grid must be Grid')
        shapes=tuple(self.shapes)
        if len(shapes)>1000 or any(not isinstance(s,Shape) for s in shapes):
            raise ValueError('up to 1000 Shape objects are supported')
        object.__setattr__(self,'shapes',shapes)

    def material_ids(self):
        """Return (nx,ny,nz) IDs; later shapes overwrite earlier ones as in SALMON."""
        xyz=np.meshgrid(*self.grid.coordinates,indexing='ij',sparse=True)
        ids=np.zeros(self.grid.counts,dtype=np.int32)
        for shape in self.shapes:
            r=[(x-c)/(s/2) for x,c,s in zip(xyz,shape.center,shape.size)]
            if shape.kind=='ellipsoid': mask=r[0]**2+r[1]**2+r[2]**2<=1
            else: mask=(np.abs(r[0])<=1)&(np.abs(r[1])<=1)&(np.abs(r[2])<=1)
            ids[mask]=shape.material
        return ids

    def slice(self,axis,index):
        if axis not in (0,1,2): raise ValueError('axis must be 0,1,2')
        return np.take(self.material_ids(),index,axis=axis)

    def boundary_contacts(self,faces):
        """Return PML face names containing matter; current SALMON uses 8 cells."""
        ids=self.material_ids();found=[]
        if len(faces)!=6: raise ValueError('six boundary faces are required')
        for i,face in enumerate(faces):
            if face!='pml': continue
            axis=i//2;indices=range(9) if i%2==0 else range(ids.shape[axis]-9,ids.shape[axis])
            if ids.shape[axis]<19: raise ValueError('PML with clear interfaces needs at least 19 cells along its axis')
            if np.any(np.take(ids,list(indices),axis=axis)!=0): found.append(('x-','x+','y-','y+','z-','z+')[i])
        return tuple(found)
