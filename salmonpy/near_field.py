"""Prescribed molecular fields from collocated spatial E/H samples.

Numerical fitting, solid harmonics, derivatives and residuals run in Fortran.
Prescribed/streamed grid input data and stateless electric response maps are provided here.
These maps do not by themselves run a live molecular/FDTD feedback driver.
"""
from dataclasses import dataclass
import ctypes as c
from pathlib import Path
import warnings
import numpy as np
from .errors import InputError, NumericalError, RuntimeFailure
from .field_input import _real
from .field_sampling import SampledField
from .sbe import _integer

__all__=['ElectricPotentialSeries','NearFieldSeries','NearFieldFrame','couple_path_samples','GridFieldSeries']


def _kernel(name):
    from ._material_handles import default_sbe_api
    api=default_sbe_api()
    try:
        probe=api.library.sp_solid_protocol;probe.argtypes=[];probe.restype=c.c_int
        if probe()!=1:raise RuntimeFailure('unsupported spherical-harmonic kernel protocol')
        if name.startswith('sp_near_'):
            near=api.library.sp_near_protocol;near.argtypes=[];near.restype=c.c_int
            if near()!=3:raise RuntimeFailure('radial-complete near-field kernel protocol 3 required')
        if name.startswith('sp_grid_'):
            grid=api.library.sp_grid_protocol;grid.argtypes=[];grid.restype=c.c_int
            if grid()!=1:raise RuntimeFailure('direct grid-field kernel protocol 1 required')
        fn=getattr(api.library,name)
    except AttributeError as exc:
        raise RuntimeFailure('install or explicitly build calculation library with spherical-harmonic support') from exc
    fn.restype=c.c_int
    signatures={'sp_solid_fit':[c.c_int,c.c_int]+[c.c_void_p]*3+[c.c_double]+[c.c_void_p]*2,
                'sp_solid_evaluate':[c.c_int,c.c_int]+[c.c_void_p]*2+[c.c_double]+[c.c_void_p]*3,
                'sp_molecular_select':[c.c_void_p]*3,
                'sp_molecular_interval_field':[c.c_void_p]*3,
                'sp_near_fit':[c.c_int,c.c_int]+[c.c_void_p]*3+[c.c_int,c.c_int]+[c.c_void_p]*3,
                'sp_near_evaluate':[c.c_int,c.c_int]+[c.c_void_p]*8,
                'sp_near_couple_samples':[c.c_int]+[c.c_void_p]*7,
                'sp_grid_fields':[c.c_int]+[c.c_void_p]*9,
                'sp_grid_potentials':[c.c_int]+[c.c_void_p]*9,
                'sp_grid_pack':[c.c_int]+[c.c_void_p]*6,
                'sp_grid_average':[c.c_int]+[c.c_void_p]*3,
                'sp_grid_electric_response':[c.c_int,c.c_int]+[c.c_void_p]*9+[c.c_double,c.c_void_p],
                'sp_grid_difference':[c.c_int]+[c.c_void_p]*3,
                'sp_grid_magnetization_current':[c.c_void_p]*4,
                'sp_grid_current_sum':[c.c_int]+[c.c_void_p]*3}
    fn.argtypes=signatures[name]
    return fn


def _ptr(value):return value.ctypes.data_as(c.c_void_p)


def _check(status):
    if status==1:raise InputError('invalid or rank-deficient electric-field expansion samples')
    if status==4:raise NumericalError('nonfinite spherical-harmonic arithmetic')
    if status:raise RuntimeFailure('spherical-harmonic allocation failed')


@dataclass(frozen=True,slots=True,eq=False)
class ElectricPotentialSeries:
    """Real orthonormal Ylm, l=1..lmax, ordered l then m=-l..l.

    V_e=sum c_lm (r/radius)^l Y_lm, E=grad V_e; no constant potential.
    Coefficients are in Hartree, positions/radius in bohr, time in atomic units.
    This regular expansion assumes no external sources inside the fitting region.
    """
    center: object
    radius: float
    lmax: int
    times: object
    coefficients: object
    relative_residual: object
    field_component: str = "total"

    def __post_init__(self):
        if self.field_component not in ('total','difference'):raise InputError('invalid field component')
        center=_real(self.center,(3,),'expansion center')
        radius=float(_real(self.radius,(),'expansion radius'))
        if radius<=0:raise InputError('positive expansion radius required')
        lmax=_integer(self.lmax,'lmax',1,8)
        times=np.asarray(self.times)
        if times.ndim!=1 or len(times)<2:raise InputError('at least two sample times required')
        times=_real(times,times.shape,'expansion times')
        if times[0]!=0. or np.any(np.diff(times)<=0):raise InputError('expansion times must start at zero and increase')
        if not np.allclose(np.diff(times),times[1],rtol=1.e-12,atol=0.):raise InputError('uniform sample timestep required')
        coefficients=_real(self.coefficients,(len(times),lmax*(lmax+2)),'expansion coefficients')
        residual=_real(self.relative_residual,(len(times),),'relative field residual')
        if np.any(residual<0):raise InputError('nonnegative field residual required')
        for key,value in (('center',center),('radius',radius),('lmax',lmax),('times',times),
                          ('coefficients',coefficients),('relative_residual',residual)):
            object.__setattr__(self,key,value)

    @classmethod
    def from_samples(cls,samples,*,center,radius,lmax=2,residual_warning=.05,reference=None):
        """Fit E sampled at common positions/times; pass field.sample_fields()['E'].

        Component-specific Yee locations must be aligned before constructing the
        SampledField. No interface extrapolation or missing resolution is inferred.
        An explicit reference subtracts aligned E samples in Fortran, for example
        a vacuum run with the same incident source. Large residuals warn and remain available to the caller; they do not stop fitting.
        """
        if not isinstance(samples,SampledField) or samples.quantity!='E' or samples.unit!='E:au' or samples.position_unit!='bohr':
            raise InputError('collocated E:au samples with bohr positions required')
        if reference is not None:
            if (not isinstance(reference,SampledField) or reference.quantity!='E' or
                    reference.unit!=samples.unit or reference.position_unit!=samples.position_unit or
                    not np.array_equal(reference.positions,samples.positions) or
                    not np.array_equal(reference.times,samples.times)):
                raise InputError('reference must have the same E units, positions and times')
        lmax=_integer(lmax,'lmax',1,8);origin=_real(center,(3,),'expansion center')
        radius=float(_real(radius,(),'expansion radius'))
        if radius<=0:raise InputError('positive expansion radius required')
        residual_warning=float(_real(residual_warning,(),'residual warning threshold'))
        if residual_warning<0:raise InputError('nonnegative residual warning threshold required')
        points=np.ascontiguousarray(samples.positions);coeff=np.empty((len(samples.times),lmax*(lmax+2)));residual=np.empty(len(samples.times))
        fn=_kernel('sp_solid_fit')
        for i,frame in enumerate(samples.values):
            field=np.ascontiguousarray(frame);stats=np.zeros(2)
            if reference is not None:
                selected=np.empty_like(field);subtract=_kernel('sp_molecular_select')
                for j in range(len(points)):
                    _check(subtract(_ptr(field[j]),_ptr(reference.values[i,j]),_ptr(selected[j])))
                field=selected
            _check(fn(len(points),lmax,_ptr(points),_ptr(field),_ptr(origin),radius,_ptr(coeff[i]),_ptr(stats)))
            residual[i]=stats[0]
        result=cls(origin,radius,lmax,samples.times,coeff,residual,"difference" if reference is not None else "total")
        if np.any(residual>residual_warning):
            warnings.warn('electric-field expansion residual exceeds threshold; inspect spatial resolution, expansion order and non-potential field components',UserWarning,stacklevel=2)
        return result

    def _evaluate(self,positions,step_index):
        step_index=_integer(step_index,'step_index',0,len(self.times)-1)
        points=np.asarray(positions)
        if points.ndim!=2 or points.shape[1]!=3 or not len(points):raise InputError('positions[N,3] required')
        points=np.ascontiguousarray(_real(points,points.shape,'positions'))
        potential=np.empty(len(points));electric=np.empty((len(points),3))
        _check(_kernel('sp_solid_evaluate')(len(points),self.lmax,_ptr(points),_ptr(self.center),self.radius,
                                           _ptr(self.coefficients[step_index]),_ptr(potential),_ptr(electric)))
        potential.flags.writeable=False;electric.flags.writeable=False
        return potential,electric

    def potential(self,positions,*,step_index):return self._evaluate(positions,step_index)[0]
    def electric_field(self,positions,*,step_index):return self._evaluate(positions,step_index)[1]

    def write(self,path):
        """Write an explicit immutable TDDFT input; never overwrite or launch a solver.

        The saved history must cover indices 0..max_steps+1. No interpolation is
        performed. Keep the file in the GS/RT input bundle as salmonpy-ylm.dat.
        """
        path=Path(path)
        with path.open('x',encoding='ascii') as stream:
            stream.write('SALMONPY_SOLID_FIELD_V1\n')
            stream.write(f'{self.lmax} {len(self.times)}\n')
            stream.write(' '.join(format(float(v),'.17e') for v in (*self.center,self.radius,self.times[1]))+'\n')
            for row in self.coefficients:stream.write(' '.join(format(float(v),'.17e') for v in row)+'\n')
            stream.write('# field_component='+self.field_component+'\n')
        return path


@dataclass(frozen=True, slots=True, eq=False)
class NearFieldSeries:
    """Spatial E/H in molecular vacuum, including radial solid-harmonic terms.

    Arrays are [time,(L+1)*(L+2)*(L+3)//6,xyz]. Order: total degree d=0..L,
    angular l=d,d-2,..., then m=-l..l; basis (r/radius)^d Y_lm. E and H use Gaussian atomic units;
    H=B in the vacuum region occupied by the molecule. No spin coupling implied.
    Optional decay_vector (inverse bohr) multiplies the basis by exp(-decay dot (r-center)).
    Center, radius, order and decay are fixed at fitting and retained in direct checkpoints.
    """
    center: object
    radius: float
    lmax: int
    times: object
    E_coefficients: object
    H_coefficients: object
    relative_residual: object = None
    decay_vector: object = (0., 0., 0.)

    def __post_init__(self):
        center = _real(self.center, (3,), 'near-field expansion center')
        decay = _real(self.decay_vector, (3,), 'decay vector in inverse bohr')
        radius = float(_real(self.radius, (), 'near-field radius'))
        if radius <= 0:
            raise InputError('positive near-field radius required')
        lmax = _integer(self.lmax, 'lmax', 0, 8)
        times = np.asarray(self.times)
        if times.ndim != 1 or len(times) < 2:
            raise InputError('at least two near-field sample times required')
        times = _real(times, times.shape, 'near-field times')
        if times[0] != 0. or np.any(np.diff(times) <= 0):
            raise InputError('near-field times must start at zero and increase')
        if not np.allclose(np.diff(times), times[1], rtol=1.e-12, atol=0.):
            raise InputError('uniform near-field timestep required')
        shape = (len(times), ((lmax+1)*(lmax+2)*(lmax+3)//6), 3)
        ce = _real(self.E_coefficients, shape, 'electric Ylm coefficients')
        ch = _real(self.H_coefficients, shape, 'magnetic Ylm coefficients')
        residual = self.relative_residual
        if residual is not None:
            residual = _real(residual, (len(times), 2), 'E/H relative residual')
            if np.any(residual < 0):
                raise InputError('nonnegative E/H residuals required')
        for key, value in (('center', center), ('radius', radius), ('lmax', lmax),
                           ('times', times), ('E_coefficients', ce), ('H_coefficients', ch),
                           ('relative_residual', residual), ('decay_vector', decay)):
            object.__setattr__(self, key, value)

    @classmethod
    def from_samples(cls, E, H, *, center=None, radius=None, lmax=1, residual_warning=.05, decay_vector=(0.,0.,0.)):
        """Fit collocated E/H; explicit center is recommended for a molecule.

        If omitted, center is the sample centroid and radius the farthest sample
        distance, both calculated in Fortran during the first fit. The same
        geometry is used for every subsequent frame. A fixed decay_vector in inverse bohr
        enables an exponential envelope; it is supplied, never inferred automatically. No field is discarded.
        """
        if (not isinstance(E, SampledField) or not isinstance(H, SampledField)
                or E.quantity != 'E' or H.quantity != 'H'
                or E.unit != 'E:au' or H.unit != 'H:au'
                or E.position_unit != 'bohr' or H.position_unit != 'bohr'
                or not np.array_equal(E.positions, H.positions)
                or not np.array_equal(E.times, H.times)):
            raise InputError('collocated E:au/H:au with matching bohr positions and times required')
        lmax = _integer(lmax, 'lmax', 0, 8)
        geometry = np.zeros(7)
        geometry[4:] = _real(decay_vector, (3,), 'decay vector in inverse bohr')
        if center is not None:
            geometry[:3] = _real(center, (3,), 'near-field center')
        if radius is not None:
            geometry[3] = float(_real(radius, (), 'near-field radius'))
            if geometry[3] <= 0:
                raise InputError('positive near-field radius required')
        warning = float(_real(residual_warning, (), 'residual warning threshold'))
        if warning < 0:
            raise InputError('nonnegative residual warning threshold required')
        points = np.ascontiguousarray(E.positions)
        ce = np.empty((len(E.times), ((lmax+1)*(lmax+2)*(lmax+3)//6), 3))
        ch = np.empty_like(ce)
        residual = np.empty((len(E.times), 2))
        fn = _kernel('sp_near_fit')
        for i in range(len(E.times)):
            fields = np.stack((E.values[i], H.values[i]))
            _check(fn(len(points), lmax, _ptr(points), _ptr(fields), _ptr(geometry),
                      int(i == 0 and center is None), int(i == 0 and radius is None),
                      _ptr(ce[i]), _ptr(ch[i]), _ptr(residual[i])))
        result = cls(geometry[:3], geometry[3], lmax, E.times, ce, ch, residual, geometry[4:])
        if np.any(residual > warning):
            warnings.warn('E/H Ylm fit residual exceeds threshold; inspect sampling and expansion order',
                          UserWarning, stacklevel=2)
        return result

    @classmethod
    def from_yee_samples(cls, E, H_left, H_right, **options):
        """Fortran-average H at n-1/2 and n+1/2 onto E at n, then fit.

        The caller supplies both brackets, including the initial condition.
        No missing endpoint is guessed or extrapolated.
        """
        if not isinstance(E,SampledField) or E.quantity!='E':
            raise InputError('electric samples required')
        for H in (H_left, H_right):
            if (not isinstance(H, SampledField) or H.quantity != 'H' or H.unit != 'H:au'
                    or H.position_unit != 'bohr' or not np.array_equal(E.positions, H.positions)
                    or len(H.times) != len(E.times)):
                raise InputError('matching magnetic half-step samples required')
        if (np.any(H_left.times >= E.times) or np.any(H_right.times <= E.times)
                or not np.allclose(.5*(H_left.times+H_right.times), E.times, rtol=1.e-12, atol=1.e-15)):
            raise InputError('H samples must bracket each E time symmetrically')
        values = np.empty_like(E.values)
        fn = _kernel('sp_molecular_interval_field')
        for i in range(len(E.times)):
            for j in range(len(E.positions)):
                _check(fn(_ptr(H_left.values[i,j]), _ptr(H_right.values[i,j]), _ptr(values[i,j])))
        H = SampledField('H', E.times, E.positions, values, 'H:au')
        return cls.from_samples(E, H, **options)

    def frame(self, step_index):
        step_index = _integer(step_index, 'step_index', 0, len(self.times)-1)
        return NearFieldFrame(self.center, self.radius, self.lmax, step_index,
                              self.times[step_index], self.E_coefficients[step_index], self.H_coefficients[step_index],
                              None if self.relative_residual is None else self.relative_residual[step_index], self.decay_vector)

    def _evaluate(self, positions, step_index):
        step_index = _integer(step_index, 'step_index', 0, len(self.times) - 1)
        points = np.asarray(positions)
        if points.ndim != 2 or points.shape[1] != 3 or len(points) == 0:
            raise InputError('positions[N,3] required')
        points = np.ascontiguousarray(_real(points, points.shape, 'near-field probe positions'))
        geometry = np.array([*self.center, self.radius, *self.decay_vector])
        E = np.empty((len(points), 3)); H = np.empty_like(E)
        scalar = np.empty(len(points)); vector = np.empty_like(E)
        _check(_kernel('sp_near_evaluate')(len(points), self.lmax, _ptr(points), _ptr(geometry),
                    _ptr(self.E_coefficients[step_index]), _ptr(self.H_coefficients[step_index]),
                    _ptr(E), _ptr(H), _ptr(scalar), _ptr(vector)))
        for value in (E, H, scalar, vector):
            value.flags.writeable = False
        return E, H, scalar, vector

    def fields(self, positions, *, step_index):
        return self._evaluate(positions, step_index)[:2]

    def potentials(self, positions, *, step_index):
        """Return electron scalar energy and A/c in the radial multipolar gauge."""
        return self._evaluate(positions, step_index)[2:]

    def write(self, path):
        """Prepare salmonpy-near-eh.dat without launching or overwriting a calculation."""
        path = Path(path)
        with path.open('x', encoding='ascii') as stream:
            stream.write('SALMONPY_NEAR_EH_V3\n')
            stream.write(f'{self.lmax} {len(self.times)}\n')
            stream.write(' '.join(format(float(v), '.17e') for v in
                                  (*self.center, self.radius, *self.decay_vector, self.times[1])) + '\n')
            stream.write('E:au H:au bohr atomic_time vacuum\n')
            for ce, ch in zip(self.E_coefficients, self.H_coefficients):
                stream.write(' '.join(format(float(v), '.17e') for v in
                                      (*ce.ravel(), *ch.ravel())) + '\n')
        return path


@dataclass(frozen=True, slots=True, eq=False)
class NearFieldFrame:
    """One spatial E/H coefficient frame, bound to center and field step number."""
    center: object
    radius: float
    lmax: int
    step_index: int
    time: float
    E_coefficients: object
    H_coefficients: object
    relative_residual: object = None
    decay_vector: object = (0., 0., 0.)

    def __post_init__(self):
        lmax = _integer(self.lmax, 'lmax', 0, 8)
        center = _real(self.center, (3,), 'near-field center')
        decay = _real(self.decay_vector, (3,), 'decay vector in inverse bohr')
        radius = float(_real(self.radius, (), 'near-field radius'))
        time = float(_real(self.time, (), 'near-field time'))
        step = _integer(self.step_index, 'step_index', 0, 2**31-2)
        if radius <= 0 or time < 0:
            raise InputError('positive radius and nonnegative time required')
        shape = (((lmax+1)*(lmax+2)*(lmax+3)//6), 3)
        ce = _real(self.E_coefficients, shape, 'electric Ylm coefficients')
        ch = _real(self.H_coefficients, shape, 'magnetic Ylm coefficients')
        residual=None if self.relative_residual is None else _real(self.relative_residual,(2,),'E/H residual')
        if residual is not None and np.any(residual<0):raise InputError('nonnegative E/H residual required')
        object.__setattr__(self,'relative_residual',residual)
        for key, value in (('center',center),('radius',radius),('lmax',lmax),('step_index',step),
                           ('time',time),('E_coefficients',ce),('H_coefficients',ch),('decay_vector',decay)):
            object.__setattr__(self,key,value)

    def write_live_input(self,path,*,max_steps,dt):
        """Prepare an initial held-field template for subsequent step_near_field calls.

        No future field is predicted. Supply both coefficient frames at each
        update. Calling ordinary step instead would retain the held template.
        """
        if self.step_index!=0 or self.time!=0.:
            raise InputError('live input starts from field step0')
        max_steps=_integer(max_steps,'max_steps',1,2**31-3)
        dt=float(_real(dt,(),'dt'))
        if dt<=0:raise InputError('positive dt required')
        path=Path(path)
        row=' '.join(format(float(v),'.17e') for v in (*self.E_coefficients.ravel(),*self.H_coefficients.ravel()))+'\n'
        with path.open('x',encoding='ascii') as stream:
            stream.write('SALMONPY_NEAR_EH_V3\n')
            stream.write(f'{self.lmax} {max_steps+2}\n')
            stream.write(' '.join(format(float(v),'.17e') for v in (*self.center,self.radius,*self.decay_vector,dt))+'\n')
            stream.write('E:au H:au bohr atomic_time vacuum\n')
            for _ in range(max_steps+2):stream.write(row)
        return path

    @classmethod
    def from_samples(cls,E,H,*,step_index,center,radius,lmax=1,residual_warning=.05,decay_vector=(0.,0.,0.)):
        """Fit one collocated E/H frame using the already chosen model center."""
        if (not isinstance(E,SampledField) or not isinstance(H,SampledField)
                or len(E.times)!=1 or not np.array_equal(E.times,H.times)
                or E.quantity!='E' or H.quantity!='H' or E.unit!='E:au' or H.unit!='H:au'
                or E.position_unit!='bohr' or H.position_unit!='bohr'
                or not np.array_equal(E.positions,H.positions)):
            raise InputError('one collocated E/H frame required')
        lmax=_integer(lmax,'lmax',0,8)
        geometry=np.array([*_real(center,(3,),'near-field center'),float(_real(radius,(),'radius')),
            *_real(decay_vector,(3,),'decay vector in inverse bohr')])
        ce=np.empty((((lmax+1)*(lmax+2)*(lmax+3)//6),3));ch=np.empty_like(ce);residual=np.empty(2)
        points=np.ascontiguousarray(E.positions);fields=np.stack((E.values[0],H.values[0]))
        _check(_kernel('sp_near_fit')(len(points),lmax,_ptr(points),_ptr(fields),_ptr(geometry),0,0,
                                    _ptr(ce),_ptr(ch),_ptr(residual)))
        threshold=float(_real(residual_warning,(),'residual warning threshold'))
        if threshold<0:raise InputError('nonnegative residual warning threshold required')
        if np.any(residual>threshold):warnings.warn('E/H Ylm fit residual exceeds threshold',UserWarning,stacklevel=2)
        return cls(geometry[:3],geometry[3],lmax,step_index,E.times[0],ce,ch,residual,geometry[4:])


def couple_path_samples(displacement, parameters, weights, E, H):
    """Integrate a supplied field along center+s*displacement, in Fortran.

    displacement: bohr[3]; parameters and quadrature weights: dimensionless[N].
    E/H: collocated samples[N,3] in E:au/H:au, molecular vacuum (B=H).
    Returns electron scalar energy (Hartree) and readonly A/c[3]. The provider
    chooses quadrature nodes/weights and checks convergence; no interpolation,
    field representation, time evolution or accuracy certification is implied.
    """
    r=_real(displacement,(3,),'path displacement in bohr')
    samples=np.asarray(parameters)
    if samples.ndim!=1 or len(samples)==0:raise InputError('nonempty path parameters[N] required')
    samples=_real(samples,samples.shape,'path parameters')
    weights=_real(weights,samples.shape,'path quadrature weights')
    electric=_real(E,(len(samples),3),'path E samples')
    magnetic=_real(H,(len(samples),3),'path H samples')
    scalar=np.empty(1);vector=np.empty(3)
    status=_kernel('sp_near_couple_samples')(len(samples),_ptr(r),_ptr(samples),_ptr(weights),
        _ptr(electric),_ptr(magnetic),_ptr(scalar),_ptr(vector))
    if status==1:raise InputError('path parameters must lie in [0,1] and quadrature weights must be nonnegative')
    if status==4:raise NumericalError('nonfinite molecular path-coupling arithmetic')
    if status:raise RuntimeFailure('molecular path coupling failed')
    vector.flags.writeable=False
    return float(scalar[0]),vector


def _grid_check(status):
    if status==1:raise InputError('invalid Cartesian grid or requested path outside the available E/H grid; supply a covering patch')
    if status==4:raise NumericalError('nonfinite Cartesian field interpolation or path integral')
    if status:raise RuntimeFailure('Cartesian field allocation or evaluation failed')


def _grid_center_h(times,left_times,right_times,left,right):
    times=np.asarray(times)
    if times.ndim!=1 or len(times)==0:raise InputError('grid-field times[T] required')
    times=_real(times,times.shape,'E times')
    lt=_real(left_times,times.shape,'left magnetic times');rt=_real(right_times,times.shape,'right magnetic times')
    if (np.any(lt>=times) or np.any(rt<=times) or not np.allclose(.5*(lt+rt),times,rtol=1e-12,atol=1e-15)):
        raise InputError('H samples must bracket E times symmetrically')
    left=np.asarray(left)
    if left.ndim<2 or left.shape[0]!=len(times) or left.shape[-1]!=3:
        raise InputError('magnetic samples[time,...,3] required')
    left=_real(left,left.shape,'left H samples');right=_real(right,left.shape,'right H samples')
    centered=np.empty_like(left)
    _grid_check(_kernel('sp_grid_average')(_integer(centered.size,'field value count',1,2**31-1),
        _ptr(left),_ptr(right),_ptr(centered)))
    return centered


@dataclass(frozen=True,slots=True,eq=False)
class GridFieldSeries:
    """Prepared direct-grid E/H history; no multipole or exponential fit.

    E/H[time,z,y,x,component] in atomic units; origin is bohr[3] for collocated
    fields or [E/H,component,axis] for staggered component origins. Spacing is
    bohr[3]. molecule_position is the FDTD position of the molecular integration
    center; center is that same point in TDDFT coordinates. No extrapolation.
    """
    origin: object
    spacing: object
    times: object
    E: object
    H: object
    molecule_position: object = None
    center: object = (0.,0.,0.)
    time_origin: float = 0.

    def __post_init__(self):
        times=np.asarray(self.times)
        if times.ndim!=1 or len(times)<2:raise InputError('at least two grid-field times required')
        times=_real(times,times.shape,'grid-field times')
        if times[0]!=0. or np.any(np.diff(times)<=0.) or not np.allclose(np.diff(times),times[1],rtol=1e-12,atol=0.):
            raise InputError('uniform grid-field times starting at zero required')
        E=np.asarray(self.E)
        if E.ndim!=5 or E.shape[0]!=len(times) or E.shape[-1]!=3 or min(E.shape[1:4])<2:
            raise InputError('E/H[time,z,y,x,3] with at least two nodes per axis required')
        E=_real(E,E.shape,'grid E:au');H=_real(self.H,E.shape,'grid H:au')
        origin=np.asarray(self.origin)
        if origin.shape==(3,):origin=np.broadcast_to(origin,(2,3,3))
        origin=_real(origin,(2,3,3),'E/H component origins in bohr')
        spacing=_real(self.spacing,(3,),'grid spacing in bohr')
        if np.any(spacing<=0):raise InputError('positive Cartesian grid spacing required')
        center=_real(self.center,(3,),'molecular integration center')
        position=_real(center if self.molecule_position is None else self.molecule_position,(3,),'molecular position in FDTD coordinates')
        time_origin=float(_real(self.time_origin,(),'source time origin in atomic units'))
        for key,value in (('origin',origin),('spacing',spacing),('times',times),('E',E),('H',H),
                          ('center',center),('molecule_position',position),('time_origin',time_origin)):
            object.__setattr__(self,key,value)

    @property
    def counts(self):return (self.E.shape[3],self.E.shape[2],self.E.shape[1])

    def window(self,*,step_index,dt):
        """Frames n,n+1 for molecular step n-1 -> n; never evolves a solver.

        Source time_origin remains provenance of the complete history. The
        returned pair uses local field-step labels required by step_grid_field.
        No interpolation or resampling of time is implicit.
        """
        n=_integer(step_index,'step_index',1,len(self.times)-2)
        dt=float(_real(dt,(),'dt'))
        if dt<=0 or not np.isclose(self.times[1],dt,rtol=1e-12,atol=0.):
            raise InputError('grid history must contain frames 0..N+1 at the molecular dt')
        return GridFieldSeries(self.origin,self.spacing,[0.,dt],self.E[n:n+2],self.H[n:n+2],
            molecule_position=self.molecule_position,center=self.center,time_origin=n*dt)

    def _evaluate(self,positions,step_index,potentials):
        step=_integer(step_index,'step_index',0,len(self.times)-1)
        points=np.asarray(positions)
        if points.ndim!=2 or points.shape[1]!=3 or len(points)==0:raise InputError('positions[N,3] in TDDFT coordinates required')
        points=_real(points,points.shape,'query positions')
        counts=np.array(self.counts,dtype=np.intc);data=np.stack((self.E[step],self.H[step]))
        first=np.empty(len(points) if potentials else (len(points),3));second=np.empty((len(points),3))
        fn=_kernel('sp_grid_potentials' if potentials else 'sp_grid_fields')
        _grid_check(fn(len(points),_ptr(counts),_ptr(self.origin),_ptr(self.spacing),
            _ptr(self.molecule_position),_ptr(self.center),_ptr(points),_ptr(data),_ptr(first),_ptr(second)))
        first.flags.writeable=False;second.flags.writeable=False
        return first,second

    def fields(self,positions,*,step_index):
        return self._evaluate(positions,step_index,False)

    def potentials(self,positions,*,step_index):
        return self._evaluate(positions,step_index,True)

    def difference(self,reference):
        """Subtract a matching reference in Fortran; no time or space resampling."""
        if not isinstance(reference,GridFieldSeries):raise InputError('GridFieldSeries reference required')
        if (self.counts!=reference.counts or self.time_origin!=reference.time_origin or
            any(not np.array_equal(getattr(self,k),getattr(reference,k)) for k in
                ('origin','spacing','times','molecule_position','center'))):
            raise InputError('grid reference geometry/time labels differ')
        result=[]
        for name in ('E','H'):
            total=getattr(self,name);own=getattr(reference,name);out=np.empty_like(total)
            _grid_check(_kernel('sp_grid_difference')(_integer(total.size,'grid sample count',1,2**31-1),_ptr(total),_ptr(own),_ptr(out)))
            result.append(out)
        return GridFieldSeries(self.origin,self.spacing,self.times,*result,
            molecule_position=self.molecule_position,center=self.center,time_origin=self.time_origin)

    def _charge_response(self,positions,old,new,volume_elements,dt,mode):
        points=np.asarray(positions)
        if points.ndim!=2 or points.shape[1]!=3 or not len(points):
            raise InputError('charge positions[N,3] required in TDDFT coordinates')
        points=_real(points,points.shape,'charge positions')
        old=_real(old,(len(points),),'old signed charge density')
        new=_real(new,(len(points),),'signed charge density')
        volumes=_real(volume_elements,(len(points),),'volume elements in bohr^3')
        if np.any(volumes<0):raise InputError('nonnegative volume elements required')
        dt=float(_real(dt,(),'response dt'))
        if dt<=0:raise InputError('positive response dt required')
        counts=np.array(self.counts,dtype=np.intc)
        result=np.empty(self.E.shape[1:])
        fn=_kernel('sp_grid_electric_response')
        _grid_check(fn(mode,len(points),_ptr(counts),_ptr(self.origin[0]),_ptr(self.spacing),
            _ptr(self.molecule_position),_ptr(self.center),_ptr(points),_ptr(volumes),_ptr(old),_ptr(new),dt,_ptr(result)))
        return _real(result,result.shape,'electric-grid response')

    def charge_polarization(self,positions,charge_density,*,volume_elements):
        """Electric-path transpose P[z,y,x,component] in e/bohr^2.

        Input is signed physical charge density e/bohr^3, not positive electron
        number density. The result belongs on E-component grids, not macro sites.
        This stateless operation does not advance FDTD or a molecular worker.
        """
        return self._charge_response(positions,charge_density,charge_density,volume_elements,1.,0)

    def charge_current(self,positions,old_charge_density,new_charge_density,*,volume_elements,dt):
        """Interval current from signed charge-density change, in atomic units.

        Positions/volumes must be identical at both endpoints. Fixed ionic charges
        cancel in the difference. This is electric polarization current; orbital
        magnetization response and a live FDTD bridge are not included.
        """
        return self._charge_response(positions,old_charge_density,new_charge_density,volume_elements,dt,1)

    def write(self,path):
        """Write prepared histories only; never launch TDDFT or overwrite files."""
        path=Path(path)
        with path.open('x',encoding='ascii') as stream:
            stream.write('SALMONPY_NEAR_GRID_V1\n')
            stream.write(' '.join(str(v) for v in (*self.counts,len(self.times)))+'\n')
            stream.write(' '.join(format(float(v),'.17e') for v in (*self.center,*self.molecule_position,self.times[1]))+'\n')
            stream.write(' '.join(format(float(v),'.17e') for v in self.origin.ravel())+'\n')
            stream.write(' '.join(format(float(v),'.17e') for v in self.spacing)+'\n')
            stream.write('E:au H:au bohr atomic_time vacuum\n')
            for E,H in zip(self.E,self.H):
                stream.write(' '.join(format(float(v),'.17e') for v in (*E.ravel(),*H.ravel()))+'\n')
            stream.write('#source_time_origin_au='+format(self.time_origin,'.17e')+'\n')
        return path

    @classmethod
    def from_samples(cls,E,H,*,counts,origin,spacing,molecule_position=None,center=(0.,0.,0.),time_origin=0.):
        """Pack complete collocated Cartesian point samples in Fortran, in any point order."""
        if (not isinstance(E,SampledField) or not isinstance(H,SampledField) or E.quantity!='E' or H.quantity!='H'
                or E.unit!='E:au' or H.unit!='H:au' or E.position_unit!='bohr' or H.position_unit!='bohr'
                or not np.array_equal(E.positions,H.positions) or not np.array_equal(E.times,H.times)):
            raise InputError('collocated E/H samples with matching times/positions in atomic units required')
        if len(counts)!=3:raise InputError('counts[nx,ny,nz] required')
        counts=np.array([_integer(v,'grid count',2,2**31-1) for v in counts],dtype=np.intc)
        origin=_real(origin,(3,),'collocated grid origin');spacing=_real(spacing,(3,),'grid spacing')
        values=np.empty((len(E.times),2,int(counts[2]),int(counts[1]),int(counts[0]),3))
        fn=_kernel('sp_grid_pack')
        for i in range(len(E.times)):
            fields=np.stack((E.values[i],H.values[i]))
            _grid_check(fn(len(E.positions),_ptr(counts),_ptr(origin),_ptr(spacing),_ptr(E.positions),_ptr(fields),_ptr(values[i])))
        return cls(origin,spacing,E.times,values[:,0],values[:,1],molecule_position,center,time_origin)

    @classmethod
    def from_yee_samples(cls,E,H_left,H_right,**options):
        """Fortran-center explicit n-1/2,n+1/2 magnetic brackets, then pack the grid."""
        if not isinstance(E,SampledField) or E.quantity!='E':raise InputError('electric samples required')
        for H in (H_left,H_right):
            if (not isinstance(H,SampledField) or H.quantity!='H' or H.unit!='H:au'
                    or H.position_unit!='bohr' or not np.array_equal(H.positions,E.positions) or len(H.times)!=len(E.times)):
                raise InputError('matching magnetic half-step samples required')
        centered=_grid_center_h(E.times,H_left.times,H_right.times,H_left.values,H_right.values)
        return cls.from_samples(E,SampledField('H',E.times,E.positions,centered,'H:au'),**options)

    @classmethod
    def from_yee_arrays(cls,origin,spacing,times,E,H_left,H_right,*,H_left_times,H_right_times,**options):
        """Center staggered structured H arrays in time, retaining component-specific spatial origins."""
        if np.asarray(E).shape!=np.asarray(H_left).shape:raise InputError('E/H structured shapes must match')
        H=_grid_center_h(times,H_left_times,H_right_times,H_left,H_right)
        return cls(origin,spacing,times,E,H,**options)
