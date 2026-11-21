"""Reference-frequency, passive subpixel material coefficients (native numerics)."""
from dataclasses import dataclass
import ctypes as c
import numpy as np
from .errors import InputError, NumericalError, RuntimeFailure
from .ld import LorentzDrude, LDModelConfig, LDTensorModelConfig, TensorLorentzDrude, _array, _positive


@dataclass(frozen=True, slots=True)
class SubpixelSmoothing:
    """Opt-in interface preparation; infer reference_omega from a pulse carrier.

    Frequencies use Gaussian atomic units. Impulse excitation requires an
    explicit reference. The approximation matches that frequency, not a band.
    """
    reference_omega: float | None = None

    def __post_init__(self):
        if self.reference_omega is not None:
            object.__setattr__(self, 'reference_omega', _positive(self.reference_omega, 'reference_omega'))

    def resolve(self, carrier=None):
        value = self.reference_omega if self.reference_omega is not None else carrier
        if value is None:
            raise InputError('subpixel smoothing needs reference_omega for an impulse without a carrier')
        return _positive(value, 'reference_omega')


def _api(name, arguments):
    from ._sbe_handles import default_api
    try:
        function = getattr(default_api().library, name)
    except AttributeError as exc:
        raise RuntimeFailure('native subpixel module unavailable; explicitly rebuild or install a matching wheel') from exc
    function.argtypes = arguments
    function.restype = c.c_int
    return function


def _is_ld_material(material):
    return isinstance(material,(LorentzDrude,LDModelConfig,LDTensorModelConfig))


def _material_data(material):
    if material is None:
        return np.eye(3), np.zeros((0, 3, 3)), np.zeros(0), np.zeros(0)
    if isinstance(material, (LDModelConfig, LDTensorModelConfig)):
        material = LorentzDrude.from_config(material)
    if isinstance(material, LorentzDrude):
        config = material.config
        if isinstance(config, LDTensorModelConfig):
            background = np.array(config.epsilon_infinity)
            strength = np.array(config.strength)
        else:
            background = np.eye(3)*config.epsilon_infinity
            strength = np.array([np.eye(3)*f for f in config.strength])
        strength = strength*np.float64(config.plasma_frequency)**2
        if not np.isfinite(strength).all():
            raise InputError('subpixel oscillator weights overflow; rescale physical input units')
        keep = np.any(strength != 0., axis=(1, 2))
        return background, strength[keep], np.array(config.resonance)[keep], np.array(config.damping)[keep]
    from .fdtd import Dielectric
    if isinstance(material, Dielectric) and material.mu == 1. and material.conductivity == 0.:
        return np.eye(3)*material.epsilon, np.zeros((0, 3, 3)), np.zeros(0), np.zeros(0)
    raise InputError('reference subpixel module requires nonmagnetic lossless background or passive LD materials')


def effective_material(left, right, *, fraction, normal, reference_omega, units):
    """Return a passive tensor LD model matching a local planar interface.

    fraction is the left phase's volume fraction, normal its outward direction.
    None denotes vacuum. The native kernel freezes complex interface field
    correctors and transforms every real background/pole by their symmetric
    energy congruence. Pole frequencies/damping are retained. No complex static
    epsilon or negative oscillator strength is introduced. This is a local
    reference-frequency approximation; inspect its spectrum away from reference.
    """
    if not (_is_ld_material(left) or _is_ld_material(right)):
        raise InputError('reference subpixel mode requires an explicitly specified LD phase; SBE/TDDFT are not transformed')
    if units != 'atomic_gaussian':
        raise InputError('declare atomic_gaussian angular frequencies explicitly')
    omega = _positive(reference_omega, 'reference_omega')
    f = _positive(fraction, 'fraction', zero=True)
    if f > 1:
        raise InputError('fraction must lie in [0,1]')
    n = _array(normal, 'normal', (3,))
    a, b = _material_data(left), _material_data(right)
    background = np.asfortranarray(np.stack([a[0], b[0]], axis=2))
    strength = np.asfortranarray(np.concatenate([a[1], b[1]], axis=0).transpose(1, 2, 0))
    resonance = np.ascontiguousarray(np.r_[a[2], b[2]])
    damping = np.ascontiguousarray(np.r_[a[3], b[3]])
    owners = np.ascontiguousarray([1]*len(a[2])+[2]*len(b[2]), dtype=np.int32)
    nq = len(resonance)
    effective_background = np.empty((3, 3), order='F')
    effective_strength = np.empty((3, 3, nq), order='F')
    target = np.empty((3, 3), dtype=np.complex128, order='F')
    matched = np.empty_like(target, order='F')
    fn = _api('sp_subpixel_mix', [c.c_int, c.c_double, c.c_double]+[c.c_void_p]*10)
    status = fn(nq, omega, f, *(x.ctypes.data for x in (n, background, strength, resonance, damping,
                      owners, effective_background, effective_strength, target, matched)))
    if status:
        raise NumericalError('reference interface is singular or invalid; choose another reference or retain staircase assignment')
    model = TensorLorentzDrude(plasma_frequency=1.,
        strength=effective_strength.transpose(2, 0, 1) if nq else np.zeros((1, 3, 3)),
        resonance=resonance if nq else [0.], damping=damping if nq else [0.],
        epsilon_infinity=effective_background)
    error = float(np.max(np.abs(target-matched)))
    def frozen(x):
        return np.frombuffer(x.tobytes(), dtype=x.dtype).reshape(x.shape)
    report = dict(reference_omega=omega, units=units, fraction=f,
                  target_epsilon=frozen(target), matched_epsilon=frozen(matched),
                  reference_error=error, pole_count=nq,
                  method='frozen complex interface correctors; symmetric energy-congruence passive LD lift',
                  scope='local planar interface at one frequency; no whole-band or grid-convergence guarantee')
    return model, report


def _plane_fraction(normal, distance, spacing):
    n = _array(normal, 'normal', (3,));h = _array(spacing, 'spacing', (3,))
    out = np.empty(1)
    fn = _api('sp_subpixel_plane_fraction', [c.c_void_p,c.c_double,c.c_void_p,c.c_void_p])
    if fn(n.ctypes.data,float(distance),h.ctypes.data,out.ctypes.data):
        raise InputError('invalid plane or voxel')
    return float(out[0])


def _shape_cut(shapes, point, spacing):
    kind = np.array([1 if s.kind == 'box' else 2 for s in shapes], dtype=np.int32)
    material = np.array([s.material for s in shapes], dtype=np.int32)
    sizes = np.asfortranarray(np.array([s.size for s in shapes], dtype=float).reshape(-1,3).T)
    centers = np.asfortranarray(np.array([s.center for s in shapes], dtype=float).reshape(-1,3).T)
    p = _array(point, 'point', (3,));h = _array(spacing, 'spacing', (3,))
    pair = np.empty(2, dtype=np.int32);f = np.empty(1);n = np.empty(3);flag = np.empty(1,dtype=np.int32)
    fn = _api('sp_subpixel_shape_cut', [c.c_int]+[c.c_void_p]*10)
    if fn(len(kind), *(x.ctypes.data for x in (kind,sizes,centers,material,p,h,pair,f,n,flag))):
        raise InputError('invalid subpixel shape input')
    return dict(materials=pair,fraction=float(f[0]),normal=n,flag=int(flag[0]))


def _profile_asset(materials, geometry, omega, rtol=1e-10, maxiter=200):
    data = [_material_data(m) for m in (None, *materials)]
    from .fdtd import _number, _vector
    rows = ['SALMONPY_SUBPIXEL_REFERENCE 2',
            _number(omega)+' '+str(len(materials))+' '+str(len(geometry.shapes))+' '+str(sum(len(x[2]) for x in data)),_number(rtol)+' '+str(maxiter)]
    for m,(background,strength,resonance,damping) in enumerate(data):
        rows.extend([f'{m} {len(resonance)} {int(_is_ld_material((None,*materials)[m]))}', _vector(background.flat)])
        for matrix,w,g in zip(strength,resonance,damping):
            rows.extend([_vector((w,g)),_vector(matrix.flat)])
    for shape in geometry.shapes:
        rows.extend([f'{1 if shape.kind=="box" else 2} {shape.material}',_vector((*shape.size,*shape.center))])
    return '\n'.join(rows)+'\n'
