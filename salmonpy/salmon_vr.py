"""File bridge to the public SALMON VR desktop application (2025 Mac release).

Independent of the native solver and optional graphics libraries. CUBE inputs
must use atomic units and represent the same fixed-ion grid. No fields are
inferred from filenames; the caller explicitly pairs files and native steps.
"""
from dataclasses import dataclass
from pathlib import Path
import json
import os
import plistlib
import shlex
import shutil
import subprocess
import sys
import numpy as np
from .outputs import read_cube, _digest

AU_TIME_FS = 0.024188843265857


def _same_grid(reference, grid):
    if reference.values.shape != grid.values.shape:
        raise ValueError('CUBE grid shape differs between frames/fields')
    for name in ('origin', 'axes', 'atoms'):
        a, b = getattr(reference, name), getattr(grid, name)
        if a.shape != b.shape or not np.allclose(a, b, rtol=0., atol=1e-10):
            raise ValueError('CUBE grid or fixed atoms differ: '+name)


def _write_cube(path, grid, values, label):
    with path.open('x', encoding='ascii') as stream:
        stream.write(label+'\nAll values in atomic units\n')
        stream.write(f'{len(grid.atoms):5d} '+' '.join(f'{x:.16e}' for x in grid.origin)+'\n')
        for count, axis in zip(values.shape, grid.axes):
            stream.write(f'{count:5d} '+' '.join(f'{x:.16e}' for x in axis)+'\n')
        for atom in grid.atoms:
            stream.write(f'{int(atom[0]):5d} '+' '.join(f'{x:.16e}' for x in atom[1:])+'\n')
        # CUBE ordering: z varies fastest. Keep lines short for desktop readers.
        for row in values.reshape(-1, values.shape[-1]):
            for start in range(0, len(row), 6):
                stream.write(' '.join(f'{x:.16e}' for x in row[start:start+6])+'\n')


@dataclass(frozen=True)
class VRExport:
    directory: Path
    application: Path | None = None

    def launch(self):
        """Explicitly launch the copied macOS app; never change OS security settings."""
        if sys.platform != 'darwin': raise RuntimeError('public desktop launch currently supports macOS only')
        if self.application is None or not self.application.is_dir():
            raise ValueError('export with application=path_to_public_app before launching')
        info = self.application/'Contents/Info.plist'
        try:
            with info.open('rb') as stream:
                metadata = plistlib.load(stream)
        except (OSError, ValueError, plistlib.InvalidFileException) as exc:
            raise ValueError('missing or invalid application Contents/Info.plist') from exc
        name = metadata.get('CFBundleExecutable') if isinstance(metadata, dict) else None
        if not isinstance(name, str) or not name or name in ('.', '..') or '/' in name or '\\' in name:
            raise ValueError('invalid application CFBundleExecutable')
        binary = self.application/'Contents/MacOS'/name
        if not binary.is_file():
            raise ValueError('application executable is missing: '+str(binary))
        if not os.access(binary, os.X_OK):
            raise PermissionError(
                'Application executable has no execute permission. Follow the public Mac manual; '
                'if you trust this application, run:\nchmod u+x '+shlex.quote(str(binary))+
                '\nThen retry launch(). No permissions or security settings were changed.')
        subprocess.run(['/usr/bin/open', '-n', str(self.application)], check=True)


def export_cubes(density_files, directory, *, steps, dt_au, ac_files=None,
                 missing_ac='error', component='x', steps_per_second=100,
                 application=None):
    """Export density-change CUBEs and optional Ac=A/c CUBEs without interpolation.

    Paths are paired in supplied order with strictly increasing native step IDs
    (0..999999). dt_au is the native time step, NOT the output-frame interval.
    Caller guarantees atomic units and density difference relative to initial GS.
    missing_ac='zero' explicitly generates display placeholders, never solver data.
    Grid consistency is checked one frame at a time. App copying/launch are opt-in.
    """
    density = tuple(Path(p).expanduser().resolve() for p in density_files)
    ids = tuple(steps)
    if not density or len(ids) != len(density): raise ValueError('one step per density CUBE is required')
    if any(isinstance(i, (bool,np.bool_)) or not isinstance(i, (int,np.integer)) or not 0<=i<=999999 for i in ids):
        raise ValueError('steps must be integers in 0..999999 for the public viewer')
    if any(b<=a for a,b in zip(ids,ids[1:])): raise ValueError('steps must be strictly increasing')
    dt = float(dt_au)
    if isinstance(dt_au, (bool,np.bool_)) or not np.isfinite(dt) or dt<=0: raise ValueError('positive finite dt_au required')
    if component not in ('x','y','z'): raise ValueError('component must be x, y or z')
    if isinstance(steps_per_second,bool) or not isinstance(steps_per_second,(int,np.integer)) or steps_per_second<=0:
        raise ValueError('steps_per_second must be a positive integer')
    if missing_ac not in ('error','zero'): raise ValueError('missing_ac must be error or zero')
    ac = None if ac_files is None else tuple(Path(p).expanduser().resolve() for p in ac_files)
    if ac is None and missing_ac!='zero': raise ValueError('A unavailable: explicitly select missing_ac="zero" for display placeholders')
    if ac is not None and len(ac)!=len(density): raise ValueError('one A CUBE per density frame is required')
    app = None if application is None else Path(application).expanduser().resolve()
    if app is not None and (app.suffix!='.app' or not (app/'Contents/MacOS').is_dir()):
        raise ValueError('application must be the public macOS .app bundle')
    dest = Path(directory).expanduser().resolve()
    if app is not None and dest.is_relative_to(app): raise ValueError('export must be outside the application bundle')
    dest.mkdir(parents=True,exist_ok=False)
    try:
        cubes=dest/'Ac_Je_drho_CUBE';cubes.mkdir()
        reference=None;frames=[]
        for index,(source,step) in enumerate(zip(density,ids)):
            grid=read_cube(source)
            if np.any(np.diag(grid.axes)<=0) or not np.allclose(grid.axes,np.diag(np.diag(grid.axes)),rtol=0.,atol=1e-12):
                raise ValueError('public viewer bridge requires a positive axis-aligned grid')
            if reference is None: reference=grid
            else: _same_grid(reference,grid)
            if np.any(grid.atoms[:,0]!=np.floor(grid.atoms[:,0])): raise ValueError('integer atomic numbers required')
            prefix=f'it{step:06d}'
            _write_cube(cubes/(prefix+'_drho.cube'),grid,grid.values,'Density change; SALMONpy export')
            if ac is None:
                field=np.zeros_like(grid.values);ac_hash=None
            else:
                agrid=read_cube(ac[index]);_same_grid(reference,agrid)
                field=agrid.values;ac_hash=_digest(ac[index])
            _write_cube(cubes/(prefix+'_Ac_'+component+'.cube'),grid,field,
                        'A unavailable; zero display placeholder' if ac is None else 'Ac=A/c; SALMONpy export')
            frames.append({'step':int(step),'time_au':float(step*dt),'density_sha256':_digest(source),'ac_sha256':ac_hash})
        note=f'dt={dt*AU_TIME_FS:.8g} fs/native step; Ac=A/c'
        if ac is None: note+='; A unavailable: zero placeholder'
        config={'CubeFilesDirectoryName':cubes.name,'InfomationText':note,
                'TimestepPerSecond':int(steps_per_second),'DrhoFileFormat':'*drho.cube',
                'AcFileFormat':'*Ac_'+component+'.cube','EnableAcGammaCorrection':False,'AcGamma':2.2}
        (dest/'config.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
        record={'schema':'salmonpy.salmon-vr/1','target':'public macOS 2025-02-28',
                'citation':'https://doi.org/10.1007/s10055-026-01319-0',
                'ac_source':'zero_placeholder' if ac is None else 'provided_cube',
                'component':component,'units':'atomic','dt_au':dt,'frames':frames,
                'viewer_interpolation':'linear between native step IDs; no export interpolation',
                'zero_placeholder_rendering':'public Mac smoke-tested; see validation report',
                'files':{p.relative_to(dest).as_posix():_digest(p) for p in cubes.iterdir()}}
        (dest/'salmonpy-vr.json').write_text(json.dumps(record,indent=2,allow_nan=False),encoding='utf-8')
        copied=None
        if app is not None:
            copied=dest/app.name;shutil.copytree(app,copied,symlinks=True)
        return VRExport(dest,copied)
    except BaseException:
        shutil.rmtree(dest)
        raise
