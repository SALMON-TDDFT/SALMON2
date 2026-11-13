"""Export a JSON frame list to the public SALMON VR desktop file layout."""
import argparse
import json
from pathlib import Path
from salmonpy.salmon_vr import export_cubes


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('frames', help='JSON: {"dt_au":0.02,"frames":[{"step":0,"density":"drho.cube","ac":"Ac.cube"}]}')
    p.add_argument('directory', help='New viewer directory')
    p.add_argument('--zero-ac', action='store_true', help='Explicit display-only placeholder when all A files are unavailable')
    p.add_argument('--application', help='Public SALMON VR.app to copy beside config.json')
    p.add_argument('--launch',action='store_true',help='Launch the copied macOS app')
    args=p.parse_args()
    if args.launch and not args.application: p.error('--launch requires --application')
    path=Path(args.frames).resolve();data=json.loads(path.read_text(encoding='utf-8'));frames=data['frames']
    supplied=[f.get('ac') is not None for f in frames]
    if any(supplied) and not all(supplied): p.error('A must be supplied for every frame or none')
    ac=[path.parent/f['ac'] for f in frames] if all(supplied) and frames else None
    result=export_cubes([path.parent/f['density'] for f in frames],args.directory,
        steps=[f['step'] for f in frames],dt_au=data['dt_au'],ac_files=ac,
        missing_ac='zero' if args.zero_ac else 'error',application=args.application)
    print(result.directory)
    if args.launch: result.launch()


if __name__=='__main__': main()
