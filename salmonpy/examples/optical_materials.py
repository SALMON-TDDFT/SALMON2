"""Offline CC0 optical library; no fitting, download, SALMON run or GS."""
import argparse
import json
from salmonpy.optical_materials import list_materials, load_material

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--material',default=None)
    parser.add_argument('--poles',type=int)
    args=parser.parse_args()
    if args.material is None:
        for row in list_materials():
            print(row['id'],row['recommended_fit'],row['fit_band_eV_actual'])
        return
    data=load_material(args.material)
    model,report=data.model(poles=args.poles)
    print(json.dumps({'metadata':data.metadata,'fit_report':report},indent=2,ensure_ascii=False))
    print('epsilon tensor at first fitted frequency:',model.permittivity_tensor(report['valid_frequency_range'][0]))
    print('tensor representation:',model.tensor_config.representation,'independent components:',model.tensor_config.independent_components)

if __name__=='__main__':main()
