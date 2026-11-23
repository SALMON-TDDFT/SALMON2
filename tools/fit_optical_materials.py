"""Explicit offline refitting of the bundled CC0 optical reference subset.

python tools/fit_optical_materials.py NEW_RESULTS [--export-catalog PACKAGE_CATALOG]
No download, SALMON job or GS collection. Full trials go to NEW_RESULTS;
export is explicit, preserving original raw data/source identity.
"""
import argparse
from dataclasses import asdict
import csv,json,warnings,hashlib
from pathlib import Path
from importlib.resources import files
import numpy as np
from salmonpy.optical_materials import load_material
from salmonpy.ld import fit_lorentz_drude,LorentzDrude


def curate(catalog):
    """Choose reusable defaults without hiding accurate-but-stiff alternatives."""
    summaries=[]
    for entry in catalog['materials']:
        for row in entry['fits']:
            report=row['report'];dt=report['recommended_dt']
            row['parameters_sha256']=None if row['parameters'] is None else hashlib.sha256(json.dumps(row['parameters'],sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()
            score=None if dt is None else float(np.multiply(dt,report['valid_frequency_range'][1]))
            row['curation']={'resolution_score':score,'minimum_resolution_score':.01,
                'pointwise_weighted_max':row.get('validation',{}).get('pointwise_weighted_max'),
                'maximum_pointwise_weighted_error':.10,
                'eligible_default':bool(row['parameters'] is not None and report['accepted'] and (score is None or score>=.01) and row.get('validation',{}).get('pointwise_weighted_max',float('inf'))<=.10),
                'scope':'heuristic time-domain resolution screening, not stability certification'}
        usable=[r for r in entry['fits'] if r['parameters'] is not None]
        preferred=[r for r in usable if r['curation']['eligible_default']]
        chosen=min(preferred,key=lambda r:r['poles']) if preferred else (min(usable,key=lambda r:r['report']['relative_rms_error']) if usable else None)
        entry['recommended_fit']=chosen['poles'] if preferred else None
        entry['best_candidate']=None if chosen is None else chosen['poles']
        entry['quality']='screened_candidate' if preferred else ('needs_refinement' if chosen else 'no_finite_fit')
        summaries.append(dict(id=entry['id'],symbol=entry['symbol'],dataset=entry['dataset'],
            energy_min_eV=entry['fit_band_eV_actual'][0],energy_max_eV=entry['fit_band_eV_actual'][1],
            recommended_poles=entry['recommended_fit'],best_candidate_poles=entry['best_candidate'],
            epsilon_infinity=None if chosen is None else chosen['parameters']['epsilon_infinity'],
            relative_rms_error=None if chosen is None else chosen['report']['relative_rms_error'],
            weighted_rms_error=None if chosen is None else chosen['report']['weighted_rms_error'],
            pointwise_weighted_max=None if chosen is None else chosen.get('validation',{}).get('pointwise_weighted_max'),
            optimizer_converged=False if chosen is None else chosen['report']['optimizer_converged'],
            accuracy_dt_hint_au=None if chosen is None else chosen['report']['recommended_dt'],
            quality=entry['quality'],doi=entry['doi'],source_sha256=entry['source_sha256']))
    catalog.setdefault('curation',{}).update(selection='smallest count meeting residual, time-resolution and measured-point screens; all candidates retained',
        resolution_screening='fastest active internal scale no more than about 10 times the calibrated omega_max; heuristic only',
        maximum_pointwise_weighted_error=.10,pointwise_scope='abs(delta epsilon)/max(1,abs(epsilon)) at supplied measurement points; not a physical accuracy guarantee')
    return summaries


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('output',type=Path)
    p.add_argument('--export-catalog',type=Path)
    p.add_argument('--materials',nargs='+',help='only refit these dataset IDs; retain other entries')
    p.add_argument('--poles',type=int,nargs='+',default=[2,4,6,8])
    p.add_argument('--optimizer',choices=['native','scipy'],default='scipy')
    p.add_argument('--starts',type=int,default=6)
    p.add_argument('--max-iterations',type=int,default=600)
    p.add_argument('--rtol',type=float,default=.05)
    p.add_argument('--minimum-damping',type=float,default=0.,help='explicit linewidth floor in atomic units; modeling assumption')
    p.add_argument('--initial-recorded',action='store_true',help='include the retained coefficients as the first start')
    p.add_argument('--recorded-settings',action='store_true',help='replay each stored trial setting, including Ag10/12/20')
    p.add_argument('--energy-range',type=float,nargs=2,default=[1.5,6.])
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=False)
    low,high=a.energy_range
    if not 0<=low<high:p.error('increasing nonnegative energy range required')
    source=files('salmonpy').joinpath('optical_data/catalog.json')
    catalog=json.loads(source.read_text(encoding='utf-8'))
    if a.materials and set(a.materials)-{e['id'] for e in catalog['materials']}:p.error('unknown material ID')
    summaries=[]
    for entry in catalog['materials']:
        if a.materials and entry['id'] not in a.materials:continue
        data=load_material(entry['id'])
        entry_low,entry_high=entry['fit_band_eV_requested'] if a.recorded_settings else (low,high)
        mask=(data.energy_eV>=entry_low)&(data.energy_eV<=entry_high)
        if np.count_nonzero(mask)<2:raise ValueError('not enough samples in fit band: '+entry['id'])
        plans=[(v['poles'],v['report']['starts'],v['report']['max_iterations'],v['report']['rtol'],v['report']['provenance'].get('initial_model')) for v in entry['fits']] if a.recorded_settings else [(v,a.starts,a.max_iterations,a.rtol,None) for v in a.poles]
        stored={v['poles']:None if v['parameters'] is None else LorentzDrude(**{k:x for k,x in v['parameters'].items() if k!='units'}) for v in entry['fits']}
        entry['fits']=[]
        entry['fit_band_eV_requested']=[entry_low,entry_high]
        entry['fit_band_eV_actual']=[float(data.energy_eV[mask][0]),float(data.energy_eV[mask][-1])]
        entry['source_energy_eV_range']=[float(data.energy_eV[0]),float(data.energy_eV[-1])]
        entry['source_samples']=len(data.omega);entry['fit_samples']=int(np.count_nonzero(mask))
        for poles,starts,max_iterations,rtol,recorded_initial in plans:
            floor=(next(v['report']['minimum_damping_au'] for v in data.fits if v['poles']==poles)
                   if a.recorded_settings else a.minimum_damping)
            initial=stored.get(poles) if a.initial_recorded else (None if recorded_initial is None else LorentzDrude(**{k:v for k,v in recorded_initial.items() if k!="units"}))
            with warnings.catch_warnings(record=True):
                warnings.simplefilter('always')
                cfg,report=fit_lorentz_drude(data.omega[mask],nk=data.nk[mask],poles=poles,
                    units='atomic_gaussian',optimizer=a.optimizer,initial=initial,drude=entry['drude'],starts=starts,max_iterations=max_iterations,
                    rtol=rtol,minimum_damping=floor,source_file=data.source_path,reference=entry['reference'],
                    metadata={'dataset_id':entry['id'],'source_commit':catalog['source_commit'],
                              'doi':entry['doi'],'conditions':entry['conditions'],'license':entry['license'],
                              'energy_eV_requested':[entry_low,entry_high]})
            record={'poles':poles,'parameters':None if cfg is None else asdict(cfg),'report':report.as_dict()}
            if cfg is not None:
                measured=data.epsilon[mask];predicted=LorentzDrude.from_config(cfg).permittivity(data.omega[mask])
                record['validation']={'pointwise_weighted_max':float(np.max(np.abs(predicted-measured)/np.maximum(1.,np.abs(measured)))),
                    'scope':'manual measured-point comparison; no experimental matching in required CI'}
            entry['fits'].append(record)
            (a.output/(entry['id']+f'-poles-{poles}.json')).write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
            print(entry['id'],poles,report.status,report.relative_rms_error,flush=True)
    catalog['curation']={'refitted_materials':a.materials,'optimizer':a.optimizer,'initial_recorded':a.initial_recorded,'recorded_settings':a.recorded_settings,'requested_band_eV':[low,high],'pole_counts':a.poles,'starts':a.starts,
        'max_iterations':a.max_iterations,'rtol':a.rtol,'minimum_damping_au':a.minimum_damping,
        'scope':'scalar nonmagnetic CC0 measurement tables; residual screening, not independent experimental or FDTD validation'}
    summaries=curate(catalog)
    text=json.dumps(catalog,ensure_ascii=False,indent=2)+'\n'
    (a.output/'catalog.json').write_text(text,encoding='utf-8')
    with (a.output/'material_list.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=list(summaries[0]),lineterminator="\n");writer.writeheader();writer.writerows(summaries)
    if a.export_catalog:a.export_catalog.write_text(text,encoding='utf-8')


if __name__=='__main__':main()
