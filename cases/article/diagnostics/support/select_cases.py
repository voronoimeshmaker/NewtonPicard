from pathlib import Path
import csv,random,math,json
SEED=20261013
FAMILIES=('cartesian','hexagonal','random')
ANOMALIES={590,1828,1627,983,820,259,1016}
def select(article,kind,out):
    sources={}
    def read(campaign):
        f=article/campaign/(article/campaign/'Result/LATEST.txt').read_text().strip()/'cases.csv';sources[campaign]=str(f)
        return [r for r in csv.DictReader(f.open()) if r['method']=='WLSQR']
    pic=read('03_picard');new={(r['mesh_family'],r['case_id']):r for r in read('02_newton')}
    rng=random.Random(SEED);chosen={}
    def add(r,reason):
        k=(r['mesh_family'],r['case_id'])
        if k not in chosen:chosen[k]=dict(r,selection_reason=reason,source=sources['03_picard'])
        else:chosen[k]['selection_reason']+=';'+reason
        if not math.isfinite(float(chosen[k]['rho_reference'])):
            chosen[k]['rho_reference']=new[k]['rho_reference'];chosen[k]['rho_source']=sources['02_newton']
    if kind=='D1':
        for r in pic:
            if float(r['rho_reference'])>=1:add(r,'rho_ge_1')
            elif r['termination'].startswith('case_exception:'):add(r,'technical_exception')
        for mesh in FAMILIES:
            for r in rng.sample([r for r in pic if r['mesh_family']==mesh and float(r['rho_reference'])<1],100):add(r,'control')
    elif kind=='D2':
        for mesh in FAMILIES:
            eligible=[r for r in pic if r['mesh_family']==mesh and math.isfinite(float(r['rho_reference']))]
            for low,high in ((0,1),(1,1.02),(1.02,1.1),(1.1,1.4),(1.4,float('inf'))):
                pool=[r for r in eligible if low<=float(r['rho_reference'])<high]
                for r in rng.sample(pool,min(6,len(pool))):add(r,'stratified_base_30')
            missing=30-sum(k[0]==mesh for k in chosen)
            for r in rng.sample([r for r in eligible if (mesh,r['case_id']) not in chosen],missing):add(r,'stratified_fill_30')
            for r in rng.sample([r for r in eligible if abs(float(r['rho_reference'])-1)<.01],10):add(r,'near_threshold')
        for r in pic:
            if float(r['rho_reference'])>=1 and r['converged_to_reference']=='1':add(r,'reference_recovery')
            if r['mesh_family']=='random' and int(r['case_id']) in ANOMALIES:add(r,'random_anomaly')
    elif kind=='D4':
        for r in pic:
            if int(r['case_id']) in ANOMALIES:add(r,'anomaly_or_same_id_control')
        for r in rng.sample([r for r in pic if r['mesh_family']=='random' and int(r['case_id']) not in ANOMALIES and abs(float(r['rho_reference'])/float(r['rho_target'])-1)<.01],20):add(r,'random_control')
    else:raise ValueError(kind)
    rows=sorted(chosen.values(),key=lambda r:(r['mesh_family'],int(r['case_id'])))
    (out/'selected.tsv').write_text(''.join(' '.join(r[k] for k in ('mesh_family','case_id','alpha','theta_degrees','c','rho_reference'))+'\n' for r in rows))
    (out/'selection.json').write_text(json.dumps({'seed':SEED,'cases':rows,'sources':sources},indent=2))
    return rows,sources
