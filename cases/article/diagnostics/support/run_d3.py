from pathlib import Path
import csv,datetime,hashlib,json,math,os,subprocess,sys,time,collections
from select_cases import select
root=Path(__file__).resolve().parents[1]; repo=root.parents[2]; article=root.parent
kind='D3';folder=root/'D3_hex_floor_location'
gate=root/'D0_gate'/(root/'D0_gate/Result/LATEST.txt').read_text().strip()
g=list(csv.DictReader((gate/'gate.csv').open()));assert len(g)==60 and all(r['pass']=='1' for r in g),'D0 required'
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');out=folder/'Result'/('run_'+stamp);out.mkdir(parents=True)
import random
sources={'L1':str(article/'branch_L1_picard/Result/run_20261010T182840Z/cases.csv'),'L3':str(article/'consolidated_refinements/branch_L3_picard/cases.csv')}
selected=[]
for level,source in sources.items():
    rows=[r for r in csv.DictReader(open(source)) if r['method']=='WLSQR']; chosen=[r for r in rows if float(r['regularisation_floor_active_iterations'])>0]
    controls=[r for r in rows if float(r['regularisation_floor_active_iterations'])==0]
    # ten controls per mesh and level provide comparisons even when all floor cases are hexagonal.
    rng=random.Random(20261013)
    for mesh in ('cartesian','hexagonal','random'):
        pool=[r for r in controls if r['mesh_family']==mesh]
        chosen+=rng.sample(pool,min(10,len(pool)))
    (out/('selected_'+level+'.tsv')).write_text(''.join(' '.join(r[k] for k in ('mesh_family','case_id','alpha','theta_degrees','c'))+'\n' for r in chosen))
    selected.extend(dict(r,level=level,source=source) for r in chosen)
(out/'selection.json').write_text(json.dumps({'seed':20261013,'cases':selected,'sources':sources},indent=2))
(out/'selected_empty.tsv').write_text('')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cmd(*args):
    r=subprocess.run(args,cwd=repo,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);return r.stdout
protected=[article/'support',article/'consolidated_refinements',repo/'GSSLib']+list(article.glob('*/Result'))
def snapshot():
    return {str(p):sha(p) for d in protected if d.exists() for p in d.rglob('*') if p.is_file()}
before=snapshot();(out/'protected_before.json').write_text(json.dumps(before,indent=2))

exe=root/'build'/('diag_'+kind)
manifest={'utc_start':stamp,'commit':cmd('git','rev-parse','HEAD').strip(),'branch':cmd('git','branch','--show-current').strip(),'git_status':cmd('git','status','--porcelain'),'executable_sha256':sha(exe),'sources':{s:sha(Path(s)) for s in sources.values()},'selection_sha256':sha(out/'selection.json'),'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),'cpu':cmd('lscpu'),'compiler':cmd('c++','--version'),'PETSc_SLEPc':cmd('pkg-config','--modversion','PETSc','SLEPc'),'ADOLC_CGAL':cmd('dpkg-query','-W','libadolc-dev','libcgal-dev'),'linked_libraries':cmd('ldd',str(exe)),'flags':(root/('build/CMakeFiles/diag_'+kind+'.dir/flags.make')).read_text(),'D0_manifest':str(gate/'manifest.txt'),'D0_executable_sha256':json.loads((gate/'manifest.txt').read_text())['executable_sha256']}
manifest['installed_libraries']=json.loads((gate/'manifest.txt').read_text())['installed_metadata_and_library_hashes']
(out/'manifest.txt').write_text(json.dumps(manifest,indent=2));start=time.monotonic()
code=0
with (out/'console.log').open('w') as log:
    for level,nx,ny in [('L1',128,43),('main',256,85),('L3',512,171),('L4',1024,341)]:
        selected_file=out/('selected_'+level+'.tsv') if level in sources else out/'selected_empty.tsv'
        proc=subprocess.Popen([str(exe),str(folder/'config.yaml'),str(selected_file),str(out),level,str(nx),str(ny)],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
        for line in proc.stdout:print(line,end='',flush=True);log.write(line);log.flush()
        rc=proc.wait();code=code or rc
        if rc:break
if code==0:
    import numpy as np
    from scipy.spatial import cKDTree
    with (out/'d3_mesh_symmetry.csv').open('w') as stream:
        writer=csv.writer(stream);writer.writerow(['level','mesh_family','cells','origin_count','min_origin_distance','max_reflection_nearest_distance'])
        for level in ('L1','main','L3','L4'):
            groups=collections.defaultdict(list)
            for r in csv.DictReader((out/('coordinates_'+level+'.csv')).open()):groups[r['mesh_family']].append((float(r['x']),float(r['y'])))
            for mesh,pts in groups.items():
                xy=np.array(pts);distance=np.linalg.norm(xy,axis=1);reflected=cKDTree(xy).query(-xy,workers=1)[0]
                writer.writerow([level,mesh,len(xy),int((distance<1e-12).sum()),format(float(distance.min()),'.17g'),format(float(reflected.max()),'.17g')])
    for name,target in [('floor','d3_floor_cells.csv'),('counts','d3_case_counts.csv')]:
        with (out/target).open('w') as dest:
            for j,level in enumerate(('L1','main','L3','L4')):
                lines=(out/(name+'_'+level+'.csv')).read_text().splitlines();dest.write('\n'.join(lines if j==0 else lines[1:])+'\n')
after=snapshot(); differences=[k for k in before.keys()|after.keys() if before.get(k)!=after.get(k)]
manifest.update(utc_end=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,exit_code=code,protected_files_changed=differences,selected_count=len(selected))
(out/'manifest.txt').write_text(json.dumps(manifest,indent=2));(folder/'Result/LATEST.txt').write_text('Result/'+out.name+'\n')
print('Result:',out,'exit:',code,'protected changes:',len(differences));sys.exit(code or bool(differences))
