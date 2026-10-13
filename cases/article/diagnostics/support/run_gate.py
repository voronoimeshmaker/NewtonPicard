from pathlib import Path
import csv,datetime,hashlib,json,math,os,subprocess,sys
root=Path(__file__).resolve().parents[1]; repo=root.parents[2]; article=root.parent
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
out=root/'D0_gate/Result'/('run_'+stamp);out.mkdir(parents=True)
source=article/'03_picard'/(article/'03_picard/Result/LATEST.txt').read_text().strip()/'cases.csv'
rows=list(csv.DictReader(source.open())); selected=[]
for mesh in ('cartesian','hexagonal','random'):
    eligible=[r for r in rows if r['method']=='WLSQR' and r['mesh_family']==mesh and math.isfinite(float(r['rho_reference'])) and float(r['rho_reference'])>0]
    selected.extend(sorted(eligible,key=lambda r:int(r['case_id']))[:20])
assert len(selected)==60,[(r['mesh_family'],r['case_id']) for r in selected]
(out/'selected.tsv').write_text(''.join(' '.join(r[k] for k in ('mesh_family','case_id','alpha','theta_degrees','c','rho_reference'))+'\n' for r in selected))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def cmd(*args):
    r=subprocess.run(args,cwd=repo,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT);return r.stdout
protected=[article/'support',article/'consolidated_refinements',repo/'GSSLib']+list(article.glob('*/Result'))
def snapshot():
    return {str(p):sha(p) for d in protected if d.exists() for p in d.rglob('*') if p.is_file()}
before=snapshot();(out/'protected_before.json').write_text(json.dumps(before,indent=2))
exe=root/'build/diag_D0'
manifest={'utc_start':stamp,'commit':cmd('git','rev-parse','HEAD').strip(),'branch':cmd('git','branch','--show-current').strip(),'git_status':cmd('git','status','--porcelain'),'executable_sha256':sha(exe),'source':str(source),'source_sha256':sha(source),'selected_cases':[(r['mesh_family'],r['case_id']) for r in selected],'OMP_NUM_THREADS':os.environ.get('OMP_NUM_THREADS'),'cpu':cmd('lscpu'),'compiler':cmd('c++','--version'),'PETSc_SLEPc':cmd('pkg-config','--modversion','PETSc','SLEPc'),'ADOLC_CGAL':cmd('dpkg-query','-W','libadolc-dev','libcgal-dev'),'linked_libraries':cmd('ldd',str(exe)),'flags':(root/'build/CMakeFiles/diag_D0.dir/flags.make').read_text()}
libs=list(Path('/usr/local/lib/cmake/VoronoiMeshMaker').glob('*'))+list(Path('/usr/local/lib/cmake/VoroGradToolkit').glob('*'))
libs+=list(Path('/usr/local/lib').glob('*Vor*'))
manifest['installed_metadata_and_library_hashes']={str(p):sha(p) for p in libs if p.is_file()}
manifest['installed_version_metadata']={str(p):p.read_text() for p in libs if 'ConfigVersion' in p.name}
(out/'manifest.txt').write_text(json.dumps(manifest,indent=2))
with (out/'console.log').open('w') as log:
    proc=subprocess.Popen([str(exe),str(root/'D0_gate/config.yaml'),str(out/'selected.tsv'),str(out/'gate.csv')],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    for line in proc.stdout:print(line,end='',flush=True);log.write(line);log.flush()
    code=proc.wait()
after=snapshot(); differences=[k for k in before.keys()|after.keys() if before.get(k)!=after.get(k)]
manifest.update(utc_end=datetime.datetime.now(datetime.timezone.utc).isoformat(),exit_code=code,protected_files_changed=differences)
(out/'manifest.txt').write_text(json.dumps(manifest,indent=2))
(root/'D0_gate/Result/LATEST.txt').write_text('Result/'+out.name+'\n')
print('Result:',out,'exit:',code,'protected changes:',len(differences))
sys.exit(code or bool(differences))
