from pathlib import Path
import csv,datetime,hashlib,json,math,os,subprocess,sys,time,collections
from select_cases import select
root=Path(__file__).resolve().parents[1]; repo=root.parents[2]; article=root.parent
kind=sys.argv[1];folder=next(root.glob(kind+'_*/'))
gate=root/'D0_gate'/(root/'D0_gate/Result/LATEST.txt').read_text().strip()
g=list(csv.DictReader((gate/'gate.csv').open()));assert len(g)==60 and all(r['pass']=='1' for r in g),'D0 required'
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');out=folder/'Result'/('run_'+stamp);out.mkdir(parents=True)
selected,sources=select(article,kind,out)
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
with (out/'console.log').open('w') as log:
    proc=subprocess.Popen([str(exe),str(folder/'config.yaml'),str(out/'selected.tsv'),str(out)],cwd=root,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    for line in proc.stdout:print(line,end='',flush=True);log.write(line);log.flush()
    code=proc.wait()
after=snapshot(); differences=[k for k in before.keys()|after.keys() if before.get(k)!=after.get(k)]
manifest.update(utc_end=datetime.datetime.now(datetime.timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-start,exit_code=code,protected_files_changed=differences,selected_count=len(selected))
(out/'manifest.txt').write_text(json.dumps(manifest,indent=2));(folder/'Result/LATEST.txt').write_text('Result/'+out.name+'\n')
print('Result:',out,'exit:',code,'protected changes:',len(differences));sys.exit(code or bool(differences))
