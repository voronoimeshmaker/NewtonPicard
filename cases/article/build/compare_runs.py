from pathlib import Path
import csv,json,collections,math,statistics
r=Path('/home/jflavio/Programas/NewtonPicard/cases/article')
def read(p):
 with p.open() as f: return [x for x in csv.DictReader(f) if None not in x and all(v is not None for v in x.values())]
hist=json.loads((r/'reference/historical_summaries.json').read_text())
for name,idx in [('02_newton',0),('03_picard',1)]:
 run=sorted((r/name/'Result').glob('run_*'))[-1]; new=read(run/'cases.csv'); old=read(Path(hist[idx]['source']).parent/'cases.csv'); key=lambda x:(x['method'],x['mesh_family'],x['case_id']); om={key(x):x for x in old}; print('\nCASE',name,'rows',len(new),'summary?',(run/'summary.csv').exists())
 for group in sorted(set((x['method'],x['mesh_family']) for x in new)):
  ns=[x for x in new if key(x)[:2]==group]; pairs=[(x,om[key(x)]) for x in ns if key(x) in om]; print('GROUP',group,len(ns),'matched',len(pairs))
  for field in ['converged','converged_to_reference','outer_iterations','krylov_iterations_total','termination']:
   dif=[(a['case_id'],a[field],b[field]) for a,b in pairs if a[field]!=b[field]]; print(field,'different',len(dif),'examples',dif[:8])
  for field in ['alpha','theta_degrees','rho_target','c','rho_reference','initial_residual','final_reference_error']:
   ds=[abs(float(a[field])-float(b[field])) for a,b in pairs if math.isfinite(float(a[field])) and math.isfinite(float(b[field]))]; print(field,'max_abs_diff',max(ds,default=0))
  print('success new/old',sum(x['converged']=='1' for x,y in pairs),sum(y['converged']=='1' for x,y in pairs),'reference new/old',sum(x['converged_to_reference']=='1' for x,y in pairs),sum(y['converged_to_reference']=='1' for x,y in pairs))
  print('median seconds new/old',statistics.median(float(x['seconds_algorithmic']) for x,y in pairs),statistics.median(float(y['seconds_algorithmic']) for x,y in pairs))
  if name=='02_newton':
   print('max final residual/error',max(float(x['final_residual']) for x in ns),max(float(x['final_reference_error']) for x in ns))
   for a,b in pairs:
    if a['converged_to_reference']!='1': print('NONREFERENCE',json.dumps(a),'OLD',json.dumps(b))