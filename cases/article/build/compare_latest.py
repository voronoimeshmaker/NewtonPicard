exec(open('/home/jflavio/Programas/NewtonPicard/cases/article/build/compare_runs.py').read().split('for name,idx')[0])
for name,idx in [('04_initial_rho',2),('05_updated_rho',3)]:
 run=sorted((r/name/'Result').glob('run_*'))[-1]; ns=read(run/'cases.csv'); old=read(Path(hist[idx]['source']).parent/'cases.csv'); key=lambda x:(x['method'],x['mesh_family'],x['case_id']); om={key(x):x for x in old}; pairs=[(x,om[key(x)]) for x in ns]; print(name,len(ns),'complete',(run/'summary.csv').exists())
 for f in ['converged','converged_to_reference','termination','outer_iterations','solver_switches']:
  print(f,'different',sum(a[f]!=b[f] for a,b in pairs))
 for m in ['WLSQR','MGG']:
  pp=[(a,b) for a,b in pairs if a['method']==m]; print(m,'n',len(pp),'converged',sum(a['converged']=='1' for a,b in pp),'reference',sum(a['converged_to_reference']=='1' for a,b in pp))
 if (run/'summary.csv').exists():
  for q in read(run/'summary.csv'):
   if q['scope']=='method_mesh' and q['method']=='WLSQR': print(q['mesh'],'median',q['median_algorithmic_seconds'])
ns=read(sorted((r/'03_picard/Result').glob('run_*'))[-1]/'cases.csv')
for label,pred in [('rho<1',lambda q:math.isfinite(float(q['rho_reference'])) and float(q['rho_reference'])<1),('rho>=1',lambda q:math.isfinite(float(q['rho_reference'])) and float(q['rho_reference'])>=1)]:
 a=[q for q in ns if q['method']=='WLSQR' and pred(q)];print('Picard',label,len(a),'reference',sum(q['converged_to_reference']=='1' for q in a),'other',sum(q['converged']=='1' and q['converged_to_reference']=='0' for q in a),'failed',sum(q['converged']=='0' for q in a))