exec(open('/home/jflavio/Programas/NewtonPicard/cases/article/build/compare_runs.py').read().split('for name,idx')[0])
run=sorted((r/'05_updated_rho/Result').glob('run_*'))[-1]; ns=read(run/'cases.csv'); old=read(Path(hist[3]['source']).parent/'cases.csv'); key=lambda x:(x['method'],x['mesh_family'],x['case_id']); om={key(x):x for x in old}
for group in sorted(set((x['method'],x['mesh_family']) for x in ns)):
 pp=[(x,om[key(x)]) for x in ns if key(x)[:2]==group]; print(group,'n',len(pp))
 for f in ['seconds_algorithmic','outer_iterations','rho_evaluations','strategy_linear_solve_calls']:
  print(f,'median new/old',statistics.median(float(a[f]) for a,b in pp),statistics.median(float(b[f]) for a,b in pp),'sum new/old',sum(float(a[f]) for a,b in pp),sum(float(b[f]) for a,b in pp))
 print('max_reference_error',max(float(a['final_reference_error']) for a,b in pp),'max_diff',max(abs(float(a['final_reference_error'])-float(b['final_reference_error'])) for a,b in pp))
 print('switches',sum(int(a['solver_switches']) for a,b in pp),'rho_diff',max(abs(float(a['rho_reference'])-float(b['rho_reference'])) for a,b in pp))
for name in ['02_newton','03_picard','04_initial_rho','05_updated_rho']:
 run=sorted((r/name/'Result').glob('run_*'))[-1]; rows=read(run/'cases.csv'); print('CHECK',name,'rows',len(rows),'duplicates',len(rows)-len({key(x) for x in rows}),'bad_converged_residual',sum(x['converged']=='1' and (not math.isfinite(float(x['final_residual'])) or float(x['final_residual'])>1e-12+1e-10*max(float(x['initial_residual']),1)) for x in rows))