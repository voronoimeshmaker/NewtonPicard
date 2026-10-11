exec(open('/home/jflavio/Programas/NewtonPicard/cases/article/build/compare_runs.py').read().split('for name,idx')[0])
n=read(sorted((r/'03_picard/Result').glob('run_*'))[-1]/'cases.csv'); old=read(Path(hist[1]['source']).parent/'cases.csv'); om={(x['method'],x['mesh_family'],x['case_id']):x for x in old}
for fam in ['cartesian','hexagonal','random']:
 pairs=[(x,om[x['method'],x['mesh_family'],x['case_id']]) for x in n if x['mesh_family']==fam and x['converged']=='1']
 for field in ['final_reference_error','final_residual']:
  worst=max(pairs,key=lambda ab:abs(float(ab[0][field])-float(ab[1][field])))
  a,b=worst; print(fam,len(pairs),field,'max_abs_diff',abs(float(a[field])-float(b[field])),'new',a[field],'old',b[field],'id',a['case_id'])
print('snapshot rows',len(n))