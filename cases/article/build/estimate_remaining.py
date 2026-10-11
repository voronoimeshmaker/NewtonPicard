from pathlib import Path
import csv,collections,json
r=Path('/home/jflavio/Programas/NewtonPicard/cases/article'); old=Path('/home/jflavio/Programas/GSMaker/Old/GSMakerLib/papers/PaperIII/CampanhasComplementares')
def read(p):
 return [q for q in csv.DictReader(p.open()) if None not in q and all(v is not None for v in q.values())]
groups=collections.defaultdict(list)
for p in old.rglob('summary.csv'):
 parts=p.parts
 labels=[s for s in parts if s in ['amp_1e-4','amp_1e-2','amp_1e-1','amp_1e0','L1','L3','L4']]
 if not labels: continue
 label=labels[-1]
 for q in read(p):
  if q.get('scope')!='method_mesh' or q.get('method')!='WLSQR':continue
  expected=200 if label.startswith('amp') else 100
  if int(q['cases'])!=expected:continue
  strategy=next((s for s in parts if s in ['newton','picard','mixed_initial','mixed_updated']),None)
  if strategy is None or (label.startswith('L') and strategy!='picard'):continue
  groups[label,strategy,q['mesh']].append((str(p),float(q['mean_case_total_seconds'])*expected))
totals=collections.defaultdict(float)
for k,v in sorted(groups.items()):
 p,t=sorted(v)[-1];totals[k[0]]+=t;print(k,'seconds',round(t),'copies',len(v))
print('TOTALS HOURS', {k:round(v/3600,3) for k,v in totals.items()})
h=json.loads((r/'reference/historical_summaries.json').read_text());orig=read(Path(h[3]['source']).parent/'cases.csv');new=read(sorted((r/'05_updated_rho/Result').glob('run_*'))[-1]/'cases.csv');key=lambda q:(q['method'],q['mesh_family'],q['case_id']);done={key(q) for q in new}; om={key(q):q for q in orig}
old_done=sum(float(om[key(q)]['seconds_case_total']) for q in new);new_done=sum(float(q['seconds_case_total']) for q in new);remaining=sum(float(q['seconds_case_total']) for q in orig if key(q) not in done)
print('CURRENT',len(new),'ratio_new_old',new_done/old_done,'remaining_historical_hours',remaining/3600,'estimated_remaining_hours',remaining/3600*new_done/old_done)