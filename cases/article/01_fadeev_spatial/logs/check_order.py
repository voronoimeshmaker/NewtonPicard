import csv,math,statistics
from pathlib import Path
r=Path('/home/jflavio/Programas/NewtonPicard/cases/article/01_fadeev_spatial/Result/run_20261009T183114Z')
rows=list(csv.DictReader((r/'fadeev_convergence.csv').open()))
for fam in ['cartesian','hexagonal','random']:
 a=[q for q in rows if q['mesh_family']==fam]; xs=[math.log(float(q['h'])) for q in a]; ys=[math.log(float(q['E2_volume_weighted'])) for q in a]; xm=statistics.mean(xs); ym=statistics.mean(ys); print(fam,sum((x-xm)*(y-ym) for x,y in zip(xs,ys))/sum((x-xm)**2 for x in xs))