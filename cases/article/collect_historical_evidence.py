from pathlib import Path
import json,csv,hashlib
root=Path('/home/jflavio/Programas/GSMaker/Old/GSMakerLib/papers/PaperIII')
out=Path('/home/jflavio/Programas/NewtonPicard/cases/article')
evidence=[]
for name in ['NewtonOnly10000','PicardOnly10000','MixedRhoInitial10000','MixedRhoEveryIteration10000']:
 for file in sorted((root/name/'Result').glob('*/summary.csv')):
  rows=list(csv.DictReader(file.open()))
  selected=[r for r in rows if r['method']=='WLSQR' and r['scope']=='method_mesh']
  evidence.append(dict(source=str(file),sha256=hashlib.sha256(file.read_bytes()).hexdigest(),wlsqr=selected))
(out/'reference/historical_summaries.json').write_text(json.dumps(evidence,indent=2))
for item in evidence:
 print(Path(item['source']).parts[-4],[(r['mesh'],r['median_algorithmic_seconds'],r['converged_to_reference']) for r in item['wlsqr']])
p=out/'AUDIT.md'
s=p.read_text();s+='''

## Correspondência com resultados históricos existentes (somente leitura)
Foi localizada a saída NewtonOnly10000/Result/run_20260920T052100Z/summary.csv.
A mediana WLSQR/cartesiana é 0.074931204 s, que corresponde a 0.074931 s citada
no artigo; hexagonal 0.1226631045 s e random 0.1913537505 s também correspondem
às tabelas arredondadas. São resultados históricos, não corridas deste pacote.
`reference/historical_summaries.json` registra os caminhos, hashes e linhas WLSQR
das quatro campanhas encontradas. O total histórico 43164 sistemas de Newton
inclui WLSQR e MGG; a soma apenas WLSQR é 21593. Ao pós-processar, explicitar
sempre se o conjunto contém uma ou ambas as reconstruções.
''';p.write_text(s)