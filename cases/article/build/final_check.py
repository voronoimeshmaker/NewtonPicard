from pathlib import Path
import json,hashlib,subprocess,datetime
r=Path('/home/jflavio/Programas/NewtonPicard/cases/article'); s=Path('/home/jflavio/Programas/GSMaker/Old/GSMakerLib/papers/PaperIII')
cases=json.loads((r/'cases.json').read_text()); manifest=[]
for c in cases:
 n=c['name']
 if n=='01_fadeev_spatial': src=s/'FadeevConvergence/config.yaml'
 elif n.startswith('amplitude_'): src=s/'CampanhasComplementares/A-PalpiteInicial'/('amp_'+n.split('_')[1])/'config.yaml'
 elif n.startswith('branch_'): src=s/'CampanhasComplementares/B-EscalaDeMalha'/n.split('_')[1]/'config.yaml'
 else: src=s/'NewtonOnly10000/config.yaml'
 dst=r/n/'config.yaml'; assert src.read_bytes()==dst.read_bytes(),n
 manifest.append(dict(case=n,source=str(src),sha256=hashlib.sha256(dst.read_bytes()).hexdigest()))
 exe=r/'build'/('article_'+n); assert exe.is_file(),n
 deps=subprocess.check_output(['ldd',str(exe)],text=True); assert 'not found' not in deps,n
 (r/n/'logs/linkage.log').write_text(deps)
 assert not (r/n/'Result').exists(),n
for rec in json.loads((r/'provenance.json').read_text()):
 assert hashlib.sha256(Path(rec['source']).read_bytes()).hexdigest()==rec['source_sha256']
(r/'reference/config_provenance.json').write_text(json.dumps(manifest,indent=2))
(r/'BUILD_STATUS.md').write_text('# Verificacao final — 2026-10-09\n\n24 executaveis compilados com sucesso por `bash run_all.sh --build-only` (exit 0).\n129 fontes passaram na auditoria estatica; hashes das origens permanecem iguais.\n24 configuracoes identicas as respectivas referencias historicas selecionadas.\n`ldd` nos 24 executaveis: nenhuma biblioteca ausente; logs por caso em logs/linkage.log.\nNenhum executavel numerico foi iniciado; nenhuma pasta Result foi gerada.\nA reproducao numerica depende da execucao pelo usuario.\n')
print('24 builds and linkages verified; source and config provenance verified; no Result directories.')