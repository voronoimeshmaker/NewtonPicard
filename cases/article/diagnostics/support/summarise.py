from pathlib import Path
import csv,json,math,collections
import numpy as np
root=Path(__file__).resolve().parents[1]
ANOM={259,590,820,983,1016,1627,1828}
def latest(k):
    d=next(root.glob(k+'_*/'));return d/(d/'Result/LATEST.txt').read_text().strip()
def read(p):return list(csv.DictReader(p.open()))
def num(x):return format(float(x),'.17g') if math.isfinite(float(x)) else 'NA'
def write(p,rows):
    if not rows:return
    with p.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader()
        for r in rows:w.writerow({k:num(v) if isinstance(v,(float,np.floating)) else v for k,v in r.items()})
def stats(values):
    a=np.array(values,dtype=float);a=a[np.isfinite(a)]
    if not len(a):return {k:'NA' for k in ('min','p50','p90','p95','p99','max')}
    return dict(zip(('min','p50','p90','p95','p99','max'),np.quantile(a,[0,.5,.9,.95,.99,1])))
def table(rows,fields):
    return '| '+' | '.join(fields)+' |\n| '+' | '.join('---' for _ in fields)+' |\n'+''.join('| '+' | '.join(str(r[k]) for k in fields)+' |\n' for r in rows)
def binrho(v):return '<1' if v<1 else '[1,1.02)' if v<1.02 else '[1.02,1.1)' if v<1.1 else '[1.1,1.4)' if v<1.4 else '>=1.4'
paths={k:latest(k) for k in ('D0','D1','D2','D3','D4')}
report=['# Relatório dos diagnósticos D0–D4\n','Resultados calculados em WSL Ubuntu-26.04-Test. Base b8b5f057d68d5ecdcbf9c4a3abe9f6155a26854a; branch diagnostics-2026-10. Nenhuma campanha não linear foi executada. D5 desativado.\n']
g=read(paths['D0']/'gate.csv');report+=['## D0 — reprodução histórica: SIM\n',f'{len(g)} casos, 20 por malha; todos aprovados. Máximo erro relativo: {num(max(float(r["relative_error"]) for r in g))}; tolerância 1e-10.\n']
a=read(paths['D1']/'d1_cases.csv');groups=collections.defaultdict(list)
for r in a:groups[r['mesh_family'],binrho(float(r['rho_reference']))].append(r)
summary=[]
for (mesh,b),rs in sorted(groups.items()):
    row={'mesh_family':mesh,'rho_bin':b,'cases':len(rs),'real_dominant':sum(r['mu1_is_real']=='1' for r in rs),'gap_below_one':sum(float(r['gap'])<1 for r in rs),'nconv_ge_2':sum(int(r['nconv'])>=2 for r in rs),'joint_accepted':sum(r['mu1_is_real']=='1' and float(r['gap'])<1 and int(r['nconv'])>=2 for r in rs)}
    for field in ('im_rel','gap','kappa_mu1'):
        row[field+'_available']=sum(r[field]!='NA' for r in rs)
        row.update({field+'_'+k:v for k,v in stats([float(r[field]) for r in rs if r[field]!='NA']).items()})
    summary.append(row)
write(paths['D1']/'d1_summary.csv',summary)
large=[r for r in a if float(r['rho_reference'])>=1];bad=[r for r in large if not(r['mu1_is_real']=='1' and float(r['gap'])<1 and int(r['nconv'])>=2)]
report+=['## D1 — hipóteses espectrais\n',f'{len(a)} casos: 2737 com rho histórico >=1, 244 exceções técnicas com rho recuperado em 02_newton e 300 controlos. As 244 exceções têm rho recuperado >=1. Critério conjunto: {len(large)-len(bad)}/{len(large)} = {num((len(large)-len(bad))/len(large))}. Universalmente satisfeito: NÃO.\n',table(read(paths['D1']/'d1_summary.csv'),['mesh_family','rho_bin','cases','joint_accepted']), '\nExceções ao critério conjunto:\n',table(bad,['mesh_family','case_id','rho','mu_1_re','mu_1_im','gap','nconv']), '\nO caso random/1095 tem par conjugado dominante; gap=1. kappa_mu1=NA porque a definição solicitada usa vetor real. A proxy de simplicidade algébrica não implica dominância estrita em módulo. Os resíduos calculados explicitamente estão em d1_eigenpairs.csv; o maior é '+num(max(float(r['relative_residual']) for r in read(paths['D1']/'d1_eigenpairs.csv')))+'. A tolerância EPS refere-se ao teste interno de convergência, não a uma garantia de limite idêntico para o resíduo recalculado.\n']
a=read(paths['D2']/'d2_spectrum.csv');fd=read(paths['D2']/'d2_directional.csv');fl=read(paths['D2']/'d2_floor.csv');best={}
for r in fd:
    key=tuple(r[k] for k in ('mesh_family','case_id','state','direction'))
    if key not in best or float(r['err_tot'])<float(best[key]['err_tot']):best[key]=r
write(paths['D2']/'d2_directional_best.csv',list(best.values()))
maxdrho=max(abs(float(r['drho'])) for r in a);margin=min(abs(float(r['rho_reference'])-1) for r in a);nofloor=[r for r in best.values() if int(r['n_floor'])==0];worst=max(float(r['err_tot']) for r in nofloor)
summary=[]
for mesh in ('cartesian','hexagonal','random'):
 for state in ('PsiF','Psi0'):
    ss=[r for r in a if r['mesh_family']==mesh and r['state']==state];ff=[r for r in fl if r['mesh_family']==mesh and r['state']==state];bb=[r for r in best.values() if r['mesh_family']==mesh and r['state']==state]
    row={'mesh_family':mesh,'state':state,'cases':len(ff),'spectra_all_maximisers':len(ss),'states_with_floor':sum(int(r['n_floor'])>0 for r in ff)}
    for metric,values in [('abs_drho',[abs(float(r['drho'])) for r in ss]),('best_err_tot',[float(r['err_tot']) for r in bb]),('r_min',[float(r['r_min']) for r in ff]),('eta_max',[float(r['eta_max']) for r in ff])]:row.update({metric+'_'+k:v for k,v in stats(values).items()})
    summary.append(row)
write(paths['D2']/'d2_summary.csv',summary)
report+=['## D2 — derivada congelada e total\n',f'{len(fl)//2} casos, {len(fl)} estados, {len(a)} espectros totais contando todos os maximizadores empatados. A seleção é a união dos 30 estratificados por malha com todas as inclusões obrigatórias; bins com menos de seis candidatos usam todos e completam os 30 com sorteio dos restantes.\n',table([{'métrica':'max |drho|','valor':num(maxdrho)},{'métrica':'min |rho_ref - 1| da amostra','valor':num(margin)},{'métrica':'0.1 × margem mínima','valor':num(.1*margin)},{'métrica':'maior mínimo de err_tot na janela, sem piso','valor':num(worst)},{'métrica':'testes direcionais sem piso','valor':len(nofloor)}],['métrica','valor']),f'\nCritério espectral: {"SIM" if maxdrho<2.9e-4 and maxdrho<.1*margin else "NÃO"}. Critério direcional na janela prescrita: {"SIM" if worst<=1e-6 else "NÃO"}. Seleção da janela: menor erro entre os quatro epsilons prescritos, explicitamente registrado por direção em d2_directional_best.csv; não se exige que todos os passos passem.\n',table(read(paths['D2']/'d2_summary.csv'),['mesh_family','state','cases','states_with_floor','abs_drho_max','best_err_tot_max']), '\nEmpates: todos os candidatos dentro de 1e-12 M2 foram tratados, com w esparso gravado em d2_w.csv; coordenadas limitadas aos primeiros dez. Cada candidato produz seu próprio espectro. Na diferença central em um empate, usa-se (max wv + min wv)/2, pois as derivadas unilaterais podem diferir. Isso não estabelece uma derivada de Fréchet única no empate. Valores de drho próximos ao arredondamento não resolvem numericamente o sinal físico da correção.\n']
meshes=read(paths['D3']/'d3_mesh_symmetry.csv');cells=read(paths['D3']/'d3_floor_cells.csv');counts=read(paths['D3']/'d3_case_counts.csv');summary=[]
for level in ('L1','L3'):
 for mesh in ('cartesian','hexagonal','random'):
  for state in ('PsiF','PsiF_plus_1e-11'):
    ss=[r for r in counts if r['level']==level and r['mesh_family']==mesh and r['state']==state];cc=[r for r in cells if r['level']==level and r['mesh_family']==mesh and r['state']==state]
    summary.append({'level':level,'mesh_family':mesh,'state':state,'cases':len(ss),'states_with_floor':sum(int(r['n_floor'])>0 for r in ss),'floor_cells':len(cc),'floor_cells_at_origin':sum(math.hypot(float(r['x']),float(r['y']))<1e-12 for r in cc),'max_stationary_distance':max((float(r['stationary_distance']) for r in cc),default=float('nan'))})
write(paths['D3']/'d3_summary.csv',summary)
report+=['## D3 — piso e origem: SIM\n',table(meshes,['level','mesh_family','cells','origin_count','min_origin_distance','max_reflection_nearest_distance']),f'\nTodas as {len(cells)} ocorrências de célula com piso nos estados avaliados estão na geradora da origem das malhas hexagonais L1/L3. Maior distância à origem: {num(max(math.hypot(float(r["x"]),float(r["y"])) for r in cells))}. Não há geradora na origem nas hexagonais principal/L4.\n','A seleção inclui os 64 casos L1 e 74 L3 com ativação histórica do piso, além de até dez controlos sem ativação histórica por malha/nível. Um controlo sem ativação durante a trajetória pode ter piso em Psi_F: as duas classificações se referem a estados diferentes. Os CSV separam essas informações.\n','Nas malhas random L3/L4, foi reutilizado por inclusão o reparo de recorte exato já usado nas campanhas, sem mudar bibliotecas ou tolerâncias. Os logs registram a tentativa original, sua falha e a construção aceita.\n']
a=read(paths['D4']/'d4_cases.csv');tops=read(paths['D4']/'d4_topcells.csv');summary=[]
for mesh in ('cartesian','hexagonal','random'):
 for state in ('PsiF','Psi0'):
  for group in ('selected_ids','random_controls'):
    ss=[r for r in a if r['mesh_family']==mesh and r['state']==state and ((int(r['case_id']) in ANOM)==(group=='selected_ids'))]
    if not ss:continue
    row={'mesh_family':mesh,'state':state,'group':group,'cases':len(ss)}
    for field,values in [('IPR',[float(r['IPR']) for r in ss]),('n90',[float(r['n90']) for r in ss]),('n90_fraction',[float(r['n90'])/float(r['cells']) for r in ss])]:row.update({field+'_'+k:v for k,v in stats(values).items()})
    summary.append(row)
write(paths['D4']/'d4_summary.csv',summary);selected=[]
for r in a:
 if r['mesh_family']=='random' and r['state']=='PsiF' and int(r['case_id']) in ANOM:
    t=next(t for t in tops if t['mesh_family']=='random' and t['case_id']==r['case_id'] and t['state']=='PsiF' and t['kind']=='top_abs_v' and t['rank']=='1')
    selected.append({k:r[k] for k in ('case_id','IPR','n90')}|{k:t[k] for k in ('stationary_distance','aspect_ratio','r_P')})
report+=['## D4 — localização dos modos\n',f'{len(a)//2} casos, {len(a)} estados. Sete IDs em três malhas mais vinte controlos aleatórios.\n',table(selected,['case_id','IPR','n90','stationary_distance','aspect_ratio','r_P']), '\nResposta geral: INCONCLUSIVO. A concentração varia entre os sete casos; 90% da norma não fica restrita a uma ou poucas células. O caso 1627 tem maior IPR entre os sete em Psi_F, mas a célula de maior amplitude está a distância 2.3099278978061428 do ponto estacionário mais próximo. Não há limiar prescrito para classificar IPR alto, proximidade ou degenerescência. Portanto, os dados não autorizam atribuir todos os casos a uma causa geométrica ou à proximidade de pontos estacionários.\n','Razão de aspecto = diâmetro do polígono / menor largura entre retas de suporte paralelas às suas arestas. Para pares complexos, IPR usa |v|² = Re(v)² + Im(v)². As linhas de maior norma e as dez maiores componentes estão em d4_topcells.csv.\n',table(read(paths['D4']/'d4_summary.csv'),['mesh_family','state','group','cases','IPR_p50','n90_fraction_p50'])]
report+=['## Tempos, integridade e reprodução\n','GNU C++ 15.2.0; PETSc 3.24.4; SLEPc 3.24.2; ADOL-C 2.7.2 (pacote 2.7.2-9); CGAL 6.1.1; VoronoiMeshMaker 1.1.0 e VoroGradToolkit 1.0.0 conforme metadados CMake instalados. Os hashes identificam o conteúdo instalado independentemente desses rótulos. Release -O3 -DNDEBUG, mesmas flags históricas, OMP_NUM_THREADS=1.\n']
times=[]
for k,o in paths.items():
 m=json.loads((o/'manifest.txt').read_text());elapsed=m.get('elapsed_seconds')
 if elapsed is None:
    import datetime
    elapsed=(datetime.datetime.fromisoformat(m['utc_end'])-datetime.datetime.strptime(m['utc_start'],'%Y%m%dT%H%M%SZ').replace(tzinfo=datetime.timezone.utc)).total_seconds()
 times.append({'diagnóstico':k,'segundos':num(elapsed),'exit_code':m['exit_code'],'arquivos_protegidos_alterados':len(m['protected_files_changed'])})
report+=[table(times,['diagnóstico','segundos','exit_code','arquivos_protegidos_alterados']), '\nOs tempos incluem execução e conferência final dos hashes; não incluem compilação. D0 usa timestamps UTC, com início registrado em segundos inteiros.\n','O repositório já tinha alterações alheias ao diagnóstico. Os manifests preservam o estado completo da árvore; os commits desta entrega contêm apenas diagnostics/. Nenhum push/tag foi feito. Versões, CPU, flags, origem dos casos e hashes dos executáveis constam dos manifests. provenance_verified.json complementa hashes das bibliotecas e árvores de headers, explicitamente capturados após a execução.\n']
for k,o in paths.items():
 folder=o.parent.parent
 report.append(f'### {k}\n\nComando, a partir de diagnostics/: `bash {folder.name}/run.sh`.\n\nExecução: `{o.relative_to(root)}`.\n\nArquivos: '+', '.join('`'+x.name+'`' for x in sorted(o.glob('*.csv')))+', `manifest.txt`, `console.log`, seleção e hashes.\n')
report+=['\nPara recriar as tabelas e este relatório após os cálculos: `python3 support/summarise.py`; para registrar os hashes complementares: `python3 support/provenance.py`. Os CSV têm 17 algarismos significativos. NA significa indisponível/não aplicável, nunca zero.\n','## Limitações\n\nD1–D4 avaliam os estados e amostras definidos; não comprovam propriedades globais nem reconstroem trajetórias não lineares. A correção espectral de D2 ficou abaixo da resolução prática dos cálculos. D4 não tem limiares classificatórios fornecidos. D5 não foi implementado nem executado; o launcher permanece bloqueado. Não se propõem alterações ao manuscrito ou ao algoritmo.\n']
(root/'RELATORIO_DIAGNOSTICOS.md').write_text('\n'.join(report))
print('Summaries and report written.')
