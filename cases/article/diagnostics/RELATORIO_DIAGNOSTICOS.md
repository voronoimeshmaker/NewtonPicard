# Relatório dos diagnósticos D0–D4

Resultados calculados em WSL Ubuntu-26.04-Test. Base b8b5f057d68d5ecdcbf9c4a3abe9f6155a26854a; branch diagnostics-2026-10. Nenhuma campanha não linear foi executada. D5 desativado.

## D0 — reprodução histórica: SIM

60 casos, 20 por malha; todos aprovados. Máximo erro relativo: 0; tolerância 1e-10.

## D1 — hipóteses espectrais

3281 casos: 2737 com rho histórico >=1, 244 exceções técnicas com rho recuperado em 02_newton e 300 controlos. As 244 exceções têm rho recuperado >=1. Critério conjunto: 2980/2981 = 0.99966454209996647. Universalmente satisfeito: NÃO.

| mesh_family | rho_bin | cases | joint_accepted |
| --- | --- | --- | --- |
| cartesian | <1 | 100 | 100 |
| cartesian | [1,1.02) | 207 | 207 |
| cartesian | [1.02,1.1) | 400 | 400 |
| cartesian | [1.1,1.4) | 400 | 400 |
| hexagonal | <1 | 100 | 100 |
| hexagonal | [1,1.02) | 209 | 209 |
| hexagonal | [1.02,1.1) | 400 | 400 |
| hexagonal | [1.1,1.4) | 402 | 402 |
| random | <1 | 100 | 100 |
| random | >=1.4 | 6 | 6 |
| random | [1,1.02) | 181 | 181 |
| random | [1.02,1.1) | 389 | 389 |
| random | [1.1,1.4) | 387 | 386 |


Exceções ao critério conjunto:

| mesh_family | case_id | rho | mu_1_re | mu_1_im | gap | nconv |
| --- | --- | --- | --- | --- | --- | --- |
| random | 1095 | 1.1692010764320628 | 1.167385270366704 | 0.065136684447029136 | 1 | 8 |


O caso random/1095 tem par conjugado dominante; gap=1. kappa_mu1=NA porque a definição solicitada usa vetor real. A proxy de simplicidade algébrica não implica dominância estrita em módulo. Os resíduos calculados explicitamente estão em d1_eigenpairs.csv; o maior é 2.9859492842493973e-11. A tolerância EPS refere-se ao teste interno de convergência, não a uma garantia de limite idêntico para o resíduo recalculado.

## D2 — derivada congelada e total

164 casos, 328 estados, 385 espectros totais contando todos os maximizadores empatados. A seleção é a união dos 30 estratificados por malha com todas as inclusões obrigatórias; bins com menos de seis candidatos usam todos e completam os 30 com sorteio dos restantes.

| métrica | valor |
| --- | --- |
| max |drho| | 2.4868995751603507e-14 |
| min |rho_ref - 1| da amostra | 9.2481364417640322e-06 |
| 0.1 × margem mínima | 9.2481364417640326e-07 |
| maior mínimo de err_tot na janela, sem piso | 7.6048978459565142e-07 |
| testes direcionais sem piso | 656 |


Critério espectral: SIM. Critério direcional na janela prescrita: SIM. Seleção da janela: menor erro entre os quatro epsilons prescritos, explicitamente registrado por direção em d2_directional_best.csv; não se exige que todos os passos passem.

| mesh_family | state | cases | states_with_floor | abs_drho_max | best_err_tot_max |
| --- | --- | --- | --- | --- | --- |
| cartesian | PsiF | 57 | 0 | 8.8817841970012523e-15 | 8.7790357707244366e-08 |
| cartesian | Psi0 | 57 | 0 | 5.5511151231257827e-15 | 5.1766039610820483e-08 |
| hexagonal | PsiF | 58 | 0 | 4.8849813083506888e-15 | 9.2770762650881307e-08 |
| hexagonal | Psi0 | 58 | 0 | 6.4392935428259079e-15 | 1.0258578680283135e-07 |
| random | PsiF | 49 | 0 | 2.4868995751603507e-14 | 7.6048978459565142e-07 |
| random | Psi0 | 49 | 0 | 7.1054273576010019e-15 | 5.0981820828230322e-07 |


Empates: todos os candidatos dentro de 1e-12 M2 foram tratados, com w esparso gravado em d2_w.csv; coordenadas limitadas aos primeiros dez. Cada candidato produz seu próprio espectro. Na diferença central em um empate, usa-se (max wv + min wv)/2, pois as derivadas unilaterais podem diferir. Isso não estabelece uma derivada de Fréchet única no empate. Valores de drho próximos ao arredondamento não resolvem numericamente o sinal físico da correção.

## D3 — piso e origem: SIM

| level | mesh_family | cells | origin_count | min_origin_distance | max_reflection_nearest_distance |
| --- | --- | --- | --- | --- | --- |
| L1 | cartesian | 5504 | 0 | 0.1171875 | 8.8817841970012523e-16 |
| L1 | hexagonal | 6247 | 1 | 2.2204460492503131e-16 | 8.8817841970012523e-16 |
| L1 | random | 5504 | 0 | 0.062418944030661887 | 0.36974002130824479 |
| main | cartesian | 21760 | 0 | 0.05859375 | 4.4408920985006262e-16 |
| main | hexagonal | 25039 | 0 | 0.050743676002994453 | 0.05859375 |
| main | random | 21760 | 0 | 0.038846035251113992 | 0.21607629137613649 |
| L3 | cartesian | 87552 | 0 | 0.029296875 | 8.8817841970012523e-16 |
| L3 | hexagonal | 100765 | 1 | 2.3592239273284576e-16 | 8.8817841970012523e-16 |
| L3 | random | 87552 | 0 | 0.017712652652418913 | 0.13879264341608535 |
| L4 | cartesian | 349184 | 0 | 0.0146484375 | 8.8817841970012523e-16 |
| L4 | hexagonal | 403259 | 0 | 0.012685919000748905 | 0.0146484375 |
| L4 | random | 349184 | 0 | 0.017343953624255753 | 0.059377598551911256 |


Todas as 316 ocorrências de célula com piso nos estados avaliados estão na geradora da origem das malhas hexagonais L1/L3. Maior distância à origem: 2.3592239273284576e-16. Não há geradora na origem nas hexagonais principal/L4.

A seleção inclui os 64 casos L1 e 74 L3 com ativação histórica do piso, além de até dez controlos sem ativação histórica por malha/nível. Um controlo sem ativação durante a trajetória pode ter piso em Psi_F: as duas classificações se referem a estados diferentes. Os CSV separam essas informações.

Nas malhas random L3/L4, foi reutilizado por inclusão o reparo de recorte exato já usado nas campanhas, sem mudar bibliotecas ou tolerâncias. Os logs registram a tentativa original, sua falha e a construção aceita.

## D4 — localização dos modos

41 casos, 82 estados. Sete IDs em três malhas mais vinte controlos aleatórios.

| case_id | IPR | n90 | stationary_distance | aspect_ratio | r_P |
| --- | --- | --- | --- | --- | --- |
| 259 | 0.0010302738120993235 | 2594 | 0.10701298570794029 | 1.4707391219175396 | 0.00046839219895107278 |
| 590 | 0.0010971468895912486 | 2222 | 0.28358138534624067 | 1.5251431811197911 | 0.00034576356933549592 |
| 820 | 0.00041631254143867401 | 4427 | 0.059598474921310048 | 1.4762057700768041 | 0.00016304765970147822 |
| 983 | 0.00027507294232954665 | 5601 | 0.35777956742058181 | 2.2020324402705103 | 0.00047195125787112801 |
| 1016 | 0.00012388006107739627 | 8873 | 0.30332053952275101 | 2.3188256905000419 | 0.00079828687214397796 |
| 1627 | 0.01510864856832415 | 2074 | 2.3099278978061428 | 1.3433907369232927 | 0.00063411894348545906 |
| 1828 | 0.00038000715535909145 | 4897 | 0.02101624724387675 | 1.5182822759491885 | 0.00044552767485105132 |


Resposta geral: INCONCLUSIVO. A concentração varia entre os sete casos; 90% da norma não fica restrita a uma ou poucas células. O caso 1627 tem maior IPR entre os sete em Psi_F, mas a célula de maior amplitude está a distância 2.3099278978061428 do ponto estacionário mais próximo. Não há limiar prescrito para classificar IPR alto, proximidade ou degenerescência. Portanto, os dados não autorizam atribuir todos os casos a uma causa geométrica ou à proximidade de pontos estacionários.

Razão de aspecto = diâmetro do polígono / menor largura entre retas de suporte paralelas às suas arestas. Para pares complexos, IPR usa |v|² = Re(v)² + Im(v)². As linhas de maior norma e as dez maiores componentes estão em d4_topcells.csv.

| mesh_family | state | group | cases | IPR_p50 | n90_fraction_p50 |
| --- | --- | --- | --- | --- | --- |
| cartesian | PsiF | selected_ids | 7 | 0.00012012281742235224 | 0.41135110294117649 |
| cartesian | Psi0 | selected_ids | 7 | 0.00012012526486619048 | 0.41135110294117649 |
| hexagonal | PsiF | selected_ids | 7 | 0.00010364553999351139 | 0.41435360837094132 |
| hexagonal | Psi0 | selected_ids | 7 | 0.00010364665262023161 | 0.41435360837094132 |
| random | PsiF | selected_ids | 7 | 0.00041631254143867401 | 0.2034466911764706 |
| random | PsiF | random_controls | 20 | 0.00012497170348891917 | 0.40955882352941175 |
| random | Psi0 | selected_ids | 7 | 0.0010242255649418148 | 0.20165441176470589 |
| random | Psi0 | random_controls | 20 | 0.0001249728503122291 | 0.40955882352941175 |

## Tempos, integridade e reprodução

GNU C++ 15.2.0; PETSc 3.24.4; SLEPc 3.24.2; ADOL-C 2.7.2 (pacote 2.7.2-9); CGAL 6.1.1; VoronoiMeshMaker 1.1.0 e VoroGradToolkit 1.0.0 conforme metadados CMake instalados. Os hashes identificam o conteúdo instalado independentemente desses rótulos. Release -O3 -DNDEBUG, mesmas flags históricas, OMP_NUM_THREADS=1.

| diagnóstico | segundos | exit_code | arquivos_protegidos_alterados |
| --- | --- | --- | --- |
| D0 | 5.5256970000000001 | 0 | 0 |
| D1 | 393.2990246249974 | 0 | 0 |
| D2 | 69.451328171999194 | 0 | 0 |
| D3 | 69.427601384006266 | 0 | 0 |
| D4 | 12.261280360005912 | 0 | 0 |


Os tempos incluem execução e conferência final dos hashes; não incluem compilação. D0 usa timestamps UTC, com início registrado em segundos inteiros.

O repositório já tinha alterações alheias ao diagnóstico. Os manifests preservam o estado completo da árvore; os commits desta entrega contêm apenas diagnostics/. Nenhum push/tag foi feito. Versões, CPU, flags, origem dos casos e hashes dos executáveis constam dos manifests. provenance_verified.json complementa hashes das bibliotecas e árvores de headers, explicitamente capturados após a execução.

### D0

Comando, a partir de diagnostics/: `bash D0_gate/run.sh`.

Execução: `D0_gate/Result/run_20261013T142042Z`.

Arquivos: `gate.csv`, `manifest.txt`, `console.log`, seleção e hashes.

### D1

Comando, a partir de diagnostics/: `bash D1_spectral_hypotheses/run.sh`.

Execução: `D1_spectral_hypotheses/Result/run_20261013T142330Z`.

Arquivos: `d1_cases.csv`, `d1_eigenpairs.csv`, `d1_summary.csv`, `manifest.txt`, `console.log`, seleção e hashes.

### D2

Comando, a partir de diagnostics/: `bash D2_frozen_vs_total/run.sh`.

Execução: `D2_frozen_vs_total/Result/run_20261013T183002Z`.

Arquivos: `d2_directional.csv`, `d2_directional_best.csv`, `d2_floor.csv`, `d2_maximisers.csv`, `d2_spectrum.csv`, `d2_summary.csv`, `d2_w.csv`, `manifest.txt`, `console.log`, seleção e hashes.

### D3

Comando, a partir de diagnostics/: `bash D3_hex_floor_location/run.sh`.

Execução: `D3_hex_floor_location/Result/run_20261013T183153Z`.

Arquivos: `coordinates_L1.csv`, `coordinates_L3.csv`, `coordinates_L4.csv`, `coordinates_main.csv`, `counts_L1.csv`, `counts_L3.csv`, `counts_L4.csv`, `counts_main.csv`, `d3_case_counts.csv`, `d3_floor_cells.csv`, `d3_mesh_symmetry.csv`, `d3_summary.csv`, `floor_L1.csv`, `floor_L3.csv`, `floor_L4.csv`, `floor_main.csv`, `manifest.txt`, `console.log`, seleção e hashes.

### D4

Comando, a partir de diagnostics/: `bash D4_random_eigvec_localisation/run.sh`.

Execução: `D4_random_eigvec_localisation/Result/run_20261013T183140Z`.

Arquivos: `d4_cases.csv`, `d4_summary.csv`, `d4_topcells.csv`, `manifest.txt`, `console.log`, seleção e hashes.


Para recriar as tabelas e este relatório após os cálculos: `python3 support/summarise.py`; para registrar os hashes complementares: `python3 support/provenance.py`. Os CSV têm 17 algarismos significativos. NA significa indisponível/não aplicável, nunca zero.

## Limitações

D1–D4 avaliam os estados e amostras definidos; não comprovam propriedades globais nem reconstroem trajetórias não lineares. A correção espectral de D2 ficou abaixo da resolução prática dos cálculos. D4 não tem limiares classificatórios fornecidos. D5 não foi implementado nem executado; o launcher permanece bloqueado. Não se propõem alterações ao manuscrito ou ao algoritmo.
