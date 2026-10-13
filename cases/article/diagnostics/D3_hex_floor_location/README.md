# D3 — localização do piso

Resposta à pergunta solicitada: SIM, nos estados e casos avaliados. Todas as células com piso estão na geradora da origem das hexagonais L1/L3; a geradora não existe nas hexagonais principal/L4. Resultados exatos em ../RELATORIO_DIAGNOSTICOS.md.

São reconstruídas as 12 malhas, com seed 12345 e dimensões das campanhas. PsiF e PsiF+1e-11 s/||s||inf usam a função de perturbação de produção. Incluem-se todos os 64 casos L1 e 74 L3 com ativação histórica do piso, e dez controlos por família/nível (seed 20261013). Controlo histórico sem piso não implica ausência de piso no estado de referência.

A distância de reflexão é obtida via scipy.spatial.cKDTree sobre todas as geradoras, consultas (-x,-y), workers=1. Origem: distância <1e-12. Piso: g<=delta² max(g). Distância estacionária: mínimo para (k*pi*cos(theta),-k*pi*sin(theta)), k=-4,...,4, consistente com a rotação de Fadeev. Norma de linha é euclidiana, com escala c do CSV.

CSV d3_mesh_symmetry: contagem, presença de origem e distâncias geométricas; d3_floor_cells: geradora, posição, distância estacionária, k, r_P=sqrt(g/M2) e norma da linha; d3_case_counts: estados inclusive os sem piso; d3_summary: agregados por nível/malha/estado. coordinates_*.csv preservam todas as posições; selection.json associa os casos aos arquivos históricos. Índices iniciam em zero. NA em estatística de conjunto vazio.

Reutiliza-se por inclusão support/mesh_repair/MeshFactory.cpp, exatamente a implementação das campanhas de refinamento concluídas. O fallback foi necessário para random L3/L4; está registrado em console.log. Não se alteram bibliotecas instaladas ou tolerâncias. principal/L4 avaliam somente geometria; não há fita AD, fatorização ou estado não linear nessas duas etapas.

Reprodução: bash D3_hex_floor_location/run.sh; depois python3 support/summarise.py. Requer scipy/numpy para as consultas e tabelas.
