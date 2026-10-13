# Diagnósticos complementares — NewtonPicard

D0–D4 concluídos. Consulte [RELATORIO_DIAGNOSTICOS.md](RELATORIO_DIAGNOSTICOS.md) para os dados e limitações. D5 está desativado e não foi executado.

Ambiente: WSL Ubuntu-26.04-Test; repositório atual /home/jflavio/Papers/NewtonPicard (o caminho antigo /home/jflavio/Programas/NewtonPicard não existe). Branch diagnostics-2026-10, base b8b5f057d68d5ecdcbf9c4a3abe9f6155a26854a. Todos os novos fontes estão nesta pasta. Nenhum resultado histórico, fonte protegido ou biblioteca instalada foi alterado.

- D0_gate: 60 reproduções históricas, requisito obrigatório.
- D1_spectral_hypotheses: 3281 casos; espectros, resíduos e vetores bilaterais.
- D2_frozen_vs_total: 164 casos em dois estados; fita AD, empates, diferenças direcionais e espectros totais.
- D3_hex_floor_location: 12 malhas e casos com piso L1/L3.
- D4_random_eigvec_localisation: 41 casos em dois estados; IPR e geometria.
- D5_branch_alignment: proteção de desativação; não implementa o diagnóstico opcional.

Em diagnostics/, execute os run.sh com bash, começando por D0. Cada execução usa build/ Release separado e Result/run_<UTC>/; LATEST.txt aponta para os resultados. Result/ e build/ são ignorados pelo Git para preservar artefatos locais sem versionar binários. Os commits versionam fontes e relatório. Não houve push nem tags.

Após os diagnósticos, execute `python3 support/summarise.py` e `python3 support/provenance.py`. NumPy/SciPy já estão disponíveis no ambiente. Não é necessário reconstruir campanhas não lineares. Todos os CSV numéricos usam 17 algarismos significativos; NA representa indisponibilidade.

D1 e D2 usam EPS bilateral para vetores esquerdos; tolerância 1e-11, maxit 3000, nev8. D0 reutiliza diretamente a rotina espectral histórica. D3 reutiliza o MeshFactory com reparo exato das campanhas de refinamento. Detalhes e fórmulas em cada README.
