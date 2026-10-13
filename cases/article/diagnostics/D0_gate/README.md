# D0 — gate de reprodutibilidade

Recalcula rho_ref em Psi_F de 20 casos WLSQR por família, usando c do CSV.
Selecção determinística: primeiros 20 case_id com rho_reference finito e positivo
em cada família; casos excepcionais não são usados para este gate.
Critério: abs(rho_new-rho_csv)/rho_csv <= 1e-10 para todos os 60 casos.
Não resolve o problema não linear. Código de produção incluído sem alterações.
selected.tsv e sources.json identificam os casos e arquivos de origem.
CSV: família, ID, c, rho histórico/novo, erro relativo, nconv, iterações EPS,
resíduo relativo máximo dos autovalores, segundos e aprovação por caso.
Resultado: pendente da execução.
