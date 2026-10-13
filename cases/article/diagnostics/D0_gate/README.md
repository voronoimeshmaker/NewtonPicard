# D0 — gate de reprodução

APROVADO: 60/60, erro relativo máximo zero. Vinte primeiros case_id WLSQR com rho_reference finito positivo, por família; random exclui ID3 indisponível e inclui ID21. Fonte: 03_picard/Result/run_20261009T191103Z/cases.csv. Usa c gravado e PsiF, sem recalibração. Tolerância |rho_new-rho_csv|/rho_csv <=1e-10.

Inclusão direta do código de produção, mesmas bibliotecas instaladas, Release, serial OMP_NUM_THREADS=1. CSV gate.csv: família, ID, c, rho histórico/reproduzido, erro relativo, convergência EPS, resíduo, tempo e pass. Os casos e hashes estão no manifest e selected.tsv. Em caso de falha o launcher devolve status não zero; os demais launchers verificam a aprovação dos 60 casos.

Reprodução, a partir de diagnostics/: bash D0_gate/run.sh.
