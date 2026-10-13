# D1 — Hipóteses espectrais
WLSQR em Psi_F. Inclui todos os rho_reference >= 1, exceções técnicas com rho recuperado em 02_newton e 100 controlos por malha (seed 20261013).
Operador de produção com shell adicional para transposição B^T A^-T; EPS Krylov-Schur bilateral, 8 pares, tol 1e-11, maxit 3000. Vetor esquerdo permite kappa para autovalor dominante real. NA se complexo.
CSV cases: partes real/imaginária e módulos dos primeiros oito pares; eigenpairs: todos os convergidos e resíduos relativos. summary: contagens e percentis por malha/bin. Simplicidade é apenas proxy, não prova.
Reprodução: bash D1_spectral_hypotheses/run.sh. Exige D0 aprovado.
