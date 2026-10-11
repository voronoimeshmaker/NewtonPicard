# Jacobiano de Newton com ADOL-C

O Newton utiliza `ADOLCJacobian`, selecionado por padrão com `newton.jacobian: adolc`.
A instalação usada neste ambiente é ADOL-C 2.7.2, disponível em `/usr`.
Adept não é uma dependência do projeto. Não há alternativa por diferenças finitas.
Configurações antigas com `jacobian: finite_difference` são rejeitadas ao construir Newton.
`fd_epsilon` foi removido da configuração C++ e não é utilizado.

O resíduo registrado é `F(psi) = (A + diagonal_shift I) psi - B(psi)`.
O stencil e a fonte aceitam escalares ativos `adouble`. A ADOL-C registra novamente
as operações a cada estado de Newton e `sparse_jac` calcula as entradas esparsas
com ColPack. Cada processo insere no PETSc apenas as linhas que possui.
A gravação ocorre fora das regiões OpenMP. Esta implementação requer PETSc com
escalares reais de precisão dupla e ADOL-C com suporte a ColPack.

A fonte padrão Fadeev é prescrita em função das coordenadas e independe de psi;
para ela, J = A + diagonal_shift I. A integração não altera a equação física. Fontes não lineares
personalizadas devem manter o tipo escalar genérico durante todo o cálculo,
sem converter `adouble` para `double`.

## Compilar e testar no Ubuntu/WSL

```sh
cd /home/jflavio/Programas/NewtonPicard
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel 4
OMP_NUM_THREADS=1 ctest --test-dir build -R '^test_adolc_jacobian$' --output-on-failure
```

Os testes em `tests/test_adolc_jacobian.cpp` verificam:

1. Fonte prescrita: o Jacobiano ADOL-C coincide com a matriz de difusão, erro < 1e-12.
2. Fonte `S_i = psi_i^3 + 0.25 sin(psi_(i+1))`: todas as entradas coincidem com
   `A_ij - V_i[3 psi_i^2 delta_ij + 0.25 cos(psi_(i+1)) delta_(j,i+1)]`,
   erro < 1e-12 em dois estados distintos. Verifica também a liberação dos tapes.
3. Newton com fonte `S_i = 1 - psi_i^3`: converge de psi=0.5 para psi=1,
   erro de solução < 1e-8.
4. A fábrica aceita ADOL-C e rejeita a antiga seleção de diferenças finitas.

Verificação distribuída do Jacobiano (dois processos):

```sh
OMP_NUM_THREADS=1 mpiexec -n 2 build/tests/test_adolc_jacobian \
  --gtest_filter='ADOLCJacobian.PrescribedSourceMatchesDiffusionMatrix:ADOLCJacobian.NonlinearCoupledSourceMatchesAnalyticDerivativeAfterRetaping'
```

O teste distribuído valida a montagem do Jacobiano; a convergência de Newton é
validada em execução serial. Os testes específicos não precisam de arquivo YAML.
O arquivo `config/default.yaml` foi restaurado a partir do projeto VoronoiGradShafranovSolver,
com malha de 16 x 16 células e `newton.jacobian: adolc`. A configuração padrão
executa Picard com o deslocamento diagonal de teste já previsto no programa.
As 14 entradas do CTest passaram; a execução padrão convergiu em uma iteração.
A validação de Newton com ADOL-C permanece no teste não linear descrito acima.

Para executar o programa a partir da pasta do projeto:

```sh
./build/gss_solver
```
## Revisão adicional

O teste `IndependentThreeCellResidualAndJacobian` usa uma malha de três células,
volumes diferentes e a matriz explícita [[2,-2,0],[-2,5,-3],[0,-3,3]], sem
usar a montagem de A como referência. Compara o resíduo e todas as derivadas
analíticas em dois estados, incluindo valores negativos e zero. Erro máximo
observado em execução serial: 4.44089e-16. Também passou com dois processos MPI.
Nenhuma diferença finita é usada nesse teste ou no Jacobiano de produção.

A inconsistência do deslocamento diagonal foi corrigida: Newton agora aplica
`solver.diagonal_shift` no resíduo, na busca de passo e no registro ADOL-C.
O Jacobiano inclui a derivada desse termo automaticamente. Dois testes adicionais
verificam o deslocamento da diagonal e a igualdade entre as soluções de Newton
e Picard, além do resíduo final. O caso cúbico não deslocado é mantido explicitamente
com `diagonal_shift = 0`.

Para executar Newton com ADOL-C:

```sh
./build/gss_solver config/newton.yaml
```

A configuração fornecida convergiu com `iterations=2` (uma atualização e a
verificação seguinte). A suíte completa passou. Quando não há convergência,
o executável agora informa o problema e retorna código 2; validado com max_iter=0.
O deslocamento mantém a regularização já usada pelo Picard; não substitui a
especificação das condições de contorno de um problema físico sem regularização.