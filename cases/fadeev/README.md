# Casos de validação com a solução analítica de Fadeev

- `newton/`: Newton, Jacobiano esparso por ADOL-C/ColPack.
- `picard/`: iteração de Picard, matriz fixa com Dirichlet.
- `common/`: implementação compartilhada da equação e geração dos relatórios.

Para epsilon=0.1, psi0=k=1:

    psi_exata(x,y) = log(cosh(y) - 0.1 cos(x)) - log(0.9)
    Laplaciano(psi) = (1.1/0.9) exp(-2 psi)

Domínio [-1,1] x [-1,1], malhas cartesianas uniformes de 8², 16² e 32² células.
Dirichlet analítico nas quatro bordas, avaliado no centro de cada face.
A transmissão na borda usa comprimento da face / distância centro-borda.
Não há deslocamento diagonal artificial. Os solvers partem do campo zero.
É uma fonte não linear avaliada na solução corrente, não uma fonte espacial
prescrita nem um vetor fabricado pela aplicação da matriz ao campo exato.

Os executáveis destes casos usam os solvers do projeto com uma política de
fonte e contorno específica. Não executar esses YAMLs com `gss_solver`: esse
programa ainda usa outra política de fonte. A implementação comum também é
utilizada por `tests/test_fadeev_solvers.cpp`, que compara os dois métodos.

Os parâmetros físicos fixos são validados; os campos genéricos de configuração
que não fazem parte deste benchmark não ativam outros modelos ou saídas.
Os relatórios são de execução serial. L1/L2 são normalizados pelo volume total;
Linf é o maior erro absoluto. A validação cobre este domínio e estes parâmetros.

Para a regressão conjunta:

```sh
cd /home/jflavio/Programas/NewtonPicard
OMP_NUM_THREADS=1 ctest --test-dir build -R '^test_fadeev_solvers$' --output-on-failure
```
