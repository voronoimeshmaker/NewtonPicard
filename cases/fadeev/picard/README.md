# Fadeev — Picard

Execute no Ubuntu/WSL:

```sh
bash /home/jflavio/Programas/NewtonPicard/cases/fadeev/picard/run.sh
```

O script compila e executa este caso. O arquivo `config.yaml` é lido pelo
executável `fadeev_picard_case`. As saídas são gravadas em `results/`:

- `config_used.yaml`: configuração utilizada.
- `convergence.csv`: erros L1, L2, Linf, resíduo final, ordem, iterações e tempo.
- `solution_N.csv`: coordenadas, solução numérica, analítica, erro e resíduo por célula.
- `residual_history_N.csv`: histórico de convergência.

As malhas são definidas por `convergence_study.n_generators_sequence`: números
quadrados crescentes. O domínio e os parâmetros analíticos são fixos neste caso
e validados pelo executável. A solução inicial é zero e diagonal_shift=0.
O caso retorna erro se não convergir, se o resíduo Linf superar 1e-10 ou se a
ordem L2 entre malhas for inferior a 1.8. A primeira malha tem ordem 0 por convenção.

Uma nova execução substitui os CSVs de mesmo nome. O histórico de Newton inclui
a avaliação inicial; o de Picard registra os resíduos após as atualizações.
Consulte `../README.md` para a equação, condições de contorno e limitações.

Todas as saídas ficam exclusivamente em `results/` nesta pasta, mesmo se o
executável for iniciado em outro diretório. Destinos externos são rejeitados.
