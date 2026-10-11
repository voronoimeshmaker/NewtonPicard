# Casos do artigo — preparação para execução pelo autor

**Não executados.** Fontes históricos conferidos e organizados em 24 pastas.
Leia `AUDIT.md` antes de interpretar resultados e `cases.json` para o inventário.
Os casos Fadeev simples em `../fadeev` permanecem intactos.

No Ubuntu-26.04-Test:

```sh
cd /home/jflavio/Programas/NewtonPicard/cases/article
bash run_all.sh --list            # somente lista; não compila nem executa
bash run_all.sh --build-only      # somente compila
bash run_all.sh                   # compila e executa TODOS, sequencialmente
```

Para executar somente um:

```sh
bash run_all.sh --case 02_newton
# ou
bash 02_newton/run.sh
```

Ordem: verificação espacial; Newton, Picard, Initial-rho, Updated-rho; amplitudes;
refinamento dos ramos de Picard. Cada campanha contém várias amostras e malhas.
Os níveis finos são grandes e podem exigir muito tempo e memória. Não foi medida
uma estimativa nesta preparação. O script não usa limite de tempo artificial nem
altera o tamanho dos casos. BUILD_JOBS controla apenas a compilação (padrão 2).
A execução usa uma thread por biblioteca, sem campanhas em paralelo.

Cada caso lê seu config.yaml e grava exclusivamente em sua pasta:
- Result/run_*: CSVs, manifesto e console do driver histórico;
- logs/: compilação e saída do lançamento;
- tmp/: temporários solicitados via TMPDIR.
A infraestrutura compartilhada de compilação fica em build/ dentro deste pacote.
O script para no primeiro erro de processo. Reexecutá-lo inicia novas corridas;
não há retomada automática nem sobrescrita intencional dos resultados antigos.
Nenhum caminho de execução depende do repositório histórico de origem.

Os programas não devem ser lançados com gss_solver. Cada pasta tem seu main.cpp;
run_all.sh usa o executável específico. Os quatro YAMLs principais têm chave
solver.nonlinear histórica, mas a estratégia efetiva é fixada no main.

Verificação estática reproduzível (sem solver):

```sh
python3 audit_static.py
```

Após executar, verifique os manifestos e estados de cada amostra. Picard pode
convergir para outro ramo ou falhar em parte da população; esses são resultados
estudados no artigo. Não confundir término do processo com recuperação da referência.

O launcher fixa a execução em um CPU permitido via taskset (ARTICLE_CPU pode
selecionar outro). Requer Python 3 com PyYAML, flock e taskset. Antes de compilar,
confere a integridade das cópias e parâmetros do protocolo. Há trava contra dois
launchers simultâneos neste pacote.

