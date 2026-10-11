# Conferência dos casos do artigo

## Estado e escopo
Preparação baseada em leitura dos fontes LaTeX do ZIP, revisão estática de código,
comparação dos arquivos copiados e compilação. **Nenhum caso, campanha ou teste
numérico foi executado nesta preparação.** Compilar não comprova reprodução dos
números do artigo. Os diretórios Result ainda devem ser gerados pelo autor.

## Origem e decisão de seleção
A pasta indicada VoronoiGradShafranovSolver contém antecessores de spectral_study
e codex_spectral_study; alguns usam diferenças finitas. Não foram selecionados
como se produzissem as campanhas finais. A pesquisa pelas sementes e parâmetros
encontrou a implementação correspondente em:
`/home/jflavio/Programas/GSMaker/Old/GSMakerLib/papers/PaperIII`.
A versão sem New é a referência histórica designada no repositório.
Os quatro cabeçalhos *10000 são idênticos (SHA256 começa com 16b1cfc812b37d59).
O nome 10000 é histórico: as configurações efetivas usam **2000 amostras**.
Correspondência de código e protocolo foi confirmada; não foi estabelecida uma
cadeia completa de hashes ligando cada célula de tabela do ZIP a um CSV histórico.

## Itens conferidos
| Artigo | Implementação preservada |
|---|---|
| Fadeev rotacionada em [-15,15] x [-5,5] | Fadeev.cpp e configurações |
| Fonte pela razão entre gradientes reconstruídos | ModifiedSource + evaluate_base_source; WLSQR/MGG externos |
| Família discreta q_c=A psi_F-b_D-c b(psi_F) | make_compensator em official.hpp/complementary.hpp |
| Calibração c=rho_target/rho_c1 por amostra | calibrate_case_scales, Cartesian/WLSQR no estado de referência |
| 2000 amostras, cinco estratos, seed 20260914 | mt19937_64; estrato case_id módulo cinco |
| Perturbação inicial 1e-4 | configuração principal e forma determinística histórica |
| ADOL-C, sem fallback automático | trace_on + fov_forward com coloração; source.derivative=adolc |
| Piso de regularização congelado | entrada adicional do tape com sementes zero; coerente com seção 2 do artigo |
| J=A-D | PreconditionedNewtonSolver atualiza a matriz com os valores de D |
| Newton inexato | GMRES, rtol=1e-2, máximo 50, passo completo; fatoração de A reutilizada |
| Critério não linear | atol+rtol*max(norm(R0),1), 5000 passos; divergência fator 1e8 |
| SLEPc | 8 autovalores, tolerância 1e-11, máximo 3000 |
| Referência alcançada | distância relativa <=1e-6, distinta da convergência de resíduo |
| Initial-rho / Updated-rho | avaliação inicial / antes de cada passo, limiar 1 |
| Verificação espacial | fonte analítica fixa, não a fonte reconstruída; 48x16,96x32,192x64 |

## Casos selecionados
- 01: tabela de convergência espacial, três famílias.
- 02–05: quatro campanhas principais, 2000 amostras por família e reconstrução.
  O driver histórico também executa MGG; as tabelas principais do artigo usam WLSQR.
- amplitude_*: quatro amplitudes citadas (1e-4,1e-2,1e-1,1), 200 amostras,
  quatro estratégias, WLSQR, três famílias. Não incluímos 0.5, ausente no texto.
- branch_L1/L3/L4_picard: níveis reportados 128x43,512x171,1024x341,
  100 amostras por família. O L2 histórico 256x85 não aparece na tabela de ramos.
  O nível fino tem 349184 células nominais; o número hexagonal é medido pelo VMM.
  Nos níveis L3/L4, os diagnósticos de rho de referência/inicial desabilitados
  no protocolo complementar original permanecem desabilitados. NÃO usar os bins
  por alvo desses arquivos como se fossem bins por rho_ref. Para a tabela de ramos,
  comparar as distâncias finais na interseção dos mesmos case_id entre níveis.

As projeções many-query são pós-processamento dos tempos das campanhas principais,
não execuções adicionais de 10^4,10^5 ou 10^6 sistemas. Estatísticas do artigo não
são copiadas como novos resultados. Exceções técnicas e convergência para outra
solução devem ser mantidas nas contagens apropriadas, nunca descartadas como sucesso.

## Alterações de preparação
- Includes relativos realocados para support/common; nenhum algoritmo alterado.
- Entradas main fixam o diretório de trabalho no caso; fontes históricas produzem
  Result/run_* dentro dele. Wrapper chama o driver histórico correto por estratégia.
- CMake independente compila uma cópia isolada da biblioteca histórica; não mistura
  seus tipos/solvers com GSSLib simplificada do NewtonPicard e não toca os originais.
- Compatibilidade de PetscErrorMessage: argumento de mensagem específica omitido
  (nullptr), mantendo a mensagem genérica. Não altera operações numéricas.
- VMM instalado exporta nomes modulares; target de compatibilidade liga a biblioteca
  legada também instalada e os alvos modulares. É uma adaptação de link, não substituição
  de geometria. A equivalência numérica da versão externa só pode ser verificada ao rodar.

`provenance.json` registra os hashes originais e preparados; `preparation.patch`
registra diferenças nos fontes copiados. Configurações principais e complementares
foram mantidas. `solver.nonlinear: picard` nos YAMLs históricos NÃO escolhe a estratégia
nesses programas: a seleção é o enum passado pelo main; preservado e documentado.

## Dependências e limitações
A compilação usa o par PETSc/SLEPc do mesmo ambiente pkg-config. Nesta preparação:
PETSc 3.24.4, SLEPc 3.24.2, ADOL-C 2.7.2. O artigo cita PETSc 3.25.99;
essa diferença impede prometer identidade binária, tempos ou reprodução bit a bit.
O backend Newton dos casos do artigo usa o driver histórico com GMRES e fonte
reconstruída; não é o Newton simplificado usado nos primeiros testes Fadeev.
As fórmulas, sementes e tolerâncias foram preservadas em vez de adaptar as tabelas
a um modelo diferente. Não há resultado numérico novo nem afirmação de validação dinâmica.

O código histórico retorna zero quando não há falha de setup, mesmo havendo casos
que não convergem ou exceções individuais. Isso faz parte do desenho das campanhas:
consulte cases.csv e manifest.txt; código zero não significa 100% de sucesso científico.


## Correspondência com resultados históricos existentes (somente leitura)
Foi localizada a saída NewtonOnly10000/Result/run_20260920T052100Z/summary.csv.
A mediana WLSQR/cartesiana é 0.074931204 s, que corresponde a 0.074931 s citada
no artigo; hexagonal 0.1226631045 s e random 0.1913537505 s também correspondem
às tabelas arredondadas. São resultados históricos, não corridas deste pacote.
`reference/historical_summaries.json` registra os caminhos, hashes e linhas WLSQR
das quatro campanhas encontradas. O total histórico 43164 sistemas de Newton
inclui WLSQR e MGG; a soma apenas WLSQR é 21593. Ao pós-processar, explicitar
sempre se o conjunto contém uma ou ambas as reconstruções.
