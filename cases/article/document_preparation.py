from pathlib import Path
import json,hashlib,difflib,subprocess
p=Path('/home/jflavio/Programas/NewtonPicard/cases/article')
manifest=json.loads((p/'provenance.json').read_text())
patches=[]
for item in manifest:
 src=Path(item['source']); dest=p/item['destination']
 assert hashlib.sha256(src.read_bytes()).hexdigest()==item['source_sha256'],str(src)
 item['prepared_sha256']=hashlib.sha256(dest.read_bytes()).hexdigest()
 if item['prepared_sha256']!=item['source_sha256']:
  patches.extend(difflib.unified_diff(src.read_text().splitlines(True),dest.read_text().splitlines(True),fromfile=str(src),tofile=item['destination']))
(p/'provenance.json').write_text(json.dumps(manifest,indent=2))
(p/'preparation.patch').write_text(''.join(patches))
(p/'AUDIT.md').write_text('''# Conferência dos casos do artigo

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
''')
(p/'README.md').write_text('''# Casos do artigo — preparação para execução pelo autor

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
''')
(p/'audit_static.py').write_text('''from pathlib import Path
import hashlib,json,subprocess
import yaml
root=Path(__file__).resolve().parent
records=json.loads((root/'provenance.json').read_text())
for rec in records:
 p=root/rec['destination']
 assert hashlib.sha256(p.read_bytes()).hexdigest()==rec['prepared_sha256'], str(p)
cases=json.loads((root/'cases.json').read_text())
assert len(cases)==24
for case in cases:
 d=root/case['name']; cfg=yaml.safe_load((d/'config.yaml').read_text())
 assert cfg['mesh']['domain']==[-15.,15.,-5.,5.]
 subprocess.run(['bash','-n',str(d/'run.sh')],check=True)
 if case['kind']!='spatial':
  assert cfg['source']['derivative']=='adolc'
  assert cfg['sampling']['seed']==20260914
  assert cfg['solver']['diagonal_shift']==0
  assert cfg['solver']['max_iter']==5000
  assert cfg['spectral']['requested_eigenpairs']==8
  assert cfg['krylov']['relative_tolerance']==0.01
  assert cfg['krylov']['max_iterations']==50
  expected=2000 if case['kind'].startswith('official:') else (100 if case['name'].startswith('branch_') else 200)
  assert cfg['sampling']['sample_count']==expected
  assert cfg['sampling']['rho_strata_edges']==[.70,.90,.98,1.02,1.10,1.30]
 assert 'current_path(ARTICLE_CASE_DIR)' in (d/'main.cpp').read_text()
subprocess.run(['bash','-n',str(root/'run_all.sh')],check=True)
print(f'Static checks passed: {len(records)} copied sources and {len(cases)} cases. No numerical program executed.')
''')
print('Audit documentation and integrity checks prepared')