# D2 — normalização congelada e derivada total

Resultado: os dois critérios foram satisfeitos na amostra; detalhes e números em ../RELATORIO_DIAGNOSTICOS.md e Result/LATEST.txt.

164 casos: união de 30 estratificados por malha com 46 recuperações obrigatórias, sete IDs aleatórios e dez casos próximos do limiar por malha. Seed 20261013; a seleção e a origem de cada parâmetro estão em selection.json. Os bins com poucos candidatos são esgotados e o grupo básico é completado com candidatos restantes. Avaliam-se PsiF e a perturbação de produção de amplitude 1e-4.

D_fr é calculado pela classe copiada de produção. A única mudança na classe é o nome e a exposição dos membros. A direção AD e_n fornece z = db/dtau. Cada maximizador fornece w = 2 delta² Gamma_Q^T grad(Psi)_Q. O shell aplica A^-1 D_fr x + A^-1 z (w^T x). A previsão de primeira ordem usa os vetores esquerdo/direito do autovalor real; é NA para par complexo. Não se forma matriz densa.

Os arquivos d2_floor.csv, d2_directional.csv e d2_spectrum.csv têm exatamente os campos identificados nos cabeçalhos. Coordenadas e índices referem-se às geradoras, índices iniciando em zero. d2_maximisers.csv lista todos os empates, coordenadas dos primeiros dez; d2_w.csv armazena entradas não nulas de cada w (entradas omitidas são zero). Empate significa M2-g <= 1e-12 M2. Cada candidato recebe seu próprio espectro; não se pressupõe derivada única. Para FD central usa-se a média das derivadas unilaterais, (max wv + min wv)/2.

A fonte de produção é reavaliada em ambos os lados da diferença central. Os quatro passos relativos são 1e-4, 1e-5, 1e-6 e 1e-7; vetores têm norma euclidiana um; random usa mt19937_64, seed 20261013+case_id, uniforme [-1,1]. O mínimo erro por direção dentre esses quatro passos é registrado em d2_directional_best.csv, sem alterar a janela. d2_summary.csv fornece máximos e percentis lineares.

Aceitação: max |drho| < 2.9e-4 e < 0.1 min |rho_ref-1|, mais melhor err_tot <=1e-6 na janela em estados sem piso. As variações espectrais observadas são da ordem do arredondamento; seu sinal não é resolvido por esta experiência.

Reprodução, a partir de diagnostics/: bash D2_frozen_vs_total/run.sh. Depois: python3 support/summarise.py. Não há integração não linear.
