# D4 — localização dos vetores dominantes

41 casos, 82 estados: sete IDs em cada família e vinte controlos random sem anomalia, seed 20261013. PsiF/Psi0 usam os mesmos estados de D2. Resposta geral: INCONCLUSIVO para atribuir uma causa única a todos os casos; a concentração varia. Comparações quantitativas estão em ../RELATORIO_DIAGNOSTICOS.md.

IPR=sum |v_i|^4/(sum |v_i|²)^2; n90 conta as maiores componentes até 90% da norma ao quadrado. Para autovalor complexo considera-se |v_i|²=Re(v_i)²+Im(v_i)². O espectro bilateral usa o mesmo shell de D1.

Razão de aspecto: diâmetro do polígono / menor largura entre retas de suporte paralelas às arestas. Área é o volume 2D usado na discretização. As coordenadas são das geradoras; distância estacionária e r_P seguem D3. Norma da linha é euclidiana do D_fr integrado, incluindo c. Índices iniciam em zero.

CSV d4_cases: espectro dominante, IPR, n90 e convergência EPS; d4_topcells: dez maiores componentes mais a célula da maior norma de linha, com geometria e métricas; d4_summary: percentis por família, estado e grupo (IDs selecionados/controlos random). NA indicaria razão de aspecto indefinida por largura não positiva. Não há limiares prescritos para chamar uma célula degenerada ou um modo localizado; não se inventam tais limiares.

Reprodução: bash D4_random_eigvec_localisation/run.sh; depois python3 support/summarise.py.
