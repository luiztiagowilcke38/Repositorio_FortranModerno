# Resolução: Aritmética de Ponto Flutuante e Erros Numéricos

## Exercícios do Livro
 -  Derive a expressão $\varepsilon_m = 2^{-(t-1)}$ para o épsilon da máquina de um sistema de ponto flutuante em base 2 com $t$ bits de mantissa (incluindo o bit implícito). Por que na prática $\varepsilon_m = 2 \cdot 2^{-(t-1)} = 2^{-(t-1)}$? Qual a diferença entre \textit{unit roundoff} e \textit{machine epsilon}?

   -  Implemente a soma compensada de Kahan em Fortran para somar $n = 10^6$ termos da série $1/k^2$ e compare o resultado com a soma direta. Use variáveis `soma\_direta`, `soma\_kahan`, `compensacao`, `erro\_arredondamento`. Compare com a solução analítica $\pi^2/6$.

   -  Para a equação quadrática $x^2 - 56x + 1 = 0$, compute as raízes pela fórmula direta e pela fórmula estável (usando a relação $x_1 x_2 = c/a$). Use `real64` e calcule os erros relativos comparando com a raiz exata. Explique o cancelamento catastrófico observado.

  