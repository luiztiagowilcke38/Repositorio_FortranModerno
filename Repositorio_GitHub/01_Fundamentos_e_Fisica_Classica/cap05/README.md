# Resolução: Álgebra Linear Fundamental em Fortran

## Exercícios do Livro
 -  Implemente uma subrotina `produto\_matvec` que compute $\mathbf{y} = \mathbf{A}\mathbf{x}$ sem usar `MATMUL`, apenas via loops duplos. Compare o desempenho com `MATMUL` e com `DGEMV` para $n = 100, 500, 1000$. Use `CPU\_TIME` e variáveis `tempo\_inicio`, `tempo\_fim`.

   -  Prove que $\|\mathbf{A}\|_2 = \sigma_{\max}(\mathbf{A})$ onde $\sigma_{\max}$ é o maior valor singular. Utilize a definição da norma induzida $\|\mathbf{A}\|_2 = \max_{\|\mathbf{x}\|_2=1} \|\mathbf{A}\mathbf{x}\|_2$.

   -  Escreva um programa que gere uma matriz de Vandermonde $V_{ij} = x_j^{i-1}$ para $n$ pontos igualmente espaçados em $[0,1]$, calcule suas normas $\|\cdot\|_1$, $\|\cdot\|_\infty$, $\|\cdot\|_F$ (implementadas em Fortran) e estime o número de condição como $\|V\|_F \|V^{-1}\|_F$ (usando DGESV para inverter).

  