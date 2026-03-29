# Resolução: Métodos Iterativos e Esparsos

## Exercícios do Livro

  - $\mathbf{r}_0 = \mathbf{b} - \mathbf{Ax}_0, \quad \mathbf{p}_0 = \mathbf{r}_0, \quad k=0$
  - Enquanto $\|\mathbf{r}_k\| > \vepsiilon$:
    \begin{itemize}
      - $\alpha_k = \frac{\mathbf{r}_k^T \mathbf{r}_k}{\mathbf{p}_k^T \mathbf{A p}_k}$
      - $\mathbf{x}_{k+1} = \mathbf{x}_k + \alpha_k \mathbf{p}_k$
      - $\mathbf{r}_{k+1} = \mathbf{r}_k - \alpha_k \mathbf{Ap}_k$
      - $\beta_k = \frac{\mathbf{r}_{k+1}^T \mathbf{r}_{k+1}}{\mathbf{r}_k^T \mathbf{r}_k}$
      - $\mathbf{p}_{k+1} = \mathbf{r}_{k+1} + \beta_k \mathbf{p}_k$
      - $k = k + 1$
    \end{itemize}
