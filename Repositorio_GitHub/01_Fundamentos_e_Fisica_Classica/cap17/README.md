# Resolução: EDPs: Volumes Finitos e CFD

## Exercícios do Livro

    -  Resolve a equação de momento usando a pressão do passo anterior $P^*$ para obter velocidades intermediárias $u^*, v^*$.
    -  Deriva-se uma equação de Poisson para a correção $P'$ baseada na divergência de $\mathbf{v}^*$:
    \begin{equation}
    \nabla^2 P' = \frac{\rho}{\Delta t} (\nabla \cdot \mathbf{v}^*)
    \end{equation}
    -  Atualiza-se a pressão $P = P^* + \alpha_P P'$ e a velocidade $\mathbf{v} = \mathbf{v}^* - \frac{\Delta t}{\rho} \nabla P'$.
