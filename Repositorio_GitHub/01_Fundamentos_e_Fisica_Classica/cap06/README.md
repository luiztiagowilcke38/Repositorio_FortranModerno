# Resolução: Métodos Diretos para Sistemas Lineares

## Exercícios do Livro

  -  Troca-se a linha $k$ pela linha $p$ onde $p = \arg\max_{i \geq k} |A_{ik}|$;
  -  Para cada linha $i > k$, calcula-se o multiplicador $m_{ik} = A_{ik}/A_{kk}$ e atualiza-se:
    $$A_{ij} \leftarrow A_{ij} - m_{ik} A_{kj}, \qquad b_i \leftarrow b_i - m_{ik} b_k, \quad j = k+1,\ldots,n.$$
