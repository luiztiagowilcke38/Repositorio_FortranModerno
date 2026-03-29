# Resolução: Transformers e Mecanismos de Atenção em Fortran

## Exercícios do Livro

    -  Por que o fator de escala $\sqrt{d_k}$ é necessário para a estabilidade do gradiente no softmax?
    -  Explique o conceito de "Multi-Head Attention" e como ele permite focar em diferentes partes da sequência simultaneamente.
    -  Converta a rotina de atenção para utilizar operações `do concurrent` visando aceleração paralela.
    -  Implemente o codificador posicional (\textit{Positional Encoding}) usando funções seno e cosseno.
    -  Verifique a complexidade computacional da atenção em função do comprimento da sequência $N$.
    -  Simule o processamento de uma pequena frase representada por vetores latentes (\textit{embeddings}).
    -  Como bibliotecas como a `stdlib` do Fortran facilitam a implementação de tensores e álgebra linear de IA?
