# Resolução: Aceleração GPU e Coarrays

## Exercícios do Livro

    -  Compare as vantagens e desvantagens de usar OpenACC vs CUDA Fortran.
    -  Explique o conceito de "latência de transferência" entre a CPU (Host) e a GPU (Device).
    -  O que são Coarrays "críticos" e como eles garantem a integridade dos dados?
    -  Reaproveite seu solver de propagação de calor e adicione diretivas OpenACC para acelerar o loop principal.
    -  Implemente uma soma de prefixos (\textit{Scan}) usando Coarrays entre 4 imagens.
    -  Meça o \textit{speedup} do seu código ao variar o tamanho do problema $N$ de $10^3$ para $10^7$ elementos em uma GPU.
    -  Como projetos de exaescala (como o supercomputador Frontier) utilizam Fortran com aceleração AMD (ROCm)?
