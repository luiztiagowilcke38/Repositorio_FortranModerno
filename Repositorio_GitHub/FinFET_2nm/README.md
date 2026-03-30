# FinFET 2nm Physical Simulator (NEGF + Poisson-Schrödinger)

Este diretório contém a implementação completa de um simulador físico para um nanotransistor FinFET de 2 nm, desenvolvido em Fortran Moderno (2008/2018). Este código acompanha o Capítulo correspondente no livro "Fortran Moderno para Computação Científica: Uma Enciclopédia Doctoral".

## Estrutura do Projeto

- **mod_constantes.f90**: Parâmetros físicos universais e dimensões geométricas do FinFET (escala de 2 nm).
- **mod_malha.f90**: Geração e gerenciamento da malha de discretização 3D.
- **mod_poisson.f90**: Solver iterativo (SOR) para a equação de Poisson 3D.
- **mod_schrodinger.f90**: Solver de Schrödinger 2D na seção transversal, com integração LAPACK (DSYEV).
- **mod_negf.f90**: Implementação do formalismo de Funções de Green fora do equilíbrio, incluindo auto-energias (Sancho-Rubio) e transmissão.
- **main_finfet.f90**: Programa principal que orquestra o loop autoconsistente (SCF) e o transporte quântico.

## Requisitos

- Compilador Fortran compatível com Fortran 2008/2018 (ex: `gfortran`).
- Bibliotecas **LAPACK** e **BLAS** instaladas no sistema.

## Como Compilar e Executar

Sugestão de comando para compilação com `gfortran`:

```bash
gfortran -O3 mod_constantes.f90 mod_malha.f90 mod_poisson.f90 mod_schrodinger.f90 mod_negf.f90 main_finfet.f90 -o simulador_finfet -llapack -lblas
```

Para executar:

```bash
./simulador_finfet
```

## Descrição Física

O simulador resolve o acoplamento entre a eletrostática (Poisson) e a mecânica quântica (Schrödinger/NEGF) para capturar efeitos de confinamento espacial e transporte balístico em dispositivos ultra-escalonados. Ideal para estudos de doutorado em nanoeletrônica computacional.
