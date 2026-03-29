# Repositório: Fortran Avançado & Biblioteca SuperDiFF v3.0

Este repositório contém o ecossistema completo de códigos para o manuscrito **"Fortran Avançado: Enciclopédia Doctoral de Modelagem e Soluções Científicas"**, de autoria de **Luiz Tiago Wilcke**.

## 🚀 Conteúdo do Repositório

### 1. Biblioteca SuperDiFF (v3.0 Turbo)
Localização: `/superdiff`

A **SuperDiFF** é uma biblioteca de engines numéricas de elite, escrita em **C puro**, projetada para rodar em hardware de alta performance (HPC). 
- **Motores EDO (ODE)**: Runge-Kutta 4 e Runge-Kutta-Fehlberg 4(5) [Adaptativo].
- **Motores EDP (PDE)**: Solvers de Calor 1D e Poisson 2D altamente otimizados com **Cache-Blocking (Tiling)**.
- **Paralelização**: Suporte total a **OpenMP** para multi-core.
- **Vetorização SIMD**: Otimização manual via **AVX/AVX-512** para processamento vetorial paralelo.
- **Cálculo Fracionário e Estocástico**: Motores para derivadas fracionárias de Grünwald-Letnikov e Equações Diferenciais Estocásticas (Euler-Maruyama).

### 2. Manuscrito e Códigos Fortran
Localização: `/livro`

Contém os fontes LaTeX (`main.tex`) e os códigos Fortran inovadores para mais de 50 capítulos, cobrindo:
- **Astrofísica Relativística**: TOV, Geodésicas de Kerr, Teukolsky e BBN.
- **Microeletrônica**: Transporte quântico DG-DD em FinFETs de 3nm.
- **Biofísica**: Hodgkin-Huxley, Dobramento de Proteínas e Hemodinâmica.
- **Econometria**: Heston, Black-Scholes Fracionário e VECM.
- **Dinâmica de Fluidos**: Navier-Stokes, Esquema de Roe, RANS e Lattice Boltzmann.

## 🛠️ Como Compilar

### Biblioteca SuperDiFF
```bash
cd superdiff
make
```
Isso gerará a biblioteca estática `libsuperdiff.a` e a compartilhada `libsuperdiff.so`.

### Livro (LaTeX)
```bash
cd livro
pdflatex main.tex
```

## 📜 Licença e Créditos
**Autor**: Luiz Tiago Wilcke
**Ano**: 2026
**Linguagem de Variáveis**: Todos os códigos utilizam nomes de variáveis rigorosamente em **Português**, seguindo o padrão pedagógico do manuscrito.

---
*Este trabalho é parte da expansão enciclopédica do conhecimento científico computacional.*
