!==================================================================================================
! MODULO: MOD_DATA_COMPRESSION (HUFFMAN CODING E HASHING)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Implementacao de algoritmos de logica de programacao para gestao de dados
! cientificos. Inclui codificacao de entropia (Huffman) e busca rapida (Hashing).
!
! REFERENCIA: Numerical Recipes in Fortran 90 / 20. Huffman Coding.
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE mod_data_compression
  USE iso_fortran_env
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: huffman_encode, hash_search

  INTEGER, PARAMETER :: DP = real64

  TYPE :: huffman_node
     INTEGER :: freq
     INTEGER :: left = -1, right = -1
  END TYPE huffman_node

CONTAINS

  !------------------------------------------------------------------------------------------------
  ! CODIFICACAO DE HUFFMAN (CONSTRUCAO DA ARVORE)
  !------------------------------------------------------------------------------------------------
  SUBROUTINE huffman_encode(n_chars, frequencies, codes)
    INTEGER, INTENT(IN)  :: n_chars, frequencies(n_chars)
    CHARACTER(LEN=*), INTENT(OUT) :: codes(n_chars)
    
    TYPE(huffman_node) :: tree(2*n_chars - 1)
    INTEGER :: i, j, n_nodes
    INTEGER :: min1, min2, idx1, idx2

    n_nodes = n_chars
    DO i = 1, n_chars
       tree(i)%freq = frequencies(i)
    END DO

    ! Construcao da arvore gulosa (Greedy)
    DO i = 1, n_chars - 1
       min1 = HUGE(1) ; min2 = HUGE(1)
       idx1 = -1 ; idx2 = -1
       
       DO j = 1, n_nodes
          IF (tree(j)%freq > 0 .AND. tree(j)%freq < min1) THEN
             min2 = min1 ; idx2 = idx1
             min1 = tree(j)%freq ; idx1 = j
          ELSE IF (tree(j)%freq > 0 .AND. tree(j)%freq < min2) THEN
             min2 = tree(j)%freq ; idx2 = j
          END IF
       END DO
       
       ! Novo no pai
       n_nodes = n_nodes + 1
       tree(n_nodes)%freq = min1 + min2
       tree(n_nodes)%left = idx1
       tree(n_nodes)%right = idx2
       ! Desativa os nos filhos para a proxima iteracao
       tree(idx1)%freq = -1 ; tree(idx2)%freq = -1
    END DO

    ! (Logica de geracao de strings binarias omitida para brevidade conceitual)
    codes = "BIN_CODE"
  END SUBROUTINE huffman_encode

  !------------------------------------------------------------------------------------------------
  ! BUSCA POR HASHING (TABELA DE ESPALHAMENTO SIMPLES)
  !------------------------------------------------------------------------------------------------
  FUNCTION hash_search(key, n_size) RESULT(idx)
    CHARACTER(LEN=*), INTENT(IN) :: key
    INTEGER, INTENT(IN) :: n_size
    INTEGER :: idx, h, i
    
    ! Funcao de hash polinomial simples
    h = 0
    DO i = 1, LEN(key)
       h = MOD(h * 31 + IACHAR(key(i:i)), n_size)
    END DO
    idx = h + 1
  END FUNCTION hash_search

END MODULE mod_data_compression
