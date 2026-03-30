MODULE ordenacao_avancada_mod
  USE iso_fortran_env, ONLY: real64, int32
  IMPLICIT NONE

CONTAINS

  SUBROUTINE heapsort_vetor(vetor_dados)
    REAL(real64), DIMENSION(:), INTENT(INOUT) :: vetor_dados
    INTEGER(int32) :: n_tam, i_p
    REAL(real64) :: aux
    
    n_tam = SIZE(vetor_dados)
    IF (n_tam < 2) RETURN
    
    DO i_p = n_tam/2, 1, -1
      CALL organizar_heap(vetor_dados, i_p, n_tam)
    END DO
    
    DO i_p = n_tam, 2, -1
      aux = vetor_dados(1); vetor_dados(1) = vetor_dados(i_p); vetor_dados(i_p) = aux
      CALL organizar_heap(vetor_dados, 1, i_p-1)
    END DO
  END SUBROUTINE heapsort_vetor

  SUBROUTINE organizar_heap(v, ind_esq, ind_dir)
    REAL(real64), DIMENSION(:), INTENT(INOUT) :: v
    INTEGER(int32), INTENT(IN) :: ind_esq, ind_dir
    INTEGER(int32) :: i, j
    REAL(real64) :: tmp
    i = ind_esq; j = 2*i; tmp = v(i)
    DO WHILE (j <= ind_dir)
      IF (j < ind_dir) THEN
        IF (v(j) < v(j+1)) j = j + 1
      END IF
      IF (tmp >= v(j)) EXIT
      v(i) = v(j); i = j; j = 2*i
    END DO
    v(i) = tmp
  END SUBROUTINE organizar_heap

  FUNCTION selecionar_kesimo_final(vetor_ref, k_desejado) RESULT(res_k)
    REAL(real64), DIMENSION(:), INTENT(IN) :: vetor_ref
    INTEGER(int32), INTENT(IN) :: k_desejado
    REAL(real64) :: res_k, pivo, troca
    REAL(real64), ALLOCATABLE :: copia(:)
    INTEGER(int32) :: n, esq, dir, i, j
    
    n = SIZE(vetor_ref); ALLOCATE(copia(n)); copia = vetor_ref
    esq = 1; dir = n
    DO WHILE (esq < dir)
      pivo = copia(k_desejado); i = esq; j = dir
      DO
        DO WHILE (copia(i) < pivo); i = i + 1; END DO
        DO WHILE (copia(j) > pivo); j = j - 1; END DO
        IF (i <= j) THEN
          troca = copia(i); copia(i) = copia(j); copia(j) = troca
          i = i + 1; j = j - 1
        END IF
        IF (i > j) EXIT
      END DO
      IF (j < k_desejado) esq = i
      IF (i > k_desejado) dir = j
    END DO
    res_k = copia(k_desejado); DEALLOCATE(copia)
  END FUNCTION selecionar_kesimo_final

END MODULE ordenacao_avancada_mod
