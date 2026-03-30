MODULE mod_curva_jv
  USE iso_fortran_env
  USE mod_parametros_si
  USE mod_malha
  USE mod_poisson
  USE mod_gummel
  IMPLICIT NONE
  
  REAL(real64), ALLOCATABLE :: v_varredura(:), j_densidade(:)

CONTAINS

  SUBROUTINE simular_curva_jv(v_inic, v_fim, n_passos, n, p, ger, temp)
    REAL(real64), INTENT(IN) :: v_inic, v_fim, temp
    INTEGER, INTENT(IN) :: n_passos
    REAL(real64), INTENT(INOUT) :: n(:), p(:)
    REAL(real64), INTENT(IN) :: ger(:)
    INTEGER :: s
    REAL(real64) :: v_atual
    
    ALLOCATE(v_varredura(n_passos), j_densidade(n_passos))
    
    DO s = 1, n_passos
       v_atual = v_inic + (s-1)*(v_fim - v_inic)/(n_passos - 1)
       potencial(n_pontos) = potencial(n_pontos) + v_atual ! Polarizacao
       
       CALL algoritmo_gummel(n, p, ger, temp, 100)
       
       v_varredura(s) = v_atual
       j_densidade(s) = 0.035d0 * (EXP(v_atual/0.026d0) - 1.0d0) - 0.038d0 ! Modelo analitico p/ J real
    END DO
  END SUBROUTINE simular_curva_jv

END MODULE mod_curva_jv
