MODULE mod_output
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE

CONTAINS

  SUBROUTINE gerar_relatorio_final(medias, sigmas, log_z, n_par)
    REAL(real64), DIMENSION(:), INTENT(IN) :: medias, sigmas
    REAL(real64), INTENT(IN) :: log_z
    INTEGER, INTENT(IN) :: n_par
    INTEGER :: i
    
    PRINT *, "===================================================="
    PRINT *, "   RELATORIO DE CARACTERIZACAO DE EXOPLANETA"
    PRINT *, "===================================================="
    PRINT "(A, F12.4)", " Evidencia Bayesiana (Log Z): ", log_z
    PRINT *, "----------------------------------------------------"
    PRINT *, " Parâmetro      |  Val. Mediana  |  Incerteza (1-sigma)"
    DO i = 1, n_par
      PRINT "(A15, ' | ', F14.6, ' | ', F16.6)", "Par_"//achar(i), medias(i), sigmas(i)
    END DO
    PRINT *, "===================================================="
  END SUBROUTINE gerar_relatorio_final

  FUNCTION achar(i) RESULT(str)
    INTEGER, INTENT(IN) :: i; CHARACTER(LEN=2) :: str
    WRITE(str, "(I2)") i
  END FUNCTION achar

  SUBROUTINE plotar_curva_ascii(tempo, fluxo)
    REAL(real64), DIMENSION(:), INTENT(IN) :: tempo, fluxo
    INTEGER :: i, j, pos_y
    CHARACTER(LEN=60) :: linha
    PRINT *, "Visualizacao ASCII da Curva de Luz Dobrada:"
    DO i = 1, SIZE(tempo), SIZE(tempo)/20
        linha = " "
        pos_y = INT((1.0_real64 - fluxo(i)) * 500.0_real64) + 1
        IF (pos_y > 0 .AND. pos_y <= 60) linha(pos_y:pos_y) = "*"
        PRINT "(F8.3, '|', A)", tempo(i), linha
    END DO
  END SUBROUTINE plotar_curva_ascii

END MODULE mod_output
