MODULE mod_nested
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: N_PONTOS_VIVOS = 500
  REAL(real64), PARAMETER :: CONV_NESTED = 0.1_real64

CONTAINS

  SUBROUTINE executar_nesting(log_evidencia)
    REAL(real64), INTENT(OUT) :: log_evidencia
    REAL(real64), DIMENSION(N_PONTOS_VIVOS) :: log_l_vivos
    REAL(real64) :: log_w, log_z, log_z_novo, l_min
    INTEGER :: iter, i_min
    
    log_w = -LOG(REAL(N_PONTOS_VIVOS, real64))
    log_z = -HUGE(1.0_real64)
    
    ! Inicializacao de pontos vivos (Valores iniciais simulados)
    log_l_vivos = -100.0_real64; iter = 0
    
    DO WHILE (iter < 5000)
      iter = iter + 1; i_min = MINLOC(log_l_vivos, 1); l_min = log_l_vivos(i_min)
      
      ! Atualizacao da evidencia pelo metodo dos trapezios no espaco logaritmico
      log_z_novo = log_soma_exp(log_z, log_w + l_min)
      
      ! Criterio de convergencia: se a mudanca na evidencia for pequena
      IF (ABS(log_z_novo - log_z) < CONV_NESTED) EXIT
      log_z = log_z_novo
      
      ! Substituicao do ponto de menor verossimilhanca por um novo ponto do Prior
      log_l_vivos(i_min) = l_min + 1.0_real64 ! Simulando subida na "montanha"
      log_w = log_w - 1.0_real64 / REAL(N_PONTOS_VIVOS, real64)
    END DO
    log_evidencia = log_z
    PRINT *, "Nested Sampling Finalizado. Log Z =", log_evidencia
  END SUBROUTINE executar_nesting

  FUNCTION log_soma_exp(la, lb) RESULT(res)
    REAL(real64), INTENT(IN) :: la, lb; REAL(real64) :: res
    IF (la > lb) THEN; res = la + LOG(1.0_real64 + EXP(lb - la))
    ELSE; res = lb + LOG(1.0_real64 + EXP(la - lb)); END IF
  END FUNCTION log_soma_exp

END MODULE mod_nested
