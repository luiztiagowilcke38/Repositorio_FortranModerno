MODULE elipticas_avancadas_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
CONTAINS

  ! Calcula sn, cn e dn simultaneamente para um modulo m
  SUBROUTINE jacobi_mestre_portugues(u_entrada, m_modulo, sn_saida, cn_saida, dn_saida)
    REAL(real64), INTENT(IN) :: u_entrada, m_modulo
    REAL(real64), INTENT(OUT) :: sn_saida, cn_saida, dn_saida
    REAL(real64) :: a_atual, b_atual, c_atual, d_aux, phi_amplitude
    REAL(real64), DIMENSION(20) :: seq_a, seq_c
    INTEGER :: idx, n_convergencia
    
    ! Casos degenerados: m=0 (circular) e m=1 (hiperbolico)
    IF (m_modulo < 1E-16_real64) THEN
      sn_saida = SIN(u_entrada); cn_saida = COS(u_entrada); dn_saida = 1.0_real64; RETURN
    ELSE IF (ABS(m_modulo - 1.0_real64) < 1E-16_real64) THEN
      sn_saida = TANH(u_entrada); cn_saida = 1.0_real64/COSH(u_entrada); dn_saida = cn_saida; RETURN
    END IF

    ! Fase AGM: Transformada de Landen para reduzir o modulo
    a_atual = 1.0_real64
    b_atual = SQRT(1.0_real64 - m_modulo)
    c_atual = SQRT(m_modulo)
    n_convergencia = 0
    
    DO idx = 1, 20
      n_convergencia = idx
      seq_a(n_convergencia) = a_atual; seq_c(n_convergencia) = c_atual
      a_atual = 0.5_real64 * (a_atual + b_atual)
      b_atual = SQRT(seq_a(idx) * b_atual)
      c_atual = 0.5_real64 * (seq_a(idx) - b_atual)
      IF (ABS(c_atual) < 1E-15_real64) EXIT
    END DO
    
    ! Fase de Retorno: Calculo da amplitude phi
    phi_amplitude = a_atual * u_entrada * (2.0_real64**n_convergencia)
    DO idx = n_convergencia, 1, -1
      ! Recursao trigonometrica profunda para reverter a transformada
      phi_amplitude = 0.5_real64 * (phi_amplitude + ASIN(seq_c(idx)/seq_a(idx) * SIN(phi_amplitude)))
    END DO
    
    ! Resultado final das funcoes elipticas
    sn_saida = SIN(phi_amplitude)
    cn_saida = COS(phi_amplitude)
    dn_saida = SQRT(1.0_real64 - m_modulo * sn_saida**2)
    
  END SUBROUTINE jacobi_mestre_portugues

END MODULE elipticas_avancadas_mod
