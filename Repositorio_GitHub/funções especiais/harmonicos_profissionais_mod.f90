MODULE harmonicos_profissionais_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  REAL(real64), PARAMETER :: VALOR_PI = 3.14159265358979323846_real64

CONTAINS

  ! Calcula o Harmônico Esférico Y_lm renormalizado
  SUBROUTINE avaliar_harmonico_esferico(grau_l, ordem_m, coord_teta, coord_fi, &
                                        parte_real, parte_imag)
    INTEGER, INTENT(IN) :: grau_l, ordem_m
    REAL(real64), INTENT(IN) :: coord_teta, coord_fi
    REAL(real64), INTENT(OUT) :: parte_real, parte_imag
    REAL(real64) :: normalizacao, leg_assoc, cosseno_teta, m_fi_local
    REAL(real64) :: termo_log_fact
    
    cosseno_teta = COS(coord_teta)
    leg_assoc = calcular_legendre_mestre(grau_l, ABS(ordem_m), cosseno_teta)
    
    ! Calculo da normalizacao via log-fatorial para evitar overflow em l > 100
    termo_log_fact = log_fatorial_seguro(grau_l - ABS(ordem_m)) - &
                     log_fatorial_seguro(grau_l + ABS(ordem_m))
    
    normalizacao = SQRT(((2*grau_l + 1) / (4.0_real64 * VALOR_PI)) * EXP(termo_log_fact))
    
    m_fi_local = REAL(ordem_m, real64) * coord_fi
    parte_real = normalizacao * leg_assoc * COS(m_fi_local)
    parte_imag = normalizacao * leg_assoc * SIN(m_fi_local)
    
    ! Ajuste de paridade de Condon-Shortley
    IF (ordem_m < 0 .AND. MOD(ABS(ordem_m), 2) /= 0) THEN
      parte_real = -parte_real; parte_imag = -parte_imag
    END IF
  END SUBROUTINE avaliar_harmonico_esferico

  ! Polinomio de Legendre P_lm via algoritmo de recorrencia estavel
  FUNCTION calcular_legendre_mestre(l, m, x_val) RESULT(polm)
    INTEGER, INTENT(IN) :: l, m
    REAL(real64), INTENT(IN) :: x_val
    REAL(real64) :: polm, p_m_m, p_m_m_mais1, sen_teta, termo_fatorial
    INTEGER :: k, ll
    
    ! Caso base dependente de (sen theta)^m
    p_m_m = 1.0_real64
    IF (m > 0) THEN
      sen_teta = SQRT((1.0_real64 - x_val) * (1.0_real64 + x_val))
      termo_fatorial = 1.0_real64
      DO k = 1, m
        p_m_m = -p_m_m * termo_fatorial * sen_teta
        termo_fatorial = termo_fatorial + 2.0_real64
      END DO
    END IF
    
    IF (l == m) THEN
      polm = p_m_m
    ELSE
      p_m_m_mais1 = x_val * (2*m + 1) * p_m_m
      IF (l == m + 1) THEN
        polm = p_m_m_mais1
      ELSE
        ! Algoritmo de tres termos para l > m+1
        DO ll = m + 2, l
          polm = (x_val * (2*ll - 1) * p_m_m_mais1 - (ll + m - 1) * p_m_m) / (ll - m)
          p_m_m = p_m_m_mais1
          p_m_m_mais1 = polm
        END DO
      END IF
    END IF
  END FUNCTION calcular_legendre_mestre

  ! Funcao Log-Fatorial para grandes argumentos
  FUNCTION log_fatorial_seguro(valor_n) RESULT(res_log)
    INTEGER, INTENT(IN) :: valor_n
    REAL(real64) :: res_log
    INTEGER :: idx
    res_log = 0.0_real64
    IF (valor_n <= 1) RETURN
    DO idx = 2, valor_n
      res_log = res_log + LOG(REAL(idx, real64))
    END DO
  END FUNCTION log_fatorial_seguro

END MODULE harmonicos_profissionais_mod
