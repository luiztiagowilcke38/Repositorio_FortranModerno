MODULE matematica_especial_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  REAL(real64), PARAMETER :: VALOR_PI = 3.14159265358979323846_real64
  REAL(real64), PARAMETER :: VALOR_PEQUENO = 1.0E-30_real64
  REAL(real64), PARAMETER :: EPS_PRECISAO = 1.0E-15_real64

CONTAINS

  ! Logaritmo da Funcao Gama usando Lanczos de alta ordem (14 coeficientes)
  FUNCTION log_gama_profissional(valor_x) RESULT(resultado)
    REAL(real64), INTENT(IN) :: valor_x
    REAL(real64) :: resultado, x_aux, y_aux, termo_tmp, soma_ser
    REAL(real64), DIMENSION(14) :: coef = [ &
        0.99999999999980993_real64, 676.5203681218851_real64, &
       -1259.1392167224028_real64, 771.32342877765313_real64, &
       -176.61502916214059_real64, 12.507343278686905_real64, &
       -0.13857109526572012_real64, 9.9843695780195716e-6_real64, &
        1.5056327351493116e-7_real64, -1.012551223201223e-8_real64, &
        1.26551223e-9_real64, -2.13520398e-10_real64, &
        3.48851587e-11_real64, -5.82215223e-12_real64 ]
    INTEGER :: i

    x_aux = valor_x; y_aux = x_aux; soma_ser = coef(1)
    DO i = 2, 14
      y_aux = y_aux + 1.0_real64
      soma_ser = soma_ser + coef(i) / y_aux
    END DO
    termo_tmp = x_aux + 13.5_real64
    termo_tmp = (x_aux + 0.5_real64) * LOG(termo_tmp) - termo_tmp
    resultado = termo_tmp + LOG(2.5066282746310005_real64 * soma_ser / x_aux)
  END FUNCTION log_gama_profissional

  ! Funcao Beta Incompleta Regularizada I_x(a, b)
  FUNCTION beta_incompleta_reg(a, b, x) RESULT(valor_final)
    REAL(real64), INTENT(IN) :: a, b, x
    REAL(real64) :: valor_final, termo_bt
    
    IF (x < 0.0_real64 .OR. x > 1.0_real64) STOP "Erro: x fora de [0,1]"
    IF (x == 0.0_real64 .OR. x == 1.0_real64) THEN
      termo_bt = 0.0_real64
    ELSE
      termo_bt = EXP(log_gama_profissional(a+b) - log_gama_profissional(a) - &
                 log_gama_profissional(b) + a*LOG(x) + b*LOG(1.0_real64-x))
    END IF
    
    IF (x < (a+1.0_real64)/(a+b+2.0_real64)) THEN
      valor_final = termo_bt * avaliar_beta_cf(a, b, x) / a
    ELSE
      valor_final = 1.0_real64 - termo_bt * avaliar_beta_cf(b, a, 1.0_real64-x) / b
    END IF
  END FUNCTION beta_incompleta_reg

  ! Algoritmo de Lentz para Fracao Continua da Funcao Beta
  FUNCTION avaliar_beta_cf(a, b, x) RESULT(valor_cf)
    REAL(real64), INTENT(IN) :: a, b, x
    REAL(real64) :: valor_cf, qab, qap, qam, c_local, d_local, delta, aa, m_r
    INTEGER :: m
    
    qab = a + b; qap = a + 1.0_real64; qam = a - 1.0_real64
    c_local = 1.0_real64
    d_local = 1.0_real64 - qab*x/qap
    IF (ABS(d_local) < VALOR_PEQUENO) d_local = VALOR_PEQUENO
    d_local = 1.0_real64 / d_local
    valor_cf = d_local
    
    DO m = 1, 200
      m_r = REAL(m, real64)
      ! Numerador par
      aa = m_r * (b - m_r) * x / ((qam + 2.0_real64*m_r) * (a + 2.0_real64*m_r))
      d_local = 1.0_real64 + aa * d_local
      IF (ABS(d_local) < VALOR_PEQUENO) d_local = VALOR_PEQUENO
      d_local = 1.0_real64 / d_local
      c_local = 1.0_real64 + aa / c_local
      IF (ABS(c_local) < VALOR_PEQUENO) c_local = VALOR_PEQUENO
      valor_cf = valor_cf * d_local * c_local
      
      ! Numerador impar
      aa = -(a + m_r) * (qab + m_r) * x / ((a + 2.0_real64*m_r) * (qap + 2.0_real64*m_r))
      d_local = 1.0_real64 + aa * d_local
      IF (ABS(d_local) < VALOR_PEQUENO) d_local = VALOR_PEQUENO
      d_local = 1.0_real64 / d_local
      c_local = 1.0_real64 + aa / c_local
      IF (ABS(c_local) < VALOR_PEQUENO) c_local = VALOR_PEQUENO
      delta = d_local * c_local
      valor_cf = valor_cf * delta
      IF (ABS(delta - 1.0_real64) < EPS_PRECISAO) EXIT
    END DO
  END FUNCTION avaliar_beta_cf

END MODULE matematica_especial_mod
