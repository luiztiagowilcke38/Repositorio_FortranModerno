MODULE mod_signal
  USE iso_fortran_env, ONLY: real64
  USE mod_constants, ONLY: PI_CONST
  IMPLICIT NONE

CONTAINS

  SUBROUTINE remover_tendencia_spline(tempo, fluxo, n, k_nos, fluxo_detrend)
    REAL(real64), DIMENSION(:), INTENT(IN) :: tempo, fluxo
    INTEGER, INTENT(IN) :: n, k_nos
    REAL(real64), DIMENSION(:), INTENT(OUT) :: fluxo_detrend
    ! Implementacao simplificada de remocao de tendencia via spline de baixa ordem
    fluxo_detrend = fluxo - (SUM(fluxo)/REAL(n, real64))
  END SUBROUTINE remover_tendencia_spline

  SUBROUTINE sigma_clipping_iterativo(dados, sigma_lim, mascara)
    REAL(real64), DIMENSION(:), INTENT(IN) :: dados
    REAL(real64), INTENT(IN) :: sigma_lim
    LOGICAL, DIMENSION(:), INTENT(OUT) :: mascara
    REAL(real64) :: media, desvio
    INTEGER :: iter, n
    n = SIZE(dados); mascara = .TRUE.
    DO iter = 1, 5
      media = SUM(dados, MASK=mascara) / REAL(COUNT(mascara), real64)
      desvio = SQRT(SUM((dados - media)**2, MASK=mascara) / REAL(COUNT(mascara), real64))
      mascara = mascara .AND. (ABS(dados - media) <= sigma_lim * desvio)
    END DO
  END SUBROUTINE sigma_clipping_iterativo

  SUBROUTINE lomb_scargle_gls(tempo, fluxo, freq_grid, potencia)
    REAL(real64), DIMENSION(:), INTENT(IN) :: tempo, fluxo
    REAL(real64), DIMENSION(:), INTENT(IN) :: freq_grid
    REAL(real64), DIMENSION(:), INTENT(OUT) :: potencia
    INTEGER :: i, j, n, m; REAL(real64) :: arg, s_ang, c_ang, yy, yc, ys, cc, ss, cs, tau
    n = SIZE(tempo); m = SIZE(freq_grid)
    ! Algoritmo GLS (Generalized Lomb-Scargle) simplificado para demonstracao
    DO j = 1, m
      yc=0; ys=0; cc=0; ss=0; cs=0
      DO i = 1, n
        arg = 2.0_real64 * PI_CONST * freq_grid(j) * tempo(i)
        c_ang = COS(arg); s_ang = SIN(arg)
        yc = yc + fluxo(i)*c_ang; ys = ys + fluxo(i)*s_ang
        cc = cc + c_ang**2; ss = ss + s_ang**2; cs = cs + c_ang*s_ang
      END DO
      potencia(j) = (yc**2 * ss + ys**2 * cc - 2.0*yc*ys*cs) / (cc*ss - cs**2)
    END DO
  END SUBROUTINE lomb_scargle_gls

  SUBROUTINE test_modulo_signal()
    REAL(real64), DIMENSION(10) :: t, f
    LOGICAL, DIMENSION(10) :: m
    INTEGER :: i
    DO i = 1, 10
      t(i) = REAL(i, real64)
      f(i) = 1.0_real64
    END DO
    f(5) = 100.0_real64 ! Outlier
    CALL sigma_clipping_iterativo(f, 3.0_real64, m)
    IF (.NOT. m(5)) THEN
      PRINT *, "[OK] Modulo de Sinal (Sigma-Clipping) validado."
    ELSE
      PRINT *, "[ERRO] Outlier nao detectado pelo sigma-clipping."
    END IF
  END SUBROUTINE test_modulo_signal

END MODULE mod_signal
