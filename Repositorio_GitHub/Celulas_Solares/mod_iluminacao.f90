MODULE mod_iluminacao
  USE iso_fortran_env
  USE mod_parametros_si
  USE mod_malha
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: n_lambda = 30
  REAL(real64) :: lambdas(n_lambda), irradiancia(n_lambda) ! Espectro AM1.5G
  
CONTAINS

  ! Calcula Geracao Otica G(x) via Beer-Lambert integrando o espectro
  SUBROUTINE calcular_geracao_am15(profundidade)
     REAL(real64), INTENT(OUT) :: profundidade(:)
     REAL(real64) :: alfa, fluxo_fotonico
     INTEGER :: i, l
     
     ! Simulação de espectro discretizado (300nm a 1200nm)
     ! Fluxo = Irradiancia / (E_foton)
     profundidade = 0.0d0
     DO l = 1, n_lambda
        alfa = 1.0d4 * EXP(- (lambdas(l)-400.0d0)/200.0d0 ) ! Simplificacao alpha(lambda)
        fluxo_fotonico = irradiancia(l) / (h_planck * c_luz / (lambdas(l)*1.0d-7))
        DO i = 1, n_pontos
           profundidade(i) = profundidade(i) + fluxo_fotonico * alfa * EXP(-alfa * x(i))
        END DO
     END DO
  END SUBROUTINE calcular_geracao_am15

  SUBROUTINE carregar_espectro_padrao()
    INTEGER :: i
    DO i = 1, n_lambda
       lambdas(i) = 300.0d0 + (i-1)*30.0d0 ! nm
       irradiancia(i) = 1.0e17 * EXP(-(lambdas(i)-550.0d0)**2 / (2.0d0*200.0d0**2))
    END DO
  END SUBROUTINE carregar_espectro_padrao

END MODULE mod_iluminacao
