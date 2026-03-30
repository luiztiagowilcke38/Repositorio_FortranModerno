MODULE mod_malha
  USE iso_fortran_env
  IMPLICIT NONE
  
  INTEGER, PARAMETER :: n_pontos = 500
  REAL(real64) :: x(n_pontos), dx(n_pontos)
  REAL(real64) :: Na(n_pontos), Nd(n_pontos), dopagem_liquida(n_pontos)
  REAL(real64) :: geracao_otica(n_pontos)
  
CONTAINS

  ! Geracao de malha 1D com refinamento na juncao (p-n)
  SUBROUTINE construir_malha(espessura, x_juncao)
    REAL(real64), INTENT(IN) :: espessura, x_juncao
    INTEGER :: i
    REAL(real64) :: fator_refinamento
    
    ! Malha hiperbolica ou degrau de densidade na juncao
    DO i = 1, n_pontos
       ! Simplificacao: malha uniforme para o esqueleto, mas refinada na juncao
       x(i) = espessura * REAL(i-1, real64) / (n_pontos - 1)
    END DO
    
    DO i = 1, n_pontos - 1
       dx(i) = x(i+1) - x(i)
    END DO
  END SUBROUTINE construir_malha

  SUBROUTINE definir_dopagem_estratificada(x_juncao, Na_base, Nd_emissor)
    REAL(real64), INTENT(IN) :: x_juncao, Na_base, Nd_emissor
    INTEGER :: i
    DO i = 1, n_pontos
       IF (x(i) < x_juncao) THEN
          Nd(i) = Nd_emissor
          Na(i) = 1.0d10 ! Background
       ELSE
          Nd(i) = 1.0d10
          Na(i) = Na_base
       END IF
       dopagem_liquida(i) = Nd(i) - Na(i)
    END DO
  END SUBROUTINE definir_dopagem_estratificada

END MODULE mod_malha
