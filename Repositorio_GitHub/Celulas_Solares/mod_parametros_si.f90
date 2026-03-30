MODULE mod_parametros_si
  USE iso_fortran_env
  IMPLICIT NONE
  
  ! Constantes Fisicas Universais
  REAL(real64), PARAMETER :: q_eletron = 1.602176634d-19    ! C
  REAL(real64), PARAMETER :: k_boltzmann = 1.380649d-23    ! J/K
  REAL(real64), PARAMETER :: eps0 = 8.8541878128d-14      ! F/cm
  REAL(real64), PARAMETER :: h_planck = 6.62607015d-34    ! J.s
  REAL(real64), PARAMETER :: c_luz = 2.99792458d10        ! cm/s
  
  ! Parametros do Silicio (Si)
  REAL(real64), PARAMETER :: eps_si = 11.7d0              ! Constante dieletrica
  REAL(real64), PARAMETER :: massa_efetiva_eletron = 1.08d0 ! m*/m0
  REAL(real64), PARAMETER :: massa_efetiva_lacuna = 0.81d0
  
  ! Coeficientes de Recombinacao Auger (Si)
  REAL(real64), PARAMETER :: C_auger_n = 2.8d-31          ! cm^6/s
  REAL(real64), PARAMETER :: C_auger_p = 9.9d-32          ! cm^6/s
  REAL(real64), PARAMETER :: B_radiativo = 4.73d-15       ! cm^3/s

CONTAINS

  ! Equacao de Varshni para o Gap de Energia Eg(T)
  REAL(real64) FUNCTION calcular_eg(temperatura)
    REAL(real64), INTENT(IN) :: temperatura
    REAL(real64), PARAMETER :: eg0 = 1.17d0, alfa = 4.73d-4, beta = 636.0d0
    calcular_eg = eg0 - (alfa * temperatura**2) / (temperatura + beta)
  END FUNCTION calcular_eg

  ! Concentracao Intrinseca ni(T)
  REAL(real64) FUNCTION calcular_ni(temperatura)
    REAL(real64), INTENT(IN) :: temperatura
    REAL(real64) :: eg, vt
    REAL(real64), PARAMETER :: Nc300 = 2.86d19, Nv300 = 3.10d19
    eg = calcular_eg(temperatura)
    vt = (k_boltzmann * temperatura) / q_eletron
    calcular_ni = SQRT(Nc300 * Nv300) * (temperatura/300.0d0)**1.5d0 * EXP(-eg/(2.0d0 * vt))
  END FUNCTION calcular_ni

  ! Modelo de Caughey-Thomas para Mobilidade dependente de dopagem
  REAL(real64) FUNCTION mobilidade_eletrons(dopagem_total, temperatura)
    REAL(real64), INTENT(IN) :: dopagem_total, temperatura
    REAL(real64), PARAMETER :: mu_min = 92.0d0, mu_max = 1360.0d0, n_ref = 1.3d17, alfa = 0.91d0
    mobilidade_eletrons = mu_min + (mu_max - mu_min) / (1.0d0 + (dopagem_total/n_ref)**alfa)
    ! Correcao simplificada de temperatura: mu ~ T^-2.5
    mobilidade_eletrons = mobilidade_eletrons * (300.0d0/temperatura)**2.5d0
  END FUNCTION mobilidade_eletrons

  REAL(real64) FUNCTION mobilidade_lacunas(dopagem_total, temperatura)
    REAL(real64), INTENT(IN) :: dopagem_total, temperatura
    REAL(real64), PARAMETER :: mu_min = 47.7d0, mu_max = 495.0d0, n_ref = 6.3d16, alfa = 0.76d0
    mobilidade_lacunas = mu_min + (mu_max - mu_min) / (1.0d0 + (dopagem_total/n_ref)**alfa)
    mobilidade_lacunas = mobilidade_lacunas * (300.0d0/temperatura)**2.2d0
  END FUNCTION mobilidade_lacunas

END MODULE mod_parametros_si
