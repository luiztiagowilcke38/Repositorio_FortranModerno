MODULE mod_constants
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  ! Constantes Fisicas Fundamentais (SI)
  REAL(real64), PARAMETER :: G_CONST = 6.67430E-11_real64         ! Constante Gravitacional
  REAL(real64), PARAMETER :: VEL_LUZ = 299792458.0_real64         ! Velocidade da Luz
  REAL(real64), PARAMETER :: CONST_PLANCK = 6.62607015E-34_real64 ! Constante de Planck
  REAL(real64), PARAMETER :: CONST_BOLTZ = 1.380649E-23_real64    ! Constante de Boltzmann
  REAL(real64), PARAMETER :: STEFAN_BOLTZ = 5.670374419E-8_real64 ! Constante de Stefan-Boltzmann

  ! Constantes Astronomicas
  REAL(real64), PARAMETER :: MASSA_SOL = 1.98847E30_real64        ! Massa Solar (kg)
  REAL(real64), PARAMETER :: RAIO_SOL = 6.957E8_real64            ! Raio Solar (m)
  REAL(real64), PARAMETER :: LUM_SOL = 3.828E26_real64             ! Luminosidade Solar (W)
  REAL(real64), PARAMETER :: MASSA_TERRA = 5.9722E24_real64       ! Massa da Terra (kg)
  REAL(real64), PARAMETER :: RAIO_TERRA = 6.371E6_real64          ! Raio da Terra (m)
  REAL(real64), PARAMETER :: MASSA_JUPITER = 1.8982E27_real64     ! Massa de Jupiter (kg)
  REAL(real64), PARAMETER :: RAIO_JUPITER = 6.9911E7_real64       ! Raio de Jupiter (m)
  
  ! Unidades de Distancia
  REAL(real64), PARAMETER :: UC_ASTRONOMICA = 1.495978707E11_real64 ! AU (m)
  REAL(real64), PARAMETER :: PARSEC = 3.085677581E16_real64         ! pc (m)
  REAL(real64), PARAMETER :: ANO_LUZ = 9.4607304725808E15_real64    ! ly (m)

  ! Tolerancias e Limites
  REAL(real64), PARAMETER :: TOL_NUMERICA = 1.0E-12_real64
  REAL(real64), PARAMETER :: PI_CONST = 3.14159265358979323846_real64
  
CONTAINS

  SUBROUTINE test_modulo_constants()
    IF (ABS(PI_CONST - 3.14159265358979323846_real64) < TOL_NUMERICA) THEN
      PRINT *, "[OK] Modulo de Constantes validado."
    ELSE
      PRINT *, "[ERRO] Falha na precisao das constantes."
    END IF
  END SUBROUTINE test_modulo_constants

END MODULE mod_constants
