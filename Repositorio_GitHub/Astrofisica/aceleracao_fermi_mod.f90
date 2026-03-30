MODULE aceleracao_fermi_mod
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  REAL(real64), PARAMETER :: VEL_CHOQUE = 5000.0E3_real64 ! 5000 km/s
  REAL(real64), PARAMETER :: VEL_LUZ = 2.9979E8_real64
  INTEGER, PARAMETER :: N_PARTICULAS = 10000

CONTAINS

  SUBROUTINE simular_espectro_fermi(energia_inicial, p_escape)
    REAL(real64), INTENT(IN) :: energia_inicial, p_escape
    REAL(real64), DIMENSION(N_PARTICULAS) :: energias_finais
    REAL(real64) :: r_rand, e_temp, ganho_medio
    INTEGER :: i, n_ciclos
    
    ganho_medio = 1.0_real64 + (4.0_real64/3.0_real64) * (VEL_CHOQUE / VEL_LUZ)
    CALL RANDOM_SEED()
    
    DO i = 1, N_PARTICULAS
      e_temp = energia_inicial; n_ciclos = 0
      DO
        CALL RANDOM_NUMBER(r_rand)
        IF (r_rand < p_escape) EXIT ! Particula escapou do choque
        
        ! Ganho de energia no ciclo de Fermi
        e_temp = e_temp * ganho_medio; n_ciclos = n_ciclos + 1
        IF (n_ciclos > 1000) EXIT ! Limite de aceleracao
      END DO
      energias_finais(i) = e_temp
    END DO
    
    CALL analisar_lei_potencia(energias_finais)
  END SUBROUTINE simular_espectro_fermi

  SUBROUTINE analisar_lei_potencia(e_vetor)
    REAL(real64), DIMENSION(:), INTENT(IN) :: e_vetor
    REAL(real64) :: e_min, e_max, de, contagem
    INTEGER :: j, n_bins = 20; REAL(real64) :: bin_low, bin_high
    
    e_min = MINVAL(e_vetor); e_max = MAXVAL(e_vetor)
    PRINT *, "Analise de Espectro (Log-Log):"
    DO j = 1, n_bins
      bin_low = e_min * (e_max/e_min)**(REAL(j-1, real64)/n_bins)
      bin_high = e_min * (e_max/e_min)**(REAL(j, real64)/n_bins)
      contagem = COUNT(e_vetor >= bin_low .AND. e_vetor < bin_high)
      IF (contagem > 0) THEN
        PRINT "(F12.3, F10.1)", LOG10(bin_low), LOG10(contagem)
      END IF
    END DO
  END SUBROUTINE analisar_lei_potencia

END MODULE aceleracao_fermi_mod
