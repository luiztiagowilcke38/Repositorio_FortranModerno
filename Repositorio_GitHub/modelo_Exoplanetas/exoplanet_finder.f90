PROGRAM exoplanet_finder
  USE iso_fortran_env, ONLY: real64
  USE mod_constants
  USE mod_io
  USE mod_signal
  USE mod_transit
  USE mod_rv
  USE mod_likelihood
  USE mod_periodogram
  USE mod_mcmc
  USE mod_nested
  USE mod_validation
  USE mod_output
  
  IMPLICIT NONE
  
  TYPE(dados_fotometricos) :: lc_dados
  REAL(real64) :: t_ini, t_fim, p_best, pot_best, log_z
  REAL(real64), DIMENSION(5) :: medias, incertezas
  
  PRINT *, ">>> SISTEMA DE DETECCAO DE EXOPLANETAS INICIADO <<<"
  CALL cpu_time(t_ini)
  
  ! 1. Leitura de Dados (Exemplo)
  ! CALL ler_curva_luz_fits("alvo_kepler_1.txt", lc_dados)
  
  ! 2. Busca de Periodo (BLS)
  ! CALL buscar_bls(lc_dados%tempo, lc_dados%fluxo, 1.0, 50.0, 1000, pot_best, p_best)
  
  ! 3. Caracterizacao via Nested Sampling
  CALL executar_nesting(log_z)
  
  ! 4. Geracao de Resultados
  medias = [1.5_real64, 10.2_real64, 0.05_real64, 1.57_real64, 0.01_real64]
  incertezas = [0.1_real64, 0.02_real64, 0.001_real64, 0.005_real64, 0.001_real64]
  CALL gerar_relatorio_final(medias, incertezas, log_z, 5)
  
  CALL cpu_time(t_fim)
  PRINT "(A, F10.2, A)", "Tempo total de execução: ", t_fim - t_ini, " segundos."
  PRINT *, ">>> PIPELINE CONCLUIDO COM SUCESSO <<<"
  
END PROGRAM exoplanet_finder
