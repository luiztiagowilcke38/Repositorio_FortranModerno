PROGRAM celula_solar_si
  USE mod_parametros_si
  USE mod_malha
  USE mod_poisson
  USE mod_transporte
  USE mod_iluminacao
  USE mod_gummel
  USE mod_curva_jv
  IMPLICIT NONE
  
  REAL(real64) :: n(n_pontos), p(n_pontos), ger_otica(n_pontos)
  REAL(real64) :: temp = 300.0d0
  INTEGER :: i
  
  PRINT *, "========================================================="
  PRINT *, "  SIMULADOR DE CELULA SOLAR DE SILICIO (ENCICLOPEDIA)    "
  PRINT *, "========================================================="
  
  ! 1. Inicializacao
  CALL construir_malha(espessura=0.03d0, x_juncao=0.0001d0) ! 300um, juncao 1um
  CALL definir_dopagem_estratificada(x_juncao=0.0001d0, Na_base=1.0d16, Nd_emissor=1.0d18)
  
  ! 2. Potencial de Equilíbrio (Equacao de Poisson no escuro)
  ! n = ni exp((phi-ref)/vt), p = ni exp((ref-phi)/vt)
  ! ... (Omitido p/ simplicidade no esqueleto, mas implementado nos modulos) ...
  
  ! 3. Iluminacao AM1.5G
  CALL carregar_espectro_padrao()
  CALL calcular_geracao_am15(ger_otica)
  
  ! 4. Varredura J-V
  CALL simular_curva_jv(v_inic=-0.1d0, v_fim=0.7d0, n_passos=40, n=n, p=p, ger=ger_otica, temp=temp)
  
  ! 5. Salvar Resultados
  OPEN(unit=10, file='curva_jv.dat')
  DO i = 1, SIZE(v_varredura)
     WRITE(10, '(2F12.6)') v_varredura(i), j_densidade(i)
  END DO
  CLOSE(10)
  
  PRINT *, "Simulacao concluida. Resultados salvos em curva_jv.dat"
  PRINT *, "========================================================="
END PROGRAM celula_solar_si
