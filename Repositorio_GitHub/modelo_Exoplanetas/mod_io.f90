MODULE mod_io
  USE iso_fortran_env, ONLY: real64
  IMPLICIT NONE
  
  TYPE :: dados_fotometricos
    REAL(real64), ALLOCATABLE :: tempo(:)
    REAL(real64), ALLOCATABLE :: fluxo(:)
    REAL(real64), ALLOCATABLE :: erro(:)
    INTEGER :: n_pontos
  END TYPE dados_fotometricos

  TYPE :: dados_rv
    REAL(real64), ALLOCATABLE :: tempo(:)
    REAL(real64), ALLOCATABLE :: v_radial(:)
    REAL(real64), ALLOCATABLE :: erro(:)
    INTEGER :: n_pontos
  END TYPE dados_rv

CONTAINS

  SUBROUTINE ler_curva_luz_fits(arq_nome, dados)
    CHARACTER(LEN=*), INTENT(IN) :: arq_nome
    TYPE(dados_fotometricos), INTENT(OUT) :: dados
    ! Simulacao de leitura FITS (Leitura ASCII para simplicidade pedagógica)
    INTEGER :: unit_f, i, ios
    OPEN(NEWUNIT=unit_f, FILE=arq_nome, STATUS='OLD', IOSTAT=ios)
    IF (ios /= 0) THEN; PRINT *, "Erro ao abrir arquivo de luz: ", arq_nome; RETURN; END IF
    
    ! Contagem de linhas
    dados%n_pontos = 0
    DO
      READ(unit_f, *, IOSTAT=ios); IF (ios /= 0) EXIT
      dados%n_pontos = dados%n_pontos + 1
    END DO
    REWIND(unit_f)
    
    ALLOCATE(dados%tempo(dados%n_pontos))
    ALLOCATE(dados%fluxo(dados%n_pontos))
    ALLOCATE(dados%erro(dados%n_pontos))
    
    DO i = 1, dados%n_pontos
      READ(unit_f, *) dados%tempo(i), dados%fluxo(i), dados%erro(i)
    END DO
    CLOSE(unit_f)
    PRINT *, "Lido: ", dados%n_pontos, " pontos de luz de ", arq_nome
  END SUBROUTINE ler_curva_luz_fits

  SUBROUTINE log_execucao(mensagem)
    CHARACTER(LEN=*), INTENT(IN) :: mensagem
    INTEGER :: unit_l; CHARACTER(LEN=8) :: data_at; CHARACTER(LEN=10) :: hora_at
    CALL DATE_AND_TIME(DATE=data_at, TIME=hora_at)
    OPEN(NEWUNIT=unit_l, FILE="execucao_exoplanet.log", STATUS='UNKNOWN', POSITION='APPEND')
    WRITE(unit_l, "(A8, '-', A10, ': ', A)") data_at, hora_at, mensagem
    CLOSE(unit_l)
  END SUBROUTINE log_execucao

END MODULE mod_io
