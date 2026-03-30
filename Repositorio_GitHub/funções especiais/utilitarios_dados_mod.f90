MODULE utilitarios_dados_mod
  USE iso_fortran_env, ONLY: int32, int64
  IMPLICIT NONE

  INTEGER(int32), DIMENSION(0:255), SAVE :: tab_crc32
  LOGICAL, SAVE :: iniciada = .FALSE.

CONTAINS

  SUBROUTINE iniciar_tabela_crc()
    INTEGER(int32) :: i, j, c
    DO i = 0, 255
      c = i
      DO j = 1, 8
        IF (IAND(c, 1_int32) /= 0) THEN
          c = IEOR(Z'EDB88320', ISHFT(c, -1))
        ELSE
          c = ISHFT(c, -1)
        END IF
      END DO
      tab_crc32(i) = c
    END DO
    iniciada = .TRUE.
  END SUBROUTINE iniciar_tabela_crc

  FUNCTION obter_crc32_buffer(dados) RESULT(crc_res)
    INTEGER(1), DIMENSION(:), INTENT(IN) :: dados
    INTEGER(int32) :: crc_res, k
    IF (.NOT. iniciada) CALL iniciar_tabela_crc()
    crc_res = NOT(0_int32)
    DO k = 1, SIZE(dados)
      crc_res = IEOR(tab_crc32(IAND(IEOR(crc_res, INT(dados(k), int32)), 255_int32)), &
                ISHFT(crc_res, -8))
    END DO
    crc_res = NOT(crc_res)
  END FUNCTION obter_crc32_buffer

  FUNCTION para_gray_v(valor) RESULT(g)
    INTEGER(int64), INTENT(IN) :: valor; INTEGER(int64) :: g
    g = IEOR(valor, ISHFT(valor, -1))
  END FUNCTION para_gray_v

  FUNCTION de_gray_v(gray) RESULT(v)
    INTEGER(int64), INTENT(IN) :: gray; INTEGER(int64) :: v, t
    t = gray; v = t
    DO WHILE (t > 0)
      t = ISHFT(t, -1); v = IEOR(v, t)
    END DO
  END FUNCTION de_gray_v

END MODULE utilitarios_dados_mod
