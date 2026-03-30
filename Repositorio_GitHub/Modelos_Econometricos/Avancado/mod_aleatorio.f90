MODULE mod_aleatorio
  USE iso_fortran_env
  IMPLICIT NONE
  INTEGER(int64) :: estado_semente = 123456789_int64

CONTAINS

  REAL(real64) FUNCTION aleatorio_uniforme()
    estado_semente = IEOR(estado_semente, ISHL(estado_semente, 13))
    estado_semente = IEOR(estado_semente, ISHR(estado_semente, 7))
    estado_semente = IEOR(estado_semente, ISHL(estado_semente, 17))
    aleatorio_uniforme = ABS(REAL(estado_semente, real64) / REAL(HUGE(estado_semente), real64))
  END FUNCTION aleatorio_uniforme

  REAL(real64) FUNCTION aleatorio_normal()
    REAL(real64) :: u1, u2
    u1 = aleatorio_uniforme()
    u2 = aleatorio_uniforme()
    aleatorio_normal = SQRT(-2.0d0 * LOG(u1)) * COS(2.0d0 * 3.14159d0 * u2)
  END FUNCTION aleatorio_normal

END MODULE mod_aleatorio
