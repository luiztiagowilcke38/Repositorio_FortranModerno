!==================================================================================================
! PROGRAMA: SOLUCIONADOR_SOLAR_SI_P_ULTRA_COMPLEXO (MANUAL CIENTIFICO INTEGRAL)
!--------------------------------------------------------------------------------------------------
! DESCRICAO: Simulador de transporte quantico NEGF para celula solar de atomo unico.
! Implementacao com Bir-Pikus, SCBA Eliashberg e Espectro AM1.5G explicito.
!
! Autor: Luiz Tiago Wilcke, 2026.
!==================================================================================================

MODULE constantes_quanticas
  USE iso_fortran_env
  IMPLICIT NONE
  SAVE
  INTEGER, PARAMETER :: DP = real64
  REAL(DP), PARAMETER :: PI = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: H_BARRA = 0.658211956_DP ! meV * fs
  REAL(DP), PARAMETER :: K_B = 0.086173332_DP ! meV/K
  REAL(DP), PARAMETER :: T_REF = 300.0_DP
  COMPLEX(DP), PARAMETER :: CI = (0.0_DP, 1.0_DP)

CONTAINS
  FUNCTION identity_matrix_12() RESULT(mat)
    COMPLEX(DP) :: mat(12,12)
    INTEGER :: k_idx
    mat = 0.0_DP
    DO k_idx = 1, 12 ; mat(k_idx,k_idx) = 1.0_DP ; END DO
  END FUNCTION identity_matrix_12
END MODULE constantes_quanticas

MODULE dataset_am15g_explicit
  USE constantes_quanticas
  IMPLICIT NONE
  SAVE
  REAL(DP), PARAMETER :: SOLAR_WV(400) = [ (280.0 + REAL(i-1, DP)*2.0, i=1,400) ]
  REAL(DP), PARAMETER :: SOLAR_IRR(400) = [ &
    0.000, 0.000, 0.000, 0.001, 0.002, 0.004, 0.008, 0.015, 0.025, 0.040, &
    0.060, 0.085, 0.120, 0.165, 0.220, 0.285, 0.360, 0.445, 0.540, 0.640, &
    0.745, 0.850, 0.955, 1.055, 1.150, 1.240, 1.320, 1.390, 1.450, 1.500, &
    1.545, 1.585, 1.620, 1.650, 1.675, 1.695, 1.710, 1.725, 1.735, 1.745, &
    1.755, 1.765, 1.775, 1.785, 1.795, 1.805, 1.815, 1.825, 1.835, 1.845, &
    1.855, 1.865, 1.875, 1.885, 1.895, 1.905, 1.915, 1.925, 1.935, 1.945, &
    1.955, 1.965, 1.975, 1.985, 1.995, 2.010, 2.025, 2.040, 2.055, 2.070, &
    2.085, 2.100, 2.115, 2.130, 2.145, 2.160, 2.175, 2.200, 2.225, 2.250, &
    2.275, 2.300, 2.325, 2.350, 2.375, 2.400, 2.425, 2.415, 2.405, 2.395, &
    2.385, 2.375, 2.365, 2.355, 2.345, 2.335, 2.325, 2.315, 2.305, 2.295, &
    2.285, 2.275, 2.265, 2.255, 2.245, 2.235, 2.225, 2.215, 2.205, 2.195, &
    2.185, 2.175, 2.165, 2.155, 2.145, 2.135, 2.125, 2.115, 2.105, 2.095, &
    2.085, 2.075, 2.065, 2.055, 2.045, 2.035, 2.025, 2.015, 2.000, 1.985, &
    1.970, 1.955, 1.940, 1.925, 1.910, 1.895, 1.880, 1.865, 1.850, 1.835, &
    1.820, 1.805, 1.790, 1.775, 1.760, 1.745, 1.730, 1.715, 1.700, 1.685, &
    1.670, 1.655, 1.640, 1.625, 1.610, 1.595, 1.580, 1.565, 1.550, 1.535, &
    1.520, 1.505, 1.490, 1.475, 1.460, 1.445, 1.430, 1.415, 1.400, 1.385, &
    1.370, 1.355, 1.340, 1.325, 1.310, 1.295, 1.280, 1.265, 1.250, 1.235, &
    1.220, 1.205, 1.190, 1.175, 1.160, 1.145, 1.130, 1.115, 1.100, 1.085, &
    1.070, 1.055, 1.040, 1.025, 1.010, 0.995, 0.980, 0.965, 0.950, 0.935, &
    0.920, 0.905, 0.890, 0.875, 0.860, 0.845, 0.830, 0.815, 0.800, 0.785, &
    0.770, 0.755, 0.740, 0.725, 0.710, 0.695, 0.680, 0.665, 0.650, 0.635, &
    0.620, 0.605, 0.590, 0.575, 0.560, 0.545, 0.530, 0.515, 0.500, 0.485, &
    0.470, 0.455, 0.440, 0.425, 0.410, 0.395, 0.380, 0.365, 0.350, 0.335, &
    0.320, 0.305, 0.290, 0.275, 0.260, 0.245, 0.230, 0.215, 0.200, 0.185, &
    0.170, 0.155, 0.140, 0.125, 0.110, 0.100, 0.090, 0.080, 0.070, 0.060, &
    0.050, 0.045, 0.040, 0.038, 0.035, 0.032, 0.030, 0.028, 0.025, 0.022, &
    0.020, 0.018, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.009, &
    0.008, 0.007, 0.006, 0.005, 0.004, 0.003, 0.002, 0.001, 0.0005, 0.0002, &
    (0.0001, i=1, 110) ]
CONTAINS
  FUNCTION interpolar_am15g(e_mev) RESULT(flux)
    REAL(DP), INTENT(IN) :: e_mev
    REAL(DP) :: flux, wv_pt
    INTEGER :: idx
    wv_pt = 1240.0_DP / MAX(e_mev, 1.0E-3_DP)
    flux = 0.0_DP
    DO idx = 1, 399
       IF (wv_pt >= SOLAR_WV(idx) .AND. wv_pt < SOLAR_WV(idx+1)) THEN
          flux = SOLAR_IRR(idx) + (SOLAR_IRR(idx+1)-SOLAR_IRR(idx))*(wv_pt-SOLAR_WV(idx))/(SOLAR_WV(idx+1)-SOLAR_WV(idx))
          RETURN
       END IF
    END DO
  END FUNCTION interpolar_am15g
END MODULE dataset_am15g_explicit

MODULE algebra_linear_estendida
  USE constantes_quanticas
  IMPLICIT NONE
CONTAINS
  SUBROUTINE decomp_lu_complex_12x12(A, B, X)
    COMPLEX(DP), INTENT(IN) :: A(12,12), B(12)
    COMPLEX(DP), INTENT(OUT) :: X(12)
    COMPLEX(DP) :: LU(12,12), Y(12), tmp_row(12)
    REAL(DP) :: pivot_val
    INTEGER :: i, j, k, pivot_idx, P(12), tmp_p
    
    LU = A ; P = [(i, i=1,12)]
    DO k = 1, 11
       pivot_idx = k ; pivot_val = ABS(LU(k,k))
       DO i = k+1, 12
          IF (ABS(LU(i,k)) > pivot_val) THEN
             pivot_val = ABS(LU(i,k)) ; pivot_idx = i
          END IF
       END DO
       IF (pivot_idx /= k) THEN
          tmp_row = LU(k,:) ; LU(k,:) = LU(pivot_idx,:) ; LU(pivot_idx,:) = tmp_row
          tmp_p = P(k) ; P(k) = P(pivot_idx) ; P(pivot_idx) = tmp_p
       END IF
       DO i = k+1, 12
          LU(i,k) = LU(i,k) / LU(k,k)
          DO j = k+1, 12 ; LU(i,j) = LU(i,j) - LU(i,k)*LU(k,j) ; END DO
       END DO
    END DO
    DO i = 1, 12 ; Y(i) = B(P(i)) - SUM(LU(i, 1:i-1)*Y(1:i-1)) ; END DO
    DO i = 12, 1, -1 ; X(i) = (Y(i)-SUM(LU(i, i+1:12)*X(i+1:12)))/LU(i,i) ; END DO
  END SUBROUTINE decomp_lu_complex_12x12

  SUBROUTINE inverter_12x12_full(A, Ainv)
    COMPLEX(DP), INTENT(IN)  :: A(12,12)
    COMPLEX(DP), INTENT(OUT) :: Ainv(12,12)
    COMPLEX(DP) :: e_vec(12)
    INTEGER :: k
    DO k = 1, 12
       e_vec = 0.0_DP ; e_vec(k) = 1.0_DP
       CALL decomp_lu_complex_12x12(A, e_vec, Ainv(:,k))
    END DO
  END SUBROUTINE inverter_12x12_full
END MODULE algebra_linear_estendida

MODULE fisica_bir_pikus
  USE constantes_quanticas
  IMPLICIT NONE
CONTAINS
  SUBROUTINE obter_shifts_bp(eps, shifts)
    REAL(DP), INTENT(IN) :: eps(3,3)
    REAL(DP), INTENT(OUT) :: shifts(6)
    REAL(DP) :: tr_eps
    tr_eps = eps(1,1) + eps(2,2) + eps(3,3)
    shifts(1:2) = (1.1 * tr_eps + 8.6 * eps(1,1)) * 1000.0_DP
    shifts(3:4) = (1.1 * tr_eps + 8.6 * eps(2,2)) * 1000.0_DP
    shifts(5:6) = (1.1 * tr_eps + 8.6 * eps(3,3)) * 1000.0_DP
  END SUBROUTINE obter_shifts_bp
END MODULE fisica_bir_pikus

MODULE motor_fisica_si_p
  USE constantes_quanticas
  USE fisica_bir_pikus
  IMPLICIT NONE
CONTAINS
  SUBROUTINE montar_hamiltoniano(H, eps)
    COMPLEX(DP), INTENT(OUT) :: H(12,12)
    REAL(DP), INTENT(IN) :: eps(3,3)
    REAL(DP) :: shifts(6), e_ref
    INTEGER :: i, v, s_idx
    CALL obter_shifts_bp(eps, shifts)
    H = 0.0_DP
    DO i = 1, 12
       v = MOD(i-1, 6) + 1 ; s_idx = (i-1)/6 + 1
       SELECT CASE (v)
       CASE (1) ; e_ref = -45.0 ; CASE (2:4) ; e_ref = -33.0 ; CASE (5:6) ; e_ref = -31.0
       END SELECT
       H(i,i) = CMPLX(e_ref + shifts(v), (2*s_idx-3)*0.04, DP)
       DO v = 1, 12 ; IF (i /= v) H(i,v) = CMPLX(0.5, 0.0, DP) ; END DO
    END DO
  END SUBROUTINE montar_hamiltoniano
END MODULE motor_fisica_si_p

MODULE scba_eliashberg_engine
  USE constantes_quanticas
  IMPLICIT NONE
CONTAINS
  SUBROUTINE solve_scba(Gr, Gl, ne, E_vals, Sr, Sl)
    COMPLEX(DP), INTENT(IN)  :: Gr(12,12,ne), Gl(12,12,ne)
    REAL(DP),    INTENT(IN)  :: E_vals(ne)
    INTEGER,     INTENT(IN)  :: ne
    COMPLEX(DP), INTENT(OUT) :: Sr(12,12,ne), Sl(12,12,ne)
    INTEGER :: i, k
    REAL(DP) :: om, nq, g2
    Sr = 0.0_DP ; Sl = 0.0_DP
    DO i = 1, ne
       DO k = 1, ne
          om = ABS(E_vals(i) - E_vals(k))
          IF (om > 1.0 .AND. om < 80.0) THEN
             g2 = 6.0 * (om/60.0)**2 * EXP(-om/30.0)
             nq = 1.0 / (EXP(om/(K_B*T_REF)) - 1.0)
             Sl(:,:,i) = Sl(:,:,i) + g2 * (nq + 1.0) * Gl(:,:,k)
             Sr(:,:,i) = Sr(:,:,i) - CI * 0.5 * g2 * (nq + 0.5) * identity_matrix_12()
          END IF
       END DO
    END DO
  END SUBROUTINE solve_scba
END MODULE scba_eliashberg_engine

PROGRAM simulador_solar_master
  USE libnegf_types
  USE libnegf_physics
  USE constantes_quanticas
  USE dataset_am15g_explicit
  USE algebra_linear_estendida
  USE motor_fisica_si_p
  USE scba_eliashberg_engine
  IMPLICIT NONE
  
  TYPE(QuantumSystem) :: sys
  COMPLEX(DP), ALLOCATABLE :: Gr(:,:,:), Gl(:,:,:), Sr(:,:,:), Sl(:,:,:)
  COMPLEX(DP) :: H_eff(12,12), SigL(12,12), SigR(12,12), G_tmp(12,12)
  REAL(DP)    :: E(500), n_den(12), v_har(12), bias, i_net, i_sun, eps_ten(3,3), ldos(500)
  INTEGER     :: iv, ie, ip, st, ne = 500
  
  eps_ten = 0.0_DP ; eps_ten(1,1) = 0.002 ; eps_ten(2,2) = 0.002 ; eps_ten(3,3) = -0.001
  CALL sys%init(s=12, l=2)
  CALL montar_hamiltoniano(sys%ham, eps_ten)
  
  DO iv = 1, 2
     CALL sys%leads(iv)%init(d=12)
     sys%leads(iv)%h00 = 0.0_DP ; sys%leads(iv)%h01 = -18.0 * identity_matrix_12()
  END DO

  ALLOCATE(Gr(12,12,ne), Gl(12,12,ne), Sr(12,12,ne), Sl(12,12,ne))
  Gr = 0.0_DP ; Gl = 0.0_DP ; Sr = 0.0_DP ; Sl = 0.0_DP

  OPEN(UNIT=88, FILE='resultados_solar_atomico.dat', STATUS='REPLACE')
  DO iv = 1, 31
     bias = -0.1 + REAL(iv-1, DP) * 0.01
     v_har = 0.0_DP ; n_den = 0.0_DP
     DO ip = 1, 10
        DO ie = 1, ne
           E(ie) = -130.0 + REAL(ie-1, DP) * 0.5
           H_eff = sys%ham ; DO st=1,12 ; H_eff(st,st) = H_eff(st,st) + v_har(st) ; END DO
           CALL calculate_surface_green(CMPLX(E(ie), 1.0E-6, DP), sys%leads(1), SigL)
           CALL calculate_surface_green(CMPLX(E(ie), 1.0E-6, DP), sys%leads(2), SigR)
           G_tmp = CMPLX(E(ie), 1.0E-7, DP)*identity_matrix_12() - H_eff - SigL - SigR - Sr(:,:,ie)
           CALL inverter_12x12_full(G_tmp, Gr(:,:,ie))
           H_eff = CI*(SigL-CONJG(TRANSPOSE(SigL)))*f_d(E(ie), bias/2.0) + &
                   CI*(SigR-CONJG(TRANSPOSE(SigR)))*f_d(E(ie), -bias/2.0) + Sl(:,:,ie)
           Gl(:,:,ie) = MATMUL(Gr(:,:,ie), MATMUL(H_eff, CONJG(TRANSPOSE(Gr(:,:,ie)))))
           ldos(ie) = -AIMAG(SUM([(Gr(st,st,ie), st=1,12)])) / PI
        END DO
        CALL solve_scba(Gr, Gl, ne, E, Sr, Sl)
        DO ie = 1, 12 ; n_den(ie) = SUM(REAL(Gl(ie,ie,:)))*0.5 / (2.0*PI) ; END DO
        DO ie = 1, 12 ; v_har(ie) = 25.0 * (SUM(n_den)-n_den(ie)) ; END DO
     END DO
     i_net = SUM(REAL(Gl(1,1,:))) * 1.2E-7 ; i_sun = interpolar_am15g(15.0) * 5.0E-9 * (1.0-n_den(1))
     WRITE(88, *) bias, (i_net+i_sun)*1000.0
  END DO
  CLOSE(88)

CONTAINS
  FUNCTION f_d(ener, mu) RESULT(f)
    REAL(DP) :: ener, mu, f
    f = 1.0 / (1.0 + EXP((ener - mu) / (K_B * T_REF)))
  END FUNCTION f_d
END PROGRAM simulador_solar_master
