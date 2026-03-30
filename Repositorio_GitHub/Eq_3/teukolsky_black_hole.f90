module mod_teukolsky_radial
    use iso_fortran_env
    implicit none
    
    real(real64), parameter :: M_BH = 1.0_real64  ! Massa do Buraco Negro
    real(real64), parameter :: A_SPIN = 0.9_real64 ! Spin (Kerr)
    real(real64), parameter :: OMEGA_W = 0.5_real64 ! Frequência da onda
    integer, parameter :: S_PARAM = -2             ! Gravitacional
    
contains

    ! Função Delta de Kerr: r^2 - 2Mr + a^2
    function get_delta(r) result(d)
        real(real64), intent(in) :: r
        real(real64) :: d
        d = r**2 - 2.0_real64*M_BH*r + A_SPIN**2
    end function get_delta

    ! Equação Radial de Teukolsky: d2R/dr2 + P(r)dR/dr + Q(r)R = 0
    ! Convertida para sistema de 1a ordem para integração RK4
    subroutine teukolsky_derivs(r, y, dy)
        real(real64), intent(in)  :: r
        complex(real64), intent(in) :: y(2) ! y(1)=R, y(2)=dR/dr
        complex(real64), intent(out) :: dy(2)
        real(real64) :: d, d_prime, k_func
        complex(real64) :: p_term, q_term
        
        d = get_delta(r)
        d_prime = 2.0_real64*r - 2.0_real64*M_BH
        k_func = (r**2 + A_SPIN**2)*OMEGA_W ! Simplificado (m=0)
        
        ! P(z) e Q(z) da Equação de Teukolsky
        p_term = (S_PARAM + 1.0_real64) * d_prime / d
        q_term = (k_func**2 - (0.0_real64, 2.0_real64)*S_PARAM*(r - M_BH)*k_func) / d + &
                 (0.0_real64, 4.0_real64)*S_PARAM*OMEGA_W*r - 10.0_real64 ! lambda_s ~ 10
        
        dy(1) = y(2)
        dy(2) = -p_term * y(2) - q_term * y(1)
    end subroutine teukolsky_derivs

end module mod_teukolsky_radial

program main_teukolsky
    use mod_teukolsky_radial
    implicit none
    complex(real64) :: y(2), k1(2), k2(2), k3(2), k4(2), ytmp(2)
    real(real64) :: r, dr
    integer :: step
    
    ! Condições Iniciais no Horizonte de Eventos (ou próximo)
    r = 2.0_real64 * M_BH + 0.1_real64 ! Logo fora do horizonte
    dr = 0.01_real64
    y(1) = (1.0_real64, 0.0_real64)    ! Amplitude inicial
    y(2) = (0.0_real64, 0.5_real64)    ! Gradiente inicial
    
    print *, "# Simulando Equação Radial de Teukolsky (Perturbações de Buracos Negros de Kerr)"
    print *, "R | Amplitude (abs(R))"
    
    do step = 1, 1000
        call teukolsky_derivs(r, y, k1)
        ytmp = y + 0.5_real64 * dr * k1
        call teukolsky_derivs(r + 0.5_real64*dr, ytmp, k2)
        ytmp = y + 0.5_real64 * dr * k2
        call teukolsky_derivs(r + 0.5_real64*dr, ytmp, k3)
        ytmp = y + dr * k3
        call teukolsky_derivs(r + dr, ytmp, k4)
        
        y = y + (dr/6.0_real64) * (k1 + 2.0_real64*k2 + 2.0_real64*k3 + k4)
        r = r + dr
        
        if (mod(step, 200) == 0) then
            print '(F8.2, A, F12.6)', r, " | ", abs(y(1))
        end if
    end do
    
    print *, "# Simulação de Teukolsky Concluída."
end program main_teukolsky
