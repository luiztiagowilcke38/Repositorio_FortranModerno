module mod_camassa_holm
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 512
    real(real64), parameter :: L = 100.0_real64
    real(real64), parameter :: DX = L / N
    real(real64), parameter :: DT = 0.01_real64
    
contains

    ! Inversão de Helmholtz: m = u - u_xx
    ! Resolve o sistema tridiagonal (1 - 1/DX^2)u_i + (1/DX^2)u_{i-1} + (1/DX^2)u_{i+1} = m_i
    subroutine helmholtz_invert(m, u)
        real(real64), intent(in)  :: m(N)
        real(real64), intent(out) :: u(N)
        real(real64) :: a(N), b(N), c(N), d(N), cp(N), dp(N)
        real(real64) :: m_diag
        integer :: i
        
        m_diag = 1.0_real64 + 2.0_real64 / (DX**2)
        b = m_diag
        a = -1.0_real64 / (DX**2)
        c = a
        d = m
        
        ! Algoritmo de Thomas (TDMA)
        cp(1) = c(1) / b(1)
        dp(1) = d(1) / b(1)
        do i = 2, N
            cp(i) = c(i) / (b(i) - a(i) * cp(i-1))
            dp(i) = (d(i) - a(i) * dp(i-1)) / (b(i) - a(i) * cp(i-1))
        end do
        
        u(N) = dp(N)
        do i = N-1, 1, -1
            u(i) = dp(i) - cp(i) * u(i+1)
        end do
    end subroutine helmholtz_invert

    ! dm/dt = -(u*mx + 2*m*ux)
    subroutine step_ch(m, u, dm)
        real(real64), intent(in)  :: m(N), u(N)
        real(real64), intent(out) :: dm(N)
        real(real64) :: mx, ux
        integer :: i, im1, ip1
        
        do i = 1, N
            im1 = mod(i - 2 + N, N) + 1
            ip1 = mod(i, N) + 1
            
            mx = (m(ip1) - m(im1)) / (2.0_real64 * DX)
            ux = (u(ip1) - u(im1)) / (2.0_real64 * DX)
            
            dm(i) = -(u(i) * mx + 2.0_real64 * m(i) * ux)
        end do
    end subroutine step_ch

end module mod_camassa_holm

program main_camassa_holm
    use mod_camassa_holm
    implicit none
    real(real64) :: m(N), u(N), dm(N)
    integer :: t, i
    
    ! Condição inicial: Multi-peakons u(x) = sum c_i exp(-|x-x_i|)
    ! Aqui m(x) = sum 2*c_i delta(x-x_i) -> aproximada por gaussianas estreitas
    m = 0.0_real64
    do i = 1, N
        m(i) = 2.0_real64 * exp(-((i-N/4)*DX)**2 / 0.5_real64) + &
               1.0_real64 * exp(-((i-N/2)*DX)**2 / 0.5_real64)
    end do
    
    print *, "# Simulando Equação de Camassa-Holm (Dinâmica de Peakons)"
    
    do t = 1, 1000
        call helmholtz_invert(m, u)
        call step_ch(m, u, dm)
        m = m + DT * dm
        
        if (mod(t, 200) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | Max Amplitude: ", maxval(u)
        end if
    end do
    
    print *, "# Simulação de Camassa-Holm Concluída."
end program main_camassa_holm
