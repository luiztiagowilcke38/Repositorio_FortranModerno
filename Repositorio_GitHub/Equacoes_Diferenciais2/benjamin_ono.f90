module mod_benjamin_ono
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 256
    real(real64), parameter :: L = 100.0_real64
    real(real64), parameter :: DX = L / N
    real(real64), parameter :: DT = 0.01_real64
    
contains

    ! Hilbert Transform simplificada via soma direta (Elite Numerical Approach)
    ! Em produção, usaria FFT: H(f) = -i * sign(k) * f(k)
    subroutine discrete_hilbert(f, hf)
        real(real64), intent(in)  :: f(N)
        real(real64), intent(out) :: hf(N)
        integer :: i, j
        real(real64) :: sum_h
        
        do i = 1, N
            sum_h = 0.0_real64
            do j = 1, N
                if (i /= j) then
                    sum_h = sum_h + f(j) / tan(3.14159265_real64 * (i - j) / N)
                end if
            end do
            hf(i) = sum_h / N
        end do
    end subroutine discrete_hilbert

    subroutine step_bo(u, du)
        real(real64), intent(in)  :: u(N)
        real(real64), intent(out) :: du(N)
        real(real64) :: ux(N), u_xx(N), h_uxx(N)
        integer :: i, im1, ip1
        
        do i = 1, N
            im1 = mod(i - 2 + N, N) + 1
            ip1 = mod(i, N) + 1
            ux(i) = (u(ip1) - u(im1)) / (2.0_real64 * DX)
            u_xx(i) = (u(ip1) - 2.0_real64*u(i) + u(im1)) / (DX**2)
        end do
        
        call discrete_hilbert(u_xx, h_uxx)
        
        ! du/dt = -u*ux - H[u_xx]
        du = -u * ux - h_uxx
    end subroutine step_bo

end module mod_benjamin_ono

program main_bo
    use mod_benjamin_ono
    implicit none
    real(real64) :: u(N)
    integer :: t, i
    
    ! Condição inicial: Soliton de Benjamin-Ono (Lorentziano)
    ! u(x,0) = 4*c / (1 + c^2 * x^2)
    do i = 1, N
        u(i) = 4.0_real64 / (1.0_real64 + (0.5_real64 * (i - N/2) * DX)**2)
    end do
    
    print *, "# Simulando Equação de Benjamin-Ono (Ondas Internas Profundas)"
    
    do t = 1, 1000
        block
            real(real64) :: du_vec(N)
            call step_bo(u, du_vec)
            u = u + DT * du_vec
        end block
        
        if (mod(t, 200) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | u_max: ", maxval(u)
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_bo
