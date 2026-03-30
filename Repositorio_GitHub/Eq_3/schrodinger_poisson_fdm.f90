module mod_schrodinger_poisson
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 64
    real(real64), parameter :: L = 20.0_real64
    real(real64), parameter :: DX = L / N
    real(real64), parameter :: DT = 0.005_real64
    
contains

    ! Solver de Poisson via Relaxação (ou SOR)
    ! Lap(phi) = |psi|^2 - mean(|psi|^2)
    subroutine solve_poisson(psi, phi)
        complex(real64), intent(in) :: psi(N, N, N)
        real(real64), intent(inout) :: phi(N, N, N)
        real(real64) :: rho(N, N, N), mean_rho
        integer :: i, j, k, im1, ip1, jm1, jp1, km1, kp1, iter
        
        rho = abs(psi)**2
        mean_rho = sum(rho) / real(N**3, real64)
        rho = rho - mean_rho
        
        do iter = 1, 100 ! Iterações de relaxação Jacobi
            do k = 1, N
                do j = 1, N
                    do i = 1, N
                        im1 = mod(i-2+N, N) + 1; ip1 = mod(i, N) + 1
                        jm1 = mod(j-2+N, N) + 1; jp1 = mod(j, N) + 1
                        km1 = mod(k-2+N, N) + 1; kp1 = mod(k, N) + 1
                        
                        phi(i, j, k) = (phi(ip1, j, k) + phi(im1, j, k) + &
                                        phi(i, jp1, k) + phi(i, jm1, k) + &
                                        phi(i, j, kp1) + phi(i, j, km1) - DX**2 * rho(i, j, k)) / 6.0_real64
                    end do
                end do
            end do
        end do
    end subroutine solve_poisson

    ! i*psi_t = -0.5*Lap(psi) + phi*psi
    subroutine step_sp(psi, phi, dpsi)
        complex(real64), intent(in)  :: psi(N, N, N)
        real(real64), intent(in)     :: phi(N, N, N)
        complex(real64), intent(out) :: dpsi(N, N, N)
        complex(real64) :: lap_psi
        integer :: i, j, k, im1, ip1, jm1, jp1, km1, kp1
        
        do k = 1, N
            do j = 1, N
                do i = 1, N
                    im1 = mod(i-2+N, N) + 1; ip1 = mod(i, N) + 1
                    jm1 = mod(j-2+N, N) + 1; jp1 = mod(j, N) + 1
                    km1 = mod(k-2+N, N) + 1; kp1 = mod(k, N) + 1
                    
                    lap_psi = (psi(ip1, j, k) + psi(im1, j, k) + &
                               psi(i, jp1, k) + psi(i, jm1, k) + &
                               psi(i, j, kp1) + psi(i, j, km1) - 6.0_real64 * psi(i, j, k)) / (DX**2)
                    
                    dpsi(i, j, k) = (0.0_real64, -1.0_real64) * (-0.5_real64 * lap_psi + phi(i, j, k) * psi(i, j, k))
                end do
            end do
        end do
    end subroutine step_sp

end module mod_schrodinger_poisson

program main_schrodinger_poisson
    use mod_schrodinger_poisson
    implicit none
    complex(real64), allocatable :: psi(:, :, :)
    real(real64), allocatable :: phi(:, :, :)
    integer :: t, i, j, k
    
    allocate(psi(N, N, N), phi(N, N, N))
    
    ! Condição inicial: Gaussiana (Halo de Matéria Escura)
    psi = (0.0_real64, 0.0_real64)
    phi = 0.0_real64
    do k = 1, N
        do j = 1, N
            do i = 1, N
                psi(i, j, k) = exp(-((i-N/2)**2 + (j-N/2)**2 + (k-N/2)**2)*DX**2 / 5.0_real64)
            end do
        end do
    end do
    
    print *, "# Simulando Equação de Schrödinger-Poisson (Fuzzy Dark Matter)"
    
    do t = 1, 1000
        call solve_poisson(psi, phi)
        block
            complex(real64) :: dpsi(N, N, N)
            call step_sp(psi, phi, dpsi)
            psi = psi + DT * dpsi
        end block
        
        if (mod(t, 200) == 0) then
            print '(A,I6,A,F12.6)', "Step: ", t, " | Central Density: ", abs(psi(N/2, N/2, N/2))**2
        end if
    end do
    
    print *, "# Simulação de Fuzzy Dark Matter Concluída."
end program main_schrodinger_poisson
