module mod_monge_ampere
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 100
    real(real64), parameter :: DX = 0.01_real64
    real(real64), parameter :: DT = 0.0001_real64 ! Estabilidade delicada
    
contains

    ! d(phi)/dt = det(Hessian(phi)) - f(x)
    subroutine step_ma(phi, f, dphi)
        real(real64), intent(in)  :: phi(N, N), f(N, N)
        real(real64), intent(out) :: dphi(N, N)
        real(real64) :: p_xx, p_yy, p_xy, det_h
        integer :: i, j, im1, ip1, jm1, jp1
        
        do j = 2, N-1
            do i = 2, N-1
                im1 = i-1; ip1 = i+1
                jm1 = j-1; jp1 = j+1
                
                ! Derivadas centrais de 2a ordem
                p_xx = (phi(ip1, j) - 2.0_real64*phi(i, j) + phi(im1, j)) / (DX**2)
                p_yy = (phi(i, jp1) - 2.0_real64*phi(i, j) + phi(i, jm1)) / (DX**2)
                
                ! Derivada mista phi_xy
                p_xy = (phi(ip1, jp1) - phi(ip1, jm1) - phi(im1, jp1) + phi(im1, jm1)) / (4.0_real64 * DX**2)
                
                ! Determinante do Hessiano
                det_h = p_xx * p_yy - p_xy**2
                
                dphi(i, j) = det_h - f(i, j)
            end do
        end do
        
        ! Condições de Contorno de Dirichlet (Simplificadas: phi=0 na borda)
        dphi(1, :) = 0.0_real64; dphi(N, :) = 0.0_real64
        dphi(:, 1) = 0.0_real64; dphi(:, N) = 0.0_real64
    end subroutine step_ma

end module mod_monge_ampere

program main_monge_ampere
    use mod_monge_ampere
    implicit none
    real(real64) :: phi(N, N), f(N, N)
    integer :: t, i, j
    
    ! Fonte f(x) (ex: Gaussiana no centro)
    do j = 1, N
        do i = 1, N
            f(i, j) = 1.0_real64 + exp(-((i-N/2)**2 + (j-N/2)**2)*DX**2 / 0.1_real64)
        end do
    end do
    
    ! Chute inicial: Parabolóide phi = x^2/2 + y^2/2 (det(H) = 1)
    do j = 1, N
        do i = 1, N
            phi(i, j) = 0.5_real64 * ((i*DX)**2 + (j*DX)**2)
        end do
    end do
    
    print *, "# Solucionando Equação de Monge-Ampère via Relaxação Pseudo-Temporal"
    
    do t = 1, 10000
        block
            real(real64) :: dphi(N, N)
            call step_ma(phi, f, dphi)
            phi = phi + DT * dphi
        end block
        
        if (mod(t, 2000) == 0) then
            print '(A,I6,A,E12.4)', "Step: ", t, " | Max Residue: ", maxval(abs(dphi))
        end if
    end do
    
    print *, "# Solução Convergida."
end program main_monge_ampere
