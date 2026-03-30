module mod_dirac_graphene
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 128
    real(real64), parameter :: L = 50.0_real64
    real(real64), parameter :: DX = L / N
    real(real64), parameter :: DT = 0.005_real64
    real(real64), parameter :: VF = 1.0_real64 ! Velocidade de Fermi unitária
    
contains

    subroutine step_dirac(psi1, psi2, dpsi1, dpsi2)
        complex(real64), intent(in)  :: psi1(N, N), psi2(N, N)
        complex(real64), intent(out) :: dpsi1(N, N), dpsi2(N, N)
        integer :: i, j, im1, ip1, jm1, jp1
        complex(real64) :: p1_x, p1_y, p2_x, p2_y
        
        do j = 1, N
            do i = 1, N
                im1 = mod(i - 2 + N, N) + 1
                ip1 = mod(i, N) + 1
                jm1 = mod(j - 2 + N, N) + 1
                jp1 = mod(j, N) + 1
                
                ! Derivadas centrais
                p1_x = (psi1(ip1, j) - psi1(im1, j)) / (2.0_real64 * DX)
                p1_y = (psi1(i, jp1) - psi1(i, jm1)) / (2.0_real64 * DX)
                p2_x = (psi2(ip1, j) - psi2(im1, j)) / (2.0_real64 * DX)
                p2_y = (psi2(i, jp1) - psi2(i, jm1)) / (2.0_real64 * DX)
                
                ! dpsi1/dt = -VF * (dx - i*dy) psi2
                ! dpsi2/dt = -VF * (dx + i*dy) psi1
                dpsi1(i, j) = -VF * (p2_x - (0.0_real64, 1.0_real64) * p2_y)
                dpsi2(i, j) = -VF * (p1_x + (0.0_real64, 1.0_real64) * p1_y)
            end do
        end do
    end subroutine step_dirac

end module mod_dirac_graphene

program main_dirac
    use mod_dirac_graphene
    implicit none
    complex(real64) :: psi1(N, N), psi2(N, N)
    complex(real64) :: k1_1(N, N), k1_2(N, N), k2_1(N, N), k2_2(N, N)
    integer :: t, i, j
    
    ! Condição inicial: Pacote de onda Gaussiano em psi1
    psi1 = (0.0_real64, 0.0_real64)
    psi2 = (0.0_real64, 0.0_real64)
    do j = 1, N
        do i = 1, N
            psi1(i, j) = exp(-((i-N/2)**2 + (j-N/2)**2) * DX**2 / 10.0_real64)
        end do
    end do
    
    print *, "# Simulando Equação de Dirac (Férmions de Dirac em Grafeno)"
    
    do t = 1, 1000
        ! Integrador simples de Euler (para demonstração)
        ! Em um modelo real de elite, usaríamos RK4 ou Split-Step
        call step_dirac(psi1, psi2, k1_1, k1_2)
        psi1 = psi1 + DT * k1_1
        psi2 = psi2 + DT * k1_2
        
        if (mod(t, 200) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | Probability: ", sum(abs(psi1)**2 + abs(psi2)**2) * DX**2
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_dirac
