module mod_gross_pitaevskii
    use iso_fortran_env
    implicit none
    
    integer, parameter :: N = 256
    real(real64), parameter :: L = 20.0_real64
    real(real64), parameter :: DX = L / N
    real(real64), parameter :: DT = 0.001_real64
    real(real64), parameter :: G_PARAM = 10.0_real64 ! Interação repulsiva
    
contains

    subroutine step_gp(psi, dpsi)
        complex(real64), intent(in)  :: psi(N)
        complex(real64), intent(out) :: dpsi(N)
        integer :: i, im1, ip1
        real(real64) :: x, v_ext
        complex(real64) :: lap_psi
        
        do i = 1, N
            x = (i - N/2) * DX
            v_ext = 0.5_real64 * x**2 ! Potencial Harmônico
            
            im1 = mod(i - 2 + N, N) + 1
            ip1 = mod(i, N) + 1
            
            lap_psi = (psi(ip1) - 2.0_real64*psi(i) + psi(im1)) / (DX**2)
            
            ! i*dpsi/dt = -0.5*Lap(psi) + V*psi + g*|psi|^2*psi
            ! dpsi/dt = -i * [-0.5*Lap(psi) + V*psi + g*|psi|^2*psi]
            dpsi(i) = (0.0_real64, -1.0_real64) * &
                      (-0.5_real64 * lap_psi + v_ext * psi(i) + G_PARAM * abs(psi(i))**2 * psi(i))
        end do
    end subroutine step_gp

    subroutine rk4_gp(psi)
        complex(real64), intent(inout) :: psi(N)
        complex(real64) :: k1(N), k2(N), k3(N), k4(N), ptmp(N)
        
        call step_gp(psi, k1)
        ptmp = psi + 0.5_real64 * DT * k1
        call step_gp(ptmp, k2)
        ptmp = psi + 0.5_real64 * DT * k2
        call step_gp(ptmp, k3)
        ptmp = psi + DT * k3
        call step_gp(ptmp, k4)
        
        psi = psi + (DT/6.0_real64) * (k1 + 2.0_real64*k2 + 2.0_real64*k3 + k4)
    end subroutine rk4_gp

end module mod_gross_pitaevskii

program main_gp
    use mod_gross_pitaevskii
    implicit none
    complex(real64) :: psi(N)
    integer :: t, i
    real(real64) :: norm
    
    ! Condição inicial: Gaussiana (aproximação do estado fundamental)
    do i = 1, N
        psi(i) = exp(-((i - N/2) * DX)**2 / 2.0_real64)
    end do
    
    ! Normalização
    norm = sqrt(sum(abs(psi)**2) * DX)
    psi = psi / norm
    
    print *, "# Simulando Equação de Gross-Pitaevskii (BEC em Armadilha Harmônica)"
    
    do t = 1, 5000
        call rk4_gp(psi)
        if (mod(t, 1000) == 0) then
            norm = sum(abs(psi)**2) * DX
            print '(A,I6,A,F8.4)', "Step: ", t, " | Norm: ", norm
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_gp
