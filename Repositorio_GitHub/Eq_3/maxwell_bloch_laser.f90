module mod_maxwell_bloch
    use iso_fortran_env
    implicit none
    
    integer, parameter :: NZ = 200
    real(real64), parameter :: DZ = 0.1_real64
    real(real64), parameter :: DT = 0.01_real64
    
    ! Parâmetros Físicos (Semi-Clássicos)
    real(real64), parameter :: KAPPA = 1.0_real64  ! Acoplamento Campo-Matéria
    real(real64), parameter :: GAMMA_P = 0.5_real64 ! Taxa de desfasagem (T2)
    real(real64), parameter :: GAMMA_D = 0.1_real64 ! Taxa de relaxação (T1)
    real(real64), parameter :: DELTA = 0.05_real64  ! Detuning (Desafinaçao)
    real(real64), parameter :: OMEGA = 1.0_real64  ! Momento de Transição
    
contains

    subroutine step_mb(e, p, d_pop, de, dp, dd)
        complex(real64), intent(in)  :: e(NZ), p(NZ)
        real(real64), intent(in)     :: d_pop(NZ)
        complex(real64), intent(out) :: de(NZ), dp(NZ)
        real(real64), intent(out)    :: dd(NZ)
        integer :: z
        complex(real64) :: pulse_grad
        
        do z = 1, NZ
            ! 1. Evolução do Campo Elétrico (Propagação)
            ! dE/dz + (1/c)dE/dt = i*kappa*P
            ! Aqui simplificamos para dE/dt = -c*dE/dz + i*kappa*P
            if (z > 1) then
                pulse_grad = (e(z) - e(z-1)) / DZ
            else
                pulse_grad = (0.0_real64, 0.0_real64)
            end if
            
            de(z) = -pulse_grad + (0.0_real64, 1.0_real64) * KAPPA * p(z)
            
            ! 2. Evolução da Polarização P
            ! dP/dt = -(i*Delta + gamma_p)*P + i*Omega*E*D
            dp(z) = -((0.0_real64, 1.0_real64) * DELTA + GAMMA_P) * p(z) + &
                    (0.0_real64, 1.0_real64) * OMEGA * e(z) * d_pop(z)
            
            ! 3. Evolução da Inversão de População D
            ! dD/dt = -gamma_d*(D - D_eq) + i*Omega*(E* * P - E * P*)
            dd(z) = -GAMMA_D * (d_pop(z) - (-1.0_real64)) + &
                    real((0.0_real64, 1.0_real64) * OMEGA * (conjg(e(z))*p(z) - e(z)*conjg(p(z))))
        end do
    end subroutine step_mb

end module mod_maxwell_bloch

program main_maxwell_bloch
    use mod_maxwell_bloch
    implicit none
    complex(real64) :: e(NZ), p(NZ)
    real(real64) :: d_pop(NZ)
    integer :: t, z
    
    ! Condições Iniciais: Átomo no estado fundamental (D = -1, P = 0)
    d_pop = -1.0_real64
    p = (0.0_real64, 0.0_real64)
    e = (0.0_real64, 0.0_real64)
    
    ! Pulso de entrada Laser em z=1
    e(1) = (1.0_real64, 0.0_real64)
    
    print *, "# Simulando Equações de Maxwell-Bloch (Interação Laser-Matéria)"
    
    do t = 1, 1000
        block
            complex(real64) :: de(NZ), dp(NZ)
            real(real64) :: dd(NZ)
            call step_mb(e, p, d_pop, de, dp, dd)
            e = e + DT * de
            p = p + DT * dp
            d_pop = d_pop + DT * dd
        end block
        
        if (mod(t, 200) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | Field at Center: ", abs(e(NZ/2))
        end if
    end do
    
    print *, "# Simulação de Propagação de Pulso Óptico Concluída."
end program main_maxwell_bloch
