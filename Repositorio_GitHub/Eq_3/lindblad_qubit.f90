module mod_lindblad
    use iso_fortran_env
    implicit none
    
    integer, parameter :: DIM = 2 ! Qubit (Sistema de 2 níveis)
    complex(real64), parameter :: H(DIM, DIM) = reshape([(.5_real64, 0._real64), (0._real64, 0._real64), &
                                                         (0._real64, 0._real64), (-.5_real64, 0._real64)], [DIM, DIM])
    real(real64), parameter :: GAMMA = 0.1_real64 ! Taxa de relaxação
    complex(real64), parameter :: L_OP(DIM, DIM) = reshape([(0._real64, 0._real64), (1._real64, 0._real64), &
                                                            (0._real64, 0._real64), (0._real64, 0._real64)], [DIM, DIM]) ! Sigma-minus
    
contains

    subroutine step_lindblad(rho, drho)
        complex(real64), intent(in)  :: rho(DIM, DIM)
        complex(real64), intent(out) :: drho(DIM, DIM)
        complex(real64) :: drho_h(DIM, DIM), drho_l(DIM, DIM)
        
        ! 1. Termo Hamiltoniano: -i[H, rho]
        drho_h = (0._real64, -1._real64) * (matmul(H, rho) - matmul(rho, H))
        
        ! 2. Termo Dissipativo: L*rho*L' - 0.5*{L'L, rho}
        drho_l = GAMMA * (matmul(matmul(L_OP, rho), conjg(transpose(L_OP))) - &
                 0.5_real64 * (matmul(matmul(conjg(transpose(L_OP)), L_OP), rho) + &
                               matmul(rho, matmul(conjg(transpose(L_OP)), L_OP))))
        
        drho = drho_h + drho_l
    end subroutine step_lindblad

    subroutine rk4_lindblad(rho, dt)
        complex(real64), intent(inout) :: rho(DIM, DIM)
        real(real64), intent(in) :: dt
        complex(real64) :: k1(DIM, DIM), k2(DIM, DIM), k3(DIM, DIM), k4(DIM, DIM), rtmp(DIM, DIM)
        
        call step_lindblad(rho, k1)
        rtmp = rho + 0.5_real64 * dt * k1
        call step_lindblad(rtmp, k2)
        rtmp = rho + 0.5_real64 * dt * k2
        call step_lindblad(rtmp, k3)
        rtmp = rho + dt * k3
        call step_lindblad(rtmp, k4)
        
        rho = rho + (dt/6.0_real64) * (k1 + 2.0_real64*k2 + 2.0_real64*k3 + k4)
    end subroutine rk4_lindblad

end module mod_lindblad

program main_lindblad
    use mod_lindblad
    implicit none
    complex(real64) :: rho(DIM, DIM)
    integer :: t
    real(real64) :: dt = 0.01_real64
    
    ! Estado inicial: Qubit excitado |1><1|
    rho = reshape([(0._real64, 0._real64), (0._real64, 0._real64), &
                   (0._real64, 0._real64), (1._real64, 0._real64)], [DIM, DIM])
    
    print *, "# Simulando Equação de Lindblad (Sistemas Quânticos Abertos)"
    print *, "T | Excitação (rho22)"
    
    do t = 1, 1000
        call rk4_lindblad(rho, dt)
        if (mod(t, 100) == 0) then
            print '(F8.2, A, F8.4)', t*dt, " | ", abs(rho(2,2))
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_lindblad
