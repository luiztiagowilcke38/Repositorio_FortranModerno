module mod_levy_jump
    use iso_fortran_env
    implicit none
    
    integer, parameter :: NS = 200
    real(real64), parameter :: S_MAX = 200.0_real64
    real(real64), parameter :: DS = S_MAX / NS
    real(real64), parameter :: T_MAX = 1.0_real64
    real(real64), parameter :: DT = 0.001_real64
    
    real(real64), parameter :: SIGMA = 0.2_real64
    real(real64), parameter :: R = 0.05_real64
    real(real64), parameter :: LAMBDA = 0.1_real64 ! Taxa de salto (Lévy Intensity)
    real(real64), parameter :: MU_J = -0.1_real64 ! Média do salto log-normal
    real(real64), parameter :: SIGMA_J = 0.1_real64 ! Volatilidade do salto
    
contains

    subroutine step_jump_diffusion(v, dv)
        real(real64), intent(in)  :: v(0:NS)
        real(real64), intent(out) :: dv(0:NS)
        integer :: i, j
        real(real64) :: s, vs, vss, jump_int
        
        do i = 1, NS-1
            s = i * DS
            vs = (v(i+1) - v(i-1)) / (2.0_real64 * DS)
            vss = (v(i+1) - 2.0_real64*v(i) + v(i-1)) / (DS**2)
            
            ! Cálculo simplificado da integral de salto (Lévy Jump)
            jump_int = 0.0_real64
            do j = 1, NS-1
                ! Distribuição Log-Normal do Salto
                jump_int = jump_int + v(j) * exp(-(log(real(j,real64)/i) - MU_J)**2 / (2.0_real64*SIGMA_J**2))
            end do
            jump_int = jump_int / (SIGMA_J * sqrt(2.0_real64*3.14159_real64))
            
            ! dV/dt = -(0.5*sigma^2*S^2*Vss + (r-lambda*k)*S*Vs - (r+lambda)*V + lambda*Integral)
            dv(i) = 0.5_real64 * SIGMA**2 * s**2 * vss + (R - 0.05_real64) * s * vs - &
                    (R + LAMBDA) * v(i) + LAMBDA * jump_int / NS
        end do
        
        dv(0) = 0.0_real64
        dv(NS) = 0.0_real64
    end subroutine step_jump_diffusion

end module mod_levy_jump

program main_levy
    use mod_levy_jump
    implicit none
    real(real64) :: v(0:NS), dv(0:NS)
    integer :: t, i
    real(real64) :: strike = 100.0_real64
    
    ! Condição final (Payoff de Opção de Compra - Call)
    do i = 0, NS
        v(i) = max(i * DS - strike, 0.0_real64)
    end do
    
    print *, "# Simulando Modelo de Merton (Jump-Diffusion / Lévy Process)"
    
    ! Integração reversa no tempo (Solução de Black-Scholes PIDE)
    do t = 1, 500
        call step_jump_diffusion(v, dv)
        v = v + DT * dv
        if (mod(t, 100) == 0) then
            print '(A,I6,A,F8.4)', "Step: ", t, " | At-The-Money Value: ", v(NS/2)
        end if
    end do
    
    print *, "# Simulação concluída."
end program main_levy
