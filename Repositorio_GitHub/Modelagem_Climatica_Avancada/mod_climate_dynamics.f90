module mod_climate_dynamics
    use iso_fortran_env
    use mod_climate_phys
    implicit none
    
contains

    subroutine advection_2d(field, u, v, dx, dy, dt, nx, ny)
        real(real64), intent(inout) :: field(nx, ny)
        real(real64), intent(in) :: u(nx, ny), v(nx, ny), dx, dy, dt
        integer, intent(in) :: nx, ny
        real(real64) :: field_new(nx, ny)
        integer :: i, j
        
        ! Esquema Upwind de Primeira Ordem (Estabilidade Numérica)
        !$omp parallel do private(i, j)
        do j = 2, ny-1
            do i = 2, nx-1
                field_new(i, j) = field(i, j) - &
                    dt * ( max(u(i,j), 0.0_real64) * (field(i,j) - field(i-1,j))/dx + &
                           min(u(i,j), 0.0_real64) * (field(i+1,j) - field(i,j))/dx + &
                           max(v(i,j), 0.0_real64) * (field(i,j) - field(i,j-1))/dy + &
                           min(v(i,j), 0.0_real64) * (field(i,j+1) - field(i,j))/dy )
            end do
        end do
        !$omp end parallel do
        
        ! Condições de Contorno Periódicas (Simplificadas)
        field_new(1, :) = field_new(nx-1, :)
        field_new(nx, :) = field_new(2, :)
        field_new(:, 1) = field_new(:, 2)
        field_new(:, ny) = field_new(:, ny-1)
        
        field = field_new
    end subroutine advection_2d

    subroutine energy_balance(state, dt, nx, ny)
        type(climate_state), intent(inout) :: state
        real(real64), intent(in) :: dt
        integer, intent(in) :: nx, ny
        real(real64) :: S_inc, R_out, epsilon
        integer :: i, j
        
        ! Modelo Radiativo de Duas Camadas Simplificado
        ! epsilon depende da concentração de CO2 (Feedback Logarítmico)
        !$omp parallel do private(i, j, S_inc, R_out, epsilon)
        do j = 1, ny
            do i = 1, nx
                ! Radiação Solar Incidente ponderada pela latitude (aproximada)
                S_inc = (S0/4.0_real64) * (1.0_real64 - state%albedo(i,j))
                
                ! Emissividade Atmosférica (Efeito Estufa)
                ! Parametrização: epsilon = eps0 + beta * log(C/C0)
                epsilon = 0.62_real64 + 0.05_real64 * log(state%co2_conc(i,j) / 280.0_real64)
                
                ! Balanço: dT/dt = (S_inc - sigma * epsilon * T^4) / C_heat
                ! Aqui usamos uma capacidade térmica efetiva simplificada (oceano misto)
                R_out = SIGMA * epsilon * (state%temp(i,j)**4)
                
                state%temp(i,j) = state%temp(i,j) + dt * (S_inc - R_out) / 1.0e8_real64 ! C_eff ~ 10^8 J/m^2/K
            end do
        end do
        !$omp end parallel do
        
        state%global_avg_temp = sum(state%temp) / (nx*ny)
    end subroutine energy_balance

end module mod_climate_dynamics
