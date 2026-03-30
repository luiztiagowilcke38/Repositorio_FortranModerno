program main_climate
    use mod_climate_phys
    use mod_climate_dynamics
    implicit none
    
    integer, parameter :: NX = 120, NY = 60 ! Resolução 3 deg Lat-Lon
    integer, parameter :: YEARS = 100
    integer, parameter :: STEPS_PER_YEAR = 365
    real(real64), parameter :: DT = 86400.0_real64 ! 1 dia
    
    type(climate_state) :: earth
    integer :: i, j, t
    real(real64) :: co2_emission_rate
    
    ! --- 1. Inicialização ---
    allocate(earth%temp(NX, NY))
    allocate(earth%co2_conc(NX, NY))
    allocate(earth%albedo(NX, NY))
    allocate(earth%u(NX, NY))
    allocate(earth%v(NX, NY))
    
    earth%temp = 288.15_real64 ! Temperatura inicial pré-industrial
    earth%co2_conc = 280.0_real64 ! CO2 pré-industrial
    earth%albedo = 0.30_real64
    earth%u = 5.0_real64 ! Ventos simplificados (Zonal)
    earth%v = 0.0_real64
    
    co2_emission_rate = 0.05_real64 ! Taxa de emissão anual (ppm/year)
    
    print *, "CENÁRIO: Aquecimento Global por Forçante de CO2"
    print *, "T0 = ", earth%temp(1,1), " K"
    print *, "CO2_0 = ", earth%co2_conc(1,1), " ppm"
    
    ! --- 2. Loop de Séculos ---
    do t = 1, YEARS * STEPS_PER_YEAR
        
        ! Dinâmica Atmosférica: Transporte de Calor e CO2
        call advection_2d(earth%temp, earth%u, earth%v, 300000.0_real64, 300000.0_real64, DT, NX, NY)
        call advection_2d(earth%co2_conc, earth%u, earth%v, 300000.0_real64, 300000.0_real64, DT, NX, NY)
        
        ! Balanço Energético com Efeito Estufa
        call energy_balance(earth, DT, NX, NY)
        
        ! Injeção de CO2 (Emissão Antropogênica Simplificada)
        earth%co2_conc = earth%co2_conc + (co2_emission_rate / STEPS_PER_YEAR)
        
        ! Feedback de Gelo-Albedo (Simplificado)
        !$omp parallel do private(i, j)
        do j = 1, NY
            do i = 1, NX
                if (earth%temp(i,j) < 273.15_real64) then
                    earth%albedo(i,j) = 0.60_real64 ! Gelo
                else
                    earth%albedo(i,j) = 0.20_real64 ! Oceano/Solo
                end if
            end do
        end do
        !$omp end parallel do
        
        ! Log de Progresso Anual
        if (mod(t, STEPS_PER_YEAR) == 0) then
            print '(A,I4,A,F8.3,A,F8.2)', "Ano: ", t/STEPS_PER_YEAR, " | T_avg: ", &
                   earth%global_avg_temp, " K | CO2: ", earth%co2_conc(1,1), " ppm"
        end if
        
    end do
    
    print *, "Simulação Concluída."
    print *, "Delta T final: ", earth%global_avg_temp - 288.15_real64, " K"

end program main_climate
