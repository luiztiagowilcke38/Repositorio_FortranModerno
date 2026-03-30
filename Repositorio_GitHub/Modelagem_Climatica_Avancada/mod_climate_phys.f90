module mod_climate_phys
    use iso_fortran_env
    implicit none
    
    ! Constantes Físicas Universais (SI)
    real(real64), parameter :: SIGMA = 5.670374e-8_real64 ! Stefan-Boltzmann [W/m^2/K^4]
    real(real64), parameter :: S0    = 1361.0_real64     ! Constante Solar [W/m^2]
    real(real64), parameter :: ALBEDO_GLOBAL = 0.30_real64 ! Albedo médio da Terra
    real(real64), parameter :: R_EARTH = 6371000.0_real64 ! Raio da Terra [m]
    real(real64), parameter :: GRAVITY = 9.80665_real64   ! Gravidade [m/s^2]
    real(real64), parameter :: CP_AIR  = 1004.6_real64    ! Calor específico do ar [J/kg/K]
    real(real64), parameter :: R_GAS   = 287.05_real64    ! Constante do ar seco [J/kg/K]
    
    ! Parâmetros do Ciclo do Carbono
    real(real64), parameter :: PPM_TO_KG = 7.76e12_real64 ! Conversão ppm -> kg CO2 (aproximado)
    
    type :: grid_type
        integer :: nx, ny, nz
        real(real64), allocatable :: lat(:), lon(:), h(:)
        real(real64) :: dx, dy, dz
    end type grid_type
    
    type :: climate_state
        real(real64), allocatable :: temp(:,:)      ! Temperatura superficial [K]
        real(real64), allocatable :: co2_conc(:,:)  ! Concentração de CO2 [ppm]
        real(real64), allocatable :: albedo(:,:)    ! Albedo local
        real(real64), allocatable :: u(:,:), v(:,:) ! Ventos superficiais [m/s]
        real(real64) :: global_avg_temp
    end type climate_state

end module mod_climate_phys
