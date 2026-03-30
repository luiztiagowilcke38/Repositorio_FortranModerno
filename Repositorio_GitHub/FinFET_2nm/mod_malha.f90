module mod_malha
    use mod_constantes
    implicit none

!==================================================================================================
! Módulo: mod_malha
! Descrição: Gerencia a discretização espacial (1D/2D) do FinFET. Define o número de pontos,
!            o passo da malha e as coordenadas x, y, z.
!==================================================================================================

    type, public :: t_malha
        integer :: nx, ny, nz
        real(dp) :: dx, dy, dz
        real(dp), allocatable :: x(:), y(:), z(:)
    end type t_malha

contains

    subroutine inicializar_malha(m, nx_in, ny_in, nz_in)
        type(t_malha), intent(out) :: m
        integer, intent(in) :: nx_in, ny_in, nz_in
        integer :: i

        m%nx = nx_in
        m%ny = ny_in
        m%nz = nz_in

        ! Discretização para o canal (comprimento L)
        m%dx = l_canal / real(m%nx - 1, dp)
        allocate(m%x(m%nx))
        do i = 1, m%nx
            m%x(i) = (i - 1) * m%dx
        end do

        ! Discretização transversal (W_fin)
        m%dy = w_fin / real(m%ny - 1, dp)
        allocate(m%y(m%ny))
        do i = 1, m%ny
            m%y(i) = (i - 1) * m%dy
        end do

        ! Discretização vertical (H_fin)
        m%dz = h_fin / real(m%nz - 1, dp)
        allocate(m%z(m%nz))
        do i = 1, m%nz
            m%z(i) = (i - 1) * m%dz
        end do
        
        print *, "--- Malha de Discretização Gerada ---"
        print *, "Pontos (nx, ny, nz): ", m%nx, m%ny, m%nz
        print *, "Passo (dx, dy, dz):  ", m%dx * 1.0d9, " nm", m%dy * 1.0d9, " nm", m%dz * 1.0d9, " nm"
        print *, "-------------------------------------"
    end subroutine inicializar_malha

    subroutine liberar_malha(m)
        type(t_malha), intent(inout) :: m
        if (allocated(m%x)) deallocate(m%x)
        if (allocated(m%y)) deallocate(m%y)
        if (allocated(m%z)) deallocate(m%z)
    end subroutine liberar_malha

end module mod_malha
