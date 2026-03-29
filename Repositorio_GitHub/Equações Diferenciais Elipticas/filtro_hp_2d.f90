!=======================================================================================
!  MODULO: mod_filtro_hp_2d
!  OBJETIVO: Extração de tendências espaciais (Econometria).
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_hp_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: ConfigHP
        real(real64) :: lambda, dx, dy; integer :: nx, ny
    end type ConfigHP
    type :: MemoriaHP
        real(real64), allocatable :: u(:,:), f(:,:), v(:,:)
    end type MemoriaHP
end module mod_hp_tipos

module mod_hp_solucionador
    use mod_hp_tipos
    implicit none
contains
    subroutine aplicar_filtro_hp(cfg, mem, tol, max_it)
        type(ConfigHP), intent(in) :: cfg; type(MemoriaHP), intent(inout) :: mem
        real(real64), intent(in) :: tol; integer, intent(in) :: max_it
        real(real64) :: h2, erro, old_u, old_v; integer :: i, j, it
        h2 = cfg%dx * cfg%dy; mem%u = mem%f; mem%v = 0.0
        do it = 1, max_it
            erro = 0.0
            !$OMP PARALLEL DO PRIVATE(i, j, old_u, old_v) REDUCTION(max:erro)
            do j = 2, cfg%ny-1; do i = 2, cfg%nx-1
                old_v = mem%v(i,j); mem%v(i,j) = 0.25*(mem%v(i+1,j)+mem%v(i-1,j)+mem%v(i,j+1)+mem%v(i,j-1)-h2*(mem%f(i,j)-mem%u(i,j))/cfg%lambda)
                old_u = mem%u(i,j); mem%u(i,j) = 0.25*(mem%u(i+1,j)+mem%u(i-1,j)+mem%u(i,j+1)+mem%u(i,j-1)-h2*mem%v(i,j))
                erro = max(erro, abs(mem%u(i,j)-old_u))
            end do; end do
            if (erro < tol) exit
        end do
    end subroutine aplicar_filtro_hp
end module mod_hp_solucionador

program hp_final
    use mod_hp_solucionador
    implicit none
    type(ConfigHP) :: cfg; type(MemoriaHP) :: mem; integer :: i,j; real(8) :: r
    cfg%nx=100; cfg%ny=100; cfg%dx=1.0; cfg%dy=1.0; cfg%lambda=1600.0
    allocate(mem%u(100,100), mem%f(100,100), mem%v(100,100))
    do j=1,100; do i=1,100; call random_number(r); mem%f(i,j)=0.01*(i+j)+0.5*sin(i/10.0)+(r-0.5)*0.2; end do; end do
    call aplicar_filtro_hp(cfg, mem, 1e-10, 50000)
    print *, "Autor: Luiz Tiago Wilcke. Tendencia Central:", mem%u(50,50)
end program hp_final
