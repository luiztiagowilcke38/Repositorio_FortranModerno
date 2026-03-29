!=======================================================================================
!  MODULO: mod_ginzburg_landau_dinamica
!  OBJETIVO: Solucionador para a Equacao de TDGL (Dinamica de vortices).
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_tdgl_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: ParametrosSuper
        real(real64) :: alfa, beta, gama, dt, dx
    end type ParametrosSuper
    type :: SistemaTDGL
        integer :: nx, ny
        complex(real64), allocatable :: psi(:,:)
    end type SistemaTDGL
end module mod_tdgl_tipos

module mod_tdgl_solucionador
    use mod_tdgl_tipos
    implicit none
contains
    subroutine evoluir_tdgl(sis, par, passos)
        type(SistemaTDGL), intent(inout) :: sis
        type(ParametrosSuper), intent(in) :: par
        integer, intent(in) :: passos
        integer :: t, i, j
        complex(real64), allocatable :: dpsi(:,:)
        complex(real64) :: lap, termo_nl
        real(real64) :: dx2inv
        allocate(dpsi(sis%nx, sis%ny)); dx2inv = 1.0_8/(par%dx**2)
        do t = 1, passos
            !$OMP PARALLEL DO PRIVATE(i, j, lap, termo_nl)
            do j = 2, sis%ny - 1; do i = 2, sis%nx - 1
                lap = (sis%psi(i+1,j)+sis%psi(i-1,j)+sis%psi(i,j+1)+sis%psi(i,j-1)-4*sis%psi(i,j))*dx2inv
                termo_nl = par%alfa*sis%psi(i,j) + par%beta*abs(sis%psi(i,j))**2*sis%psi(i,j)
                dpsi(i, j) = (1.0_8/par%gama)*(lap - termo_nl)
            end do; end do
            sis%psi = sis%psi + par%dt * dpsi
            sis%psi(1,:)=sis%psi(2,:); sis%psi(sis%nx,:)=sis%psi(sis%nx-1,:)
            sis%psi(:,1)=sis%psi(:,2); sis%psi(:,sis%ny)=sis%psi(:,sis%ny-1)
        end do
    end subroutine evoluir_tdgl
end module mod_tdgl_solucionador

program tdgl_final
    use mod_tdgl_solucionador
    implicit none
    type(SistemaTDGL) :: amos
    type(ParametrosSuper) :: par
    par%alfa=-1.0; par%beta=1.0; par%gama=1.0; par%dt=0.01; par%dx=0.5
    amos%nx=128; amos%ny=128; allocate(amos%psi(128,128)); amos%psi=(0.1, 0.0)
    call evoluir_tdgl(amos, par, 10000)
    print *, "Autor: Luiz Tiago Wilcke. Densidade:", abs(amos%psi(64,64))**2
end program tdgl_final
