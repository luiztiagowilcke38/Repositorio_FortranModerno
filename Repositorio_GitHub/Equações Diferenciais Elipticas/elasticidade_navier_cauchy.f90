!=======================================================================================
!  MODULO: mod_elasticidade_navier_cauchy
!  OBJETIVO: Solucionador de alto desempenho para o sistema de Navier-Cauchy.
!  AUTOR: Luiz Tiago Wilcke
!=======================================================================================

module mod_elasticidade_tipos
    use, intrinsic :: iso_fortran_env
    implicit none
    type :: PropriedadesMaterial
        real(real64) :: E, nu, lambda, mu
    contains
        procedure :: calcular_lame
    end type PropriedadesMaterial
    type :: CampoElasticidade
        integer :: nx, ny
        real(real64) :: dx, dy
        real(real64), allocatable :: u(:,:), v(:,:), fx(:,:), fy(:,:)
        real(real64), allocatable :: sigma_xx(:,:), sigma_yy(:,:), sigma_xy(:,:)
    end type CampoElasticidade
contains
    subroutine calcular_lame(this)
        class(PropriedadesMaterial), intent(inout) :: this
        this%mu = this%E / (2.0_8 * (1.0_8 + this%nu))
        this%lambda = (this%E * this%nu) / ((1.0_8 + this%nu) * (1.0_8 - 2.0_8 * this%nu))
    end subroutine calcular_lame
end module mod_elasticidade_tipos

module mod_elasticidade_solucionador
    use mod_elasticidade_tipos
    implicit none
contains
    subroutine solver_navier_cauchy(campo, mat, tolerancia, max_it)
        type(CampoElasticidade), intent(inout) :: campo
        type(PropriedadesMaterial), intent(in) :: mat
        real(real64), intent(in) :: tolerancia
        integer, intent(in) :: max_it
        real(real64) :: dx2inv, dy2inv, dxdy4inv, erro_local, c1, c2, c3, denom_u, denom_v, u_ant, v_ant, du_cruz, dv_cruz
        integer :: i, j, iter
        dx2inv = 1.0_8/(campo%dx**2); dy2inv = 1.0_8/(campo%dy**2); dxdy4inv = 1.0_8/(4.0_8*campo%dx*campo%dy)
        c1 = mat%lambda + 2.0_8*mat%mu; c2 = mat%mu; c3 = mat%lambda + mat%mu
        denom_u = 2.0_8*(c1*dx2inv + c2*dy2inv); denom_v = 2.0_8*(c2*dx2inv + c1*dy2inv)
        do iter = 1, max_it
            erro_local = 0.0_8
            !$OMP PARALLEL DO PRIVATE(i, j, u_ant, v_ant, du_cruz, dv_cruz) REDUCTION(max:erro_local)
            do j = 2, campo%ny - 1
                do i = 2, campo%nx - 1
                    u_ant = campo%u(i, j); v_ant = campo%v(i, j)
                    dv_cruz = (campo%v(i+1,j+1)-campo%v(i-1,j+1)-campo%v(i+1,j-1)+campo%v(i-1,j-1))*dxdy4inv
                    du_cruz = (campo%u(i+1,j+1)-campo%u(i-1,j+1)-campo%u(i+1,j-1)+campo%u(i-1,j-1))*dxdy4inv
                    campo%u(i,j) = (campo%fx(i,j) + c1*(campo%u(i+1,j)+campo%u(i-1,j))*dx2inv + c2*(campo%u(i,j+1)+campo%u(i,j-1))*dy2inv + c3*dv_cruz)/denom_u
                    campo%v(i,j) = (campo%fy(i,j) + c2*(campo%v(i+1,j)+campo%v(i-1,j))*dx2inv + c1*(campo%v(i,j+1)+campo%v(i,j-1))*dy2inv + c3*du_cruz)/denom_v
                    erro_local = max(erro_local, abs(campo%u(i,j)-u_ant), abs(campo%v(i,j)-v_ant))
                end do
            end do
            if (erro_local < tolerancia) exit
        end do
    end subroutine solver_navier_cauchy
end module mod_elasticidade_solucionador

program elasticidade_final
    use mod_elasticidade_solucionador
    implicit none
    type(CampoElasticidade) :: viga
    type(PropriedadesMaterial) :: mat
    mat%E = 210e9; mat%nu = 0.3; call mat%calcular_lame()
    viga%nx = 100; viga%ny = 100; viga%dx = 0.02; viga%dy = 0.005
    allocate(viga%u(100,100), viga%v(100,100), viga%fx(100,100), viga%fy(100,100))
    viga%u=0; viga%v=0; viga%fx=0; viga%fy=-77000.0
    call solver_navier_cauchy(viga, mat, 1e-9_8, 20000)
    print *, "Autor: Luiz Tiago Wilcke. Deflexao:", viga%v(50,50)
end program elasticidade_final
