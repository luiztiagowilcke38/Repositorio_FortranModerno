module mod_poisson
    use mod_constantes
    use mod_malha
    implicit none

!==================================================================================================
! Módulo: mod_poisson
! Descrição: Resolve a equação de Poisson 3D: Grad(eps * Grad(phi)) = -rho.
!            Utiliza o método de Diferenças Finitas e o algoritmo de Gauss-Seidel 
!            com relaxação (SOR) para eficiência.
!==================================================================================================

contains

    subroutine resolver_poisson(m, phi, rho, convergido)
        type(t_malha), intent(in) :: m
        real(dp), intent(inout) :: phi(m%nx, m%ny, m%nz)
        real(dp), intent(in) :: rho(m%nx, m%ny, m%nz)
        logical, intent(out) :: convergido
        
        real(dp) :: phi_new, residuo, max_residuo
        real(dp), parameter :: omega = 1.85d0  ! Fator de sobre-relaxação
        real(dp), parameter :: tol = 1.0d-6
        integer :: i, j, k, iter
        integer, parameter :: max_iter = 5000

        max_residuo = 0.0d0
        convergido = .false.

        do iter = 1, max_iter
            max_residuo = 0.0d0
            do i = 2, m%nx - 1
                do j = 2, m%ny - 1
                    do k = 2, m%nz - 1
                        ! Discretização 3D da equação de Poisson
                        phi_new = ( (phi(i+1,j,k) + phi(i-1,j,k)) / m%dx**2 + &
                                    (phi(i,j+1,k) + phi(i,j-1,k)) / m%dy**2 + &
                                    (phi(i,j,k+1) + phi(i,j,k-1)) / m%dz**2 + &
                                    rho(i,j,k) / eps_si ) / &
                                  ( 2.0d0 / m%dx**2 + 2.0d0 / m%dy**2 + 2.0d0 / m%dz**2 )
                        
                        residuo = abs(phi_new - phi(i,j,k))
                        if (residuo > max_residuo) max_residuo = residuo
                        
                        ! Atualização com SOR
                        phi(i,j,k) = phi(i,j,k) + omega * (phi_new - phi(i,j,k))
                    end do
                end do
            end do
            
            if (max_residuo < tol) then
                convergido = .true.
                exit
            end if
        end do
        
        print *, "Poisson: Itéritos = ", iter, " Residuo Final = ", max_residuo
    end subroutine resolver_poisson

    ! Condição de Contorno de Dirichlet (Potencial de Gate e Bias)
    subroutine aplicar_contornos_poisson(m, phi, vg, vd, vs)
        type(t_malha), intent(in) :: m
        real(dp), intent(inout) :: phi(m%nx, m%ny, m%nz)
        real(dp), intent(in) :: vg, vd, vs

        ! Contato Source (i=1)
        phi(1,:,:) = vs
        ! Contato Drain (i=nx)
        phi(m%nx,:,:) = vd
        
        ! Gate (Interface óxido-semicondutor nas faces laterais e superior)
        ! Face Superior (k=nz)
        phi(:, :, m%nz) = vg
        ! Faces Laterais (j=1 e j=ny)
        phi(:, 1, :) = vg
        phi(:, m%ny, :) = vg
        
        ! Base do Fin (Buried Oxide ou Substrato - isolante)
        ! Aqui poderíamos usar Neumann (grad(phi)=0), mas para simplificar:
        phi(:, :, 1) = 0.0d0 ! Referência
    end subroutine aplicar_contornos_poisson

end module mod_poisson
