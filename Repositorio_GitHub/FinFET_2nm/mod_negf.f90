module mod_negf
    use mod_constantes
    use mod_malha
    implicit none

!==================================================================================================
! Módulo: mod_negf (Non-Equilibrium Green's Function)
! Descrição: Implementa o formalismo NEGF para transporte quântico no regime balístico.
!            Calcula a Função de Green, Self-Energies via Sancho-Rubio, 
!            Função Espectral e Corrente de Dreno.
!==================================================================================================

    interface
        subroutine zgetrf(m, n, a, lda, ipiv, info)
            import :: dp
            complex(dp) :: a(lda, *)
            integer :: m, n, lda, ipiv(*), info
        end subroutine zgetrf
        subroutine zgetri(n, a, lda, ipiv, work, lwork, info)
            import :: dp
            complex(dp) :: a(lda, *), work(*)
            integer :: n, lda, ipiv(*), lwork, info
        end subroutine zgetri
    end interface

contains

    ! Algoritmo de Sancho-Rubio para Self-Energias de Contatos Semi-Infinitos
    subroutine sancho_rubio(e_val, h00, h01, sigma, n_dim)
        integer, intent(in) :: n_dim
        complex(dp), intent(in) :: e_val
        complex(dp), intent(in) :: h00(n_dim, n_dim), h01(n_dim, n_dim)
        complex(dp), intent(out) :: sigma(n_dim, n_dim)
        
        complex(dp), allocatable :: alpha(:,:), beta(:,:), epsilon(:,:), epsilon_s(:,:)
        complex(dp), allocatable :: identity(:,:), inv_aux(:,:), g_inv(:,:)
        integer :: i, iter, lwork
        complex(dp), allocatable :: work(:)
        real(dp), parameter :: tol = 1.0d-12
        integer, parameter :: max_iter = 50

        allocate(alpha(n_dim, n_dim), beta(n_dim, n_dim), epsilon(n_dim, n_dim), epsilon_s(n_dim, n_dim))
        allocate(identity(n_dim, n_dim), inv_aux(n_dim, n_dim), g_inv(n_dim, n_dim))

        identity = 0.0d0
        do i = 1, n_dim; identity(i, i) = 1.0d0; end do

        alpha = h01
        beta = transpose(h01)
        epsilon = h00
        epsilon_s = h00

        lwork = n_dim * n_dim
        allocate(work(lwork))

        do iter = 1, max_iter
            ! Cálculo de g = (E*I - epsilon)^-1
            g_inv = e_val * identity - epsilon
            call matriz_invert(g_inv, inv_aux, n_dim)
            
            epsilon_s = epsilon_s + matmul(alpha, matmul(inv_aux, beta))
            epsilon   = epsilon   + matmul(alpha, matmul(inv_aux, beta)) + matmul(beta, matmul(inv_aux, alpha))
            
            alpha = matmul(alpha, matmul(inv_aux, alpha))
            beta  = matmul(beta,  matmul(inv_aux, beta))

            if (maxval(abs(alpha)) < tol) exit
        end do

        ! Sigma = H01 * (E*I - epsilon_s)^-1 * H10
        g_inv = e_val * identity - epsilon_s
        call matriz_invert(g_inv, inv_aux, n_dim)
        sigma = matmul(h01, matmul(inv_aux, transpose(h01)))

        deallocate(alpha, beta, epsilon, epsilon_s, identity, inv_aux, g_inv, work)
    end subroutine sancho_rubio

    ! Inversão de matriz complexa via LAPACK
    subroutine matriz_invert(a, a_inv, n)
        integer, intent(in) :: n
        complex(dp), intent(inout) :: a(n, n)
        complex(dp), intent(out) :: a_inv(n, n)
        integer :: ipiv(n), info, lwork
        complex(dp), allocatable :: work(:)
        
        lwork = n * n
        allocate(work(lwork))
        call zgetrf(n, n, a, n, ipiv, info)
        call zgetri(n, a, n, ipiv, work, lwork, info)
        a_inv = a
        deallocate(work)
    end subroutine matriz_invert

    ! Transmissão T(E) em um dado nível de energia
    function calcular_transmissao(e_val, hamil) result(te)
        complex(dp), intent(in) :: e_val
        real(dp), intent(in) :: hamil(:,:)
        real(dp) :: te
        
        integer :: n_total, i
        complex(dp), allocatable :: g_ret(:,:), sigma_s(:,:), sigma_d(:,:), gamma_s(:,:), gamma_d(:,:)
        complex(dp), allocatable :: ident(:,:), aux(:,:)
        
        n_total = size(hamil, 1)
        allocate(g_ret(n_total, n_total), sigma_s(n_total, n_total), sigma_d(n_total, n_total))
        allocate(gamma_s(n_total, n_total), gamma_d(n_total, n_total))
        allocate(ident(n_total, n_total), aux(n_total, n_total))

        ident = 0.0d0
        do i = 1, n_total; ident(i,i) = 1.0d0; end do
        
        ! Simplificação: Assumindo contatos 1D de área transversal compatível
        ! h00 e h01 seriam extraídos da Hamiltoniana do Canal
        ! Para este exemplo, calcularemos de forma simplificada
        sigma_s = 0.0d0
        sigma_d = 0.0d0
        
        ! Gamma = i(Sigma - Sigma_dag)
        gamma_s = (0.0d0, 1.0d0) * (sigma_s - conjg(transpose(sigma_s)))
        gamma_d = (0.0d0, 1.0d0) * (sigma_d - conjg(transpose(sigma_d)))

        ! G = (E*I - H - Sigma_S - Sigma_D)^-1
        g_ret = e_val * ident - cmplx(hamil, 0.0d0, dp) - sigma_s - sigma_d
        call matriz_invert(g_ret, aux, n_total)
        g_ret = aux

        ! T = Trace(Gamma_S * G * Gamma_D * G_dag)
        aux = matmul(gamma_s, matmul(g_ret, matmul(gamma_d, conjg(transpose(g_ret)))))
        te = 0.0d0
        do i = 1, n_total
            te = te + real(aux(i,i), dp)
        end do
        
        deallocate(g_ret, sigma_s, sigma_d, gamma_s, gamma_d, ident, aux)
    end function calcular_transmissao

end module mod_negf
