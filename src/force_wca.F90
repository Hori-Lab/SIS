subroutine force_wca(forces)

   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : wca_sigma, wca_eps, wca_mp, nwca
   use var_top, only : nmp

   implicit none

   real(PREC), intent(inout) :: forces(3,nmp)
  
   integer :: iwca, imp1, imp2
   real(PREC) :: v12(3), d2, for(3)
   real(PREC) :: wca_sigma_2, coef

   wca_sigma_2 = wca_sigma**2
   coef = 12.0e0_PREC * wca_eps / wca_sigma_2

   !$omp do private(imp1,imp2,v12,d2,for)
   do iwca = 1, nwca

      imp1 = wca_mp(1, iwca)
      imp2 = wca_mp(2, iwca)

      v12(:) = pbc_vec_d(xyz(:,imp1), xyz(:,imp2))
      
      d2 = dot_product(v12, v12)

      if (d2 >= wca_sigma_2) cycle

      d2 = wca_sigma_2 / d2

      for(:) = coef * (d2**7 - d2**4) * v12(:)

      forces(:, imp1) = forces(:, imp1) + for(:)
      forces(:, imp2) = forces(:, imp2) - for(:)
   enddo
   !$omp end do nowait

end subroutine force_wca
