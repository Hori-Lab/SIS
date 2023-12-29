subroutine energy_pull(irep, Epull)

   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : npull_CF, pull_CF_pairs, pull_CF_forces

   implicit none
  
   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Epull

   integer :: ipair, imp1, imp2
   real(PREC) :: v(3)
   real(PREC) :: e_pull

   e_pull = 0.0e0_PREC

   !!$omp parallel do private(imp1, imp2, v) reduction(-:e_pull)
   !!! There will be only a few (most likely only one) pair, so not use OpenMP.
   do ipair = 1, npull_CF
      imp1 = pull_CF_pairs(1, ipair)
      imp2 = pull_CF_pairs(2, ipair)

      v = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))

      e_pull = e_pull - dot_product(v, pull_CF_forces(1:3,ipair))
   enddo
   !!$omp end parallel do

   Epull = e_pull

end subroutine energy_pull
