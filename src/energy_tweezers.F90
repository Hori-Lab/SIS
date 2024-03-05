subroutine energy_tweezers(irep, Etwz)

   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : ntwz_DCF, twz_DCF_pairs, twz_DCF_forces

   implicit none
  
   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Etwz

   integer :: ipair, imp1, imp2
   real(PREC) :: v(3)
   real(PREC) :: e_twz

   e_twz = 0.0e0_PREC

   !!$omp parallel do private(imp1, imp2, v) reduction(-:e_twz)
   !!! There will be only a few (most likely only one) pair, so not use OpenMP.
   do ipair = 1, ntwz_DCF
      imp1 = twz_DCF_pairs(1, ipair)
      imp2 = twz_DCF_pairs(2, ipair)

      v = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))

      e_twz = e_twz - dot_product(v, twz_DCF_forces(1:3,ipair,irep))
   enddo
   !!$omp end parallel do

   Etwz = e_twz

end subroutine energy_tweezers
