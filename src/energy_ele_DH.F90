subroutine energy_ele_DH(irep, Eele)

   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz, lambdaD, diele_dTcoef, temp_independent
   use var_potential, only : ele_mp, nele, ele_cutoff, ele_coef

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(out) :: Eele

   integer :: iele
   real(PREC) :: dist, rk
   real(PREC) :: e_ele, rcdist

   rcdist = 1.0_PREC / lambdaD(irep)
   e_ele = 0.0_PREC

   !$omp parallel do private(dist) reduction(+:e_ele)
   do iele = 1, nele(irep)

      dist = norm2(pbc_vec_d(xyz(:, ele_mp(1, iele, irep), irep), &
                             xyz(:, ele_mp(2, iele, irep), irep)))

      if (dist > ele_cutoff(irep)) cycle
        
      if (temp_independent == 0) then
         e_ele = e_ele + ele_coef(irep) /dist*exp(-dist*rcdist)
      else
         rk = dist * rcdist
         e_ele = e_ele + ele_coef(irep)/dist*exp(-rk) &
                * (-(1.0_PREC + 0.5_PREC*rk) * diele_dTcoef(irep))

      endif

   end do
   !$omp end parallel do

   Eele = e_ele

end subroutine energy_ele_DH
