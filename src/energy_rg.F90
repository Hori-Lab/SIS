subroutine energy_rg(irep, Erg)

   use const, only : PREC
   use var_top, only : nmp
   use var_state, only : xyz
   use var_potential, only : bias_rg_k, bias_rg_0

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Erg

   integer :: imp
   real(PREC) :: rg, s2, v(3), cm(3)

   ! Center of mass
   cm(:) = sum(xyz(:,:,irep), 2) / real(nmp, kind=PREC)

   s2 = 0.0_PREC
   do imp = 1, nmp
      v(:) = xyz(:,imp,irep) - cm(:)
      s2 = s2 + dot_product(v,v)
   enddo

   rg = sqrt(s2 / real(nmp, kind=PREC))

   Erg = 0.5_PREC * bias_rg_k * (rg - bias_rg_0)**2

end subroutine energy_rg
