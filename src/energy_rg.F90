subroutine energy_rg(irep, Erg)

   use const, only : PREC
   use const_idx, only : POTT
   use var_top, only : nmp
   use var_state, only : xyz, rg
   use var_potential, only : bias_rg_k, bias_rg_0, bias_rg_pott

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Erg

   integer :: imp
   real(PREC) :: s2, v(3), cm(3)

   ! Center of mass
   cm(:) = sum(xyz(:,:,irep), 2) / real(nmp, kind=PREC)

   s2 = 0.0_PREC
   do imp = 1, nmp
      v(:) = xyz(:,imp,irep) - cm(:)
      s2 = s2 + dot_product(v,v)
   enddo

   rg(irep) = sqrt(s2 / real(nmp, kind=PREC))

   if (bias_rg_pott == POTT%FLATBOTTOMED) then
      if (rg(irep) <= bias_rg_0) then
         Erg = 0.0_PREC
         return
      endif
   endif

   Erg = 0.5_PREC * bias_rg_k * (rg(irep) - bias_rg_0)**2

end subroutine energy_rg
