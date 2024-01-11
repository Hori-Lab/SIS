subroutine force_rg(irep, forces)

   use const, only : PREC
   use var_state, only : xyz, rg
   use var_top, only : nmp
   use var_potential, only : bias_rg_k, bias_rg_0

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: imp
   real(PREC) :: s2, v(3), cm(3), factor

   ! Center of mass
   cm(:) = sum(xyz(:,:,irep), 2) / real(nmp, kind=PREC)

   ! Rg
   s2 = 0.0_PREC
   do imp = 1, nmp
      v(:) = xyz(:, imp, irep) - cm(:)
      s2 = s2 + dot_product(v,v)
   enddo

   rg(irep) = sqrt(s2 / real(nmp, kind=PREC))

   factor = - bias_rg_k * (rg(irep) - bias_rg_0) / (rg(irep) * nmp)

   do imp = 1, nmp
      forces(:, imp) = forces(:, imp) + factor * (xyz(:, imp, irep) - cm(:))
   enddo

end subroutine force_rg
