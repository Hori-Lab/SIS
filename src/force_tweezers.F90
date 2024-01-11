subroutine force_tweezers(irep, forces)

   use const, only : PREC
   use var_potential, only : ntwz_DCF, twz_DCF_pairs, twz_DCF_forces
   use var_top, only : nmp

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3,nmp)

   integer :: ipair, imp1, imp2

   do ipair = 1, ntwz_DCF

      imp1 = twz_DCF_pairs(1, ipair)
      imp2 = twz_DCF_pairs(2, ipair)

      forces(:, imp1) = forces(:, imp1) + twz_DCF_forces(:, ipair)
      forces(:, imp2) = forces(:, imp2) - twz_DCF_forces(:, ipair)
   enddo

end subroutine force_tweezers
