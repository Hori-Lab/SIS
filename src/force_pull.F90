subroutine force_pull(irep, forces)

   use const, only : PREC
   use var_potential, only : npull_CF, pull_CF_pairs, pull_CF_forces
   use var_top, only : nmp

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3,nmp)
  
   integer :: ipair, imp1, imp2
   
   do ipair = 1, npull_CF

      imp1 = pull_CF_pairs(1, ipair)
      imp2 = pull_CF_pairs(2, ipair)

      forces(:, imp1) = forces(:, imp1) + pull_CF_forces(:, ipair)
      forces(:, imp2) = forces(:, imp2) - pull_CF_forces(:, ipair)
   enddo

end subroutine force_pull
