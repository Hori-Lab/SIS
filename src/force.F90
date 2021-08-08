subroutine force()

   use const
   use const_idx, only : ENE
   use var_state, only : forces

   implicit none

   forces(:,:) = 0.0e0_PREC

   call force_bond()
   call force_angl()
   call force_bp()

end subroutine force
