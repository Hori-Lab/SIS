subroutine energy_restraint(irep, Erest)

   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : flg_rest_sigb, nrest_sigb, &
                             rest_sigb_id, rest_sigb_rcut, rest_sigb_para

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Erest

   integer :: imp, jmp, irest
   real(PREC) :: r, eps, d, s

   if (flg_rest_sigb) then

      do irest = 1, nrest_sigb
         imp = rest_sigb_id(1, irest)
         jmp = rest_sigb_id(2, irest)

         r = norm2(pbc_vec_d(xyz(:, imp, irep), &
                             xyz(:, jmp, irep)))

         if (r >= rest_sigb_rcut(irest)) cycle

         eps = rest_sigb_para(1, irest)
         d   = rest_sigb_para(2, irest)
         s   = rest_sigb_para(3, irest)

         Erest = Erest - 0.5 * eps * (1.0_PREC - tanh((r - d) / s))
      enddo

   endif

endsubroutine energy_restraint
