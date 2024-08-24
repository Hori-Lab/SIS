subroutine force_restraint(irep, forces)

   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_top, only : nmp
   use var_state, only : xyz
   use var_potential, only : flg_rest_sigb, nrest_sigb, &
                             rest_sigb_id, rest_sigb_rcut, rest_sigb_para

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: imp, jmp, irest
   real(PREC) :: r, eps, d, s, coef, v(3)

   if (flg_rest_sigb) then
      do irest = 1, nrest_sigb
         imp = rest_sigb_id(1, irest)
         jmp = rest_sigb_id(2, irest)

         v(:) = pbc_vec_d(xyz(:, imp, irep), xyz(:, jmp, irep))
         r = norm2(v)

         if (r >= rest_sigb_rcut(irest)) cycle

         eps = rest_sigb_para(1, irest)
         d   = rest_sigb_para(2, irest)
         s   = rest_sigb_para(3, irest)

         coef = 0.5_PREC * eps / s * (1.0_PREC - tanh((r-d)/s)**2) / r

         forces(:, imp) = forces(:, imp) - coef * v(:)
         !forces(:, jmp) = forces(:, jmp) + coef * v(:)
      enddo

   endif

end subroutine force_restraint
