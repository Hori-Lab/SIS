subroutine force_tweezers(irep, forces)

   use const, only : PREC
   use var_potential, only : ntwz_DCF, twz_DCF_pairs, twz_DCF_forces, &
                             ntwz_FR, twz_FR_pairs, twz_FR_k, twz_FR_init, twz_FR_velo
   use var_state, only : xyz, istep
   use var_top, only : nmp
   use pbc, only : pbc_vec_d

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3,nmp)

   integer :: ipair, imp1, imp2
   real(PREC) :: v(3)

   do ipair = 1, ntwz_DCF

      imp1 = twz_DCF_pairs(1, ipair)
      imp2 = twz_DCF_pairs(2, ipair)

      forces(:, imp1) = forces(:, imp1) + twz_DCF_forces(:, ipair, irep)
      forces(:, imp2) = forces(:, imp2) - twz_DCF_forces(:, ipair, irep)
   enddo

   do ipair = 1, ntwz_FR
      imp1 = twz_FR_pairs(1, ipair)
      imp2 = twz_FR_pairs(2, ipair)

      v(:) = twz_FR_init(:, 1, ipair) + twz_FR_velo(:, 1, ipair) * istep
      v(:) = pbc_vec_d(v(:), xyz(:, imp1, irep))
      forces(:, imp1) = forces(:, imp1) + twz_FR_k(1, ipair) * v(:)

      v(:) = twz_FR_init(:, 2, ipair) + twz_FR_velo(:, 2, ipair) * istep
      v(:) = pbc_vec_d(v(:), xyz(:, imp2, irep))
      forces(:, imp2) = forces(:, imp2) + twz_FR_k(2, ipair) * v(:)
   enddo

end subroutine force_tweezers
