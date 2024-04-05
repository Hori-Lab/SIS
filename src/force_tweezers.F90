subroutine force_tweezers(irep, forces)

   use const, only : PREC
   use const_phys, only : JOUL2KCAL_MOL
   use var_io, only : flg_out_twz, hdl_twz
   use var_potential, only : ntwz_DCF, twz_DCF_pairs, twz_DCF_forces, &
                             ntwz_FR, twz_FR_pairs, twz_FR_k, twz_FR_init, twz_FR_velo
   use var_state, only : xyz, istep, nstep_save
   use var_top, only : nmp
   use pbc, only : pbc_vec_d

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3,nmp)

   integer :: ipair, imp1, imp2
   real(PREC) :: d, p1(3), p2(3), v1(3), v2(3)
   logical :: flg_step_save

   flg_step_save = .False.
   if (flg_out_twz) then
      if (mod(istep, nstep_save) == 0) flg_step_save = .True.
   endif

   !! Dual-trap Constant Force
   do ipair = 1, ntwz_DCF

      imp1 = twz_DCF_pairs(1, ipair)
      imp2 = twz_DCF_pairs(2, ipair)

      forces(:, imp1) = forces(:, imp1) + twz_DCF_forces(:, ipair, irep)
      forces(:, imp2) = forces(:, imp2) - twz_DCF_forces(:, ipair, irep)
   enddo

   !! Force Ramp
   do ipair = 1, ntwz_FR
      imp1 = twz_FR_pairs(1, ipair)
      imp2 = twz_FR_pairs(2, ipair)

      ! position of the trap
      p1(:) = twz_FR_init(:, 1, ipair) + twz_FR_velo(:, 1, ipair) * istep
      ! distance between the trap and the particle
      v1(:) = pbc_vec_d(p1(:), xyz(:, imp1, irep))
      ! force
      forces(:, imp1) = forces(:, imp1) + twz_FR_k(1, ipair) * v1(:)

      ! position of the trap
      p2(:) = twz_FR_init(:, 2, ipair) + twz_FR_velo(:, 2, ipair) * istep
      ! distance between the trap and the particle
      v2(:) = pbc_vec_d(p2(:), xyz(:, imp2, irep))
      ! force
      forces(:, imp2) = forces(:, imp2) + twz_FR_k(2, ipair) * v2(:)

      if (flg_step_save) then
         d = norm2(xyz(:, imp2, irep) - xyz(:, imp1, irep))
         write(hdl_twz, '(i13,x,i3,6(x,g13.8))') istep, ipair, norm2(p2-p1), d, norm2(v1), norm2(v2), &
                  twz_FR_k(1, ipair) * norm2(v1) / (JOUL2KCAL_MOL*1.0e-22), &
                  twz_FR_k(2, ipair) * norm2(v2) / (JOUL2KCAL_MOL*1.0e-22) ! output in pN
      endif
   enddo

end subroutine force_tweezers
