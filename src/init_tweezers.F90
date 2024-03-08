subroutine init_tweezers

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const, only : PREC
   use const_idx, only : REPT
   use var_replica, only : flg_repvar, nrep_proc, rep2val, irep2grep
   use var_potential, only : ntwz_DCF, twz_DCF_direction, twz_DCF_forces

   implicit none

   integer :: irep

   allocate(twz_DCF_forces(3, ntwz_DCF, nrep_proc))
   twz_DCF_forces(:,:,:) = 0.0_PREC

   if (flg_repvar(REPT%TWZDCF)) then

      ! Currently only ntwz_DCF = 1 is allowed for REMD.
      if (ntwz_DCF /= 1) then
         print '(2a)', 'Error: ntwz_DCF must be 1 in f-REMD.'
         call sis_abort()
      endif

      twz_DCF_direction(1:3, 1) = twz_DCF_direction(1:3, 1) / norm2(twz_DCF_direction(1:3, 1))

      print '(a)', '# Normalised force vector is,'
      print '(a,f8.2,a,f8.2,a,f8.2,a)', '# (', twz_DCF_direction(1,1), ',', &
                     twz_DCF_direction(2,1), ',', twz_DCF_direction(3,1), ')'

      do irep = 1, nrep_proc
         twz_DCF_forces(:, 1, irep) = rep2val(irep2grep(irep), REPT%TWZDCF) * twz_DCF_direction(:, 1)
      enddo

   else
      do irep = 1, nrep_proc
         twz_DCF_forces(:, :, irep) = twz_DCF_direction(:, :)
      enddo

   endif

   print *
   flush(output_unit)

endsubroutine init_tweezers
