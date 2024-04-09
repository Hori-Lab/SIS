subroutine init_tweezers

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const, only : PREC
   use const_idx, only : REPT, RSTBLK
   use var_replica, only : flg_repvar, nrep_proc, rep2val, irep2grep
   use var_potential, only : ntwz_DCF, twz_DCF_direction, twz_DCF_forces, &
                             ntwz_FR, twz_FR_pairs, twz_FR_init, twz_FR_speed, twz_FR_velo
   use var_state, only : xyz, restarted

   implicit none

   integer :: irep, ipair
   integer :: imp1, imp2
   integer :: rst_status
   real(PREC) :: v(3)

   print '(a)', 'Set tweezers parameters'

   if (ntwz_DCF > 0) then
      allocate(twz_DCF_forces(3, ntwz_DCF, nrep_proc))
      twz_DCF_forces(:,:,:) = 0.0_PREC

      if (flg_repvar(REPT%TWZDCF)) then

         ! Currently only ntwz_DCF = 1 is allowed for REMD.
         if (ntwz_DCF /= 1) then
            print '(a)', 'Error: ntwz_DCF must be 1 in f-REMD.'
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
   endif

   if (ntwz_FR > 0) then

      print '(a)', '## Force Ramp mode'

      allocate(twz_FR_init(3, 2, ntwz_FR))
      allocate(twz_FR_velo(3, 2, ntwz_FR))

      if (restarted) then
         call read_rst(RSTBLK%TWZ, rst_status)
      endif

      do ipair = 1, ntwz_FR
         imp1 = twz_FR_pairs(1, ipair)
         imp2 = twz_FR_pairs(2, ipair)

         if (.not. restarted) then
            twz_FR_init(:, 1, ipair) = xyz(:, imp1, 1)
            twz_FR_init(:, 2, ipair) = xyz(:, imp2, 1)
         endif

         v(:) = twz_FR_init(:, 2, ipair) - twz_FR_init(:, 1, ipair)
         twz_FR_velo(:, 1, ipair) = -v(:) / norm2(v) * twz_FR_speed(1, ipair)
         twz_FR_velo(:, 2, ipair) = v(:) / norm2(v) * twz_FR_speed(2, ipair)

         print '(a,x,i5,i5,i5)', '## pair, imp1, imp2: ', ipair, imp1, imp2
         print '(a,x,3(f14.3))', '## initial position 2: ', twz_FR_init(:, 2, ipair)
         print '(a,x,3(f14.3))', '## initial position 1: ', twz_FR_init(:, 1, ipair)
         print '(a,x,g14.8)', '## speed 1: ', twz_FR_speed(1, ipair)
         print '(a,x,g14.8)', '## speed 2: ', twz_FR_speed(2, ipair)
         print '(a,3(x,g14.8))', '## velocity 1: ', twz_FR_velo(:, 1, ipair)
         print '(a,3(x,g14.8))', '## velocity 2: ', twz_FR_velo(:, 2, ipair)
      enddo

      print *
      flush(output_unit)
   endif

endsubroutine init_tweezers
