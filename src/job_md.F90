subroutine job_md()

   use, intrinsic :: iso_fortran_env, Only : iostat_end, INT64
   use const
   use const_phys, only : KCAL2JOUL, N_AVO, PI, BOLTZ_KCAL_MOL
   use const_idx, only : ENE, SEQT, RSTBLK, BPT, REPT
   use progress, only : progress_init, progress_update, wall_time_sec
   use pbc, only : pbc_box, set_pbc_size, flg_pbc
   use var_top, only : nmp, seq, mass, lmp_mp, ichain_mp
   use var_state, only : restarted, flg_bp_energy, &
                         viscosity_Pas, xyz,  energies, dt, velos, accels, tempK, &
                         nstep, nstep_save, nstep_save_rst, stop_wall_time_sec, fix_com_origin, &
                         nl_margin, Ekinetic, &
                         flg_variable_box, variable_box_step, variable_box_change, &
                         opt_anneal, nanneal, anneal_tempK, anneal_step, &
                         istep, ianneal, istep_anneal_next
   use var_io, only : flg_progress, step_progress, hdl_dcd, hdl_out, cfile_dcd, hdl_rep
   use var_potential, only : stage_sigma, wca_sigma, bp_paras, bp_cutoff_energy, bp_cutoff_dist, &
                             ele_cutoff, flg_stage
   use var_replica, only : nrep_all, nrep_proc, flg_replica, rep2val, irep2grep, rep2lab, &
                           nstep_rep_exchange, nstep_rep_save, nrep_all
   use var_parallel, only : myrank
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   integer :: i, irep, imp, icol
   real(PREC) :: tK
   real(PREC) :: dxyz(3)
   real(PREC) :: xyz_move(3, nmp, nrep_proc)
   type(file_dcd) :: fdcd(nrep_proc)
   real(PREC) :: fric(nmp), radius
   real(PREC) :: v, c1, c2, md_coef(3, nmp), md_coef_rep(1, nmp, nrep_proc)
   real(PREC) :: rnd_bm(3, nmp)
   real(PREC) :: accels_pre(3)
   real(PREC) :: d2, d2max, d2max_2nd
   real(PREC), allocatable :: forces(:, :)
   real(PREC), allocatable :: replica_energies(:, :)
   logical :: flg_stop
   character(len=37) :: out_fmt

   ! Function
   real(PREC) :: rnd_boxmuller

   ! Format in .out file
   if (flg_replica) then
      write(out_fmt, '(a23,i2,a12)') '(i10, 1x, i4, 1x, f6.2,', ENE%ELE+2, '(1x, g13.6))'
   else
      write(out_fmt, '(a15,i2,a12)') '(i10, 1x, f6.2,', ENE%ELE+2, '(1x, g13.6))'
   endif

   allocate(mass(nmp))
   allocate(accels(3, nmp, nrep_proc))
   allocate(velos(3, nmp, nrep_proc))
   allocate(forces(3, nmp))
   allocate(energies(0:ENE%MAX, nrep_proc))
   allocate(Ekinetic(nrep_proc))
   allocate(replica_energies(2, nrep_all))

   ! set PBC box
   if (flg_pbc) then
      do irep = 1, nrep_proc
         fdcd(irep)%flg_unitcell = .True.
         fdcd(irep)%box(:) = pbc_box(:)
      enddo
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Set up constants for the dynamics
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   radius = 10.0
   v = viscosity_Pas * sqrt(1.0e3_PREC / KCAL2JOUL) * N_AVO * 1.0e-20_PREC
   fric(:) = 6.0e0_PREC * PI * v * radius
   print '(a)', 'Langevin dynamics parameters'
   print '(a,f10.3)', 'Stokes radius = ', radius
   print '(a,f10.3)', 'Viscosity = ', v
   print '(a,f10.3)', 'Friction coefficient = ', fric(1)

   !v = 0.5 ! ps^(-1)
   !write(*,*) 'v =', v

   print '(a)', 'mass(A) = 328.212'
   print '(a)', 'mass(G) = 344.212'
   print '(a)', 'mass(C) = 304.182'
   print '(a)', 'mass(U) = 305.164'
   print *
   flush(6)

   ! Coefficients that do not depend on replicas
   do imp = 1, nmp
      ! Mass hard coded
      select case (seq(lmp_mp(imp), ichain_mp(imp)))
         case (SEQT%A)
            mass(imp) = 328.212
         case (SEQT%G)
            mass(imp) = 344.212
         case (SEQT%C)
            mass(imp) = 304.182
         case (SEQT%U)
            mass(imp) = 305.164
         case default
            write(*,*) 'Error: Unknown seq type, ', seq(lmp_mp(imp), ichain_mp(imp)), ' imp=',imp
            call sis_abort()
      endselect
   enddo

   call set_md_coef()


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! set up the initial state
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (restarted) then
      call read_rst(RSTBLK%VELO)
      call read_rst(RSTBLK%ACCEL)

   else

      if (flg_stage) call check_int_stage()

      ! Set up initial velocities by Maxwell–Boltzmann distribution
      do irep = 1, nrep_proc
         if (flg_replica) then
            tK = rep2val(irep2grep(irep), REPT%TEMP)
         else
            tK = tempK
         endif
         do imp = 1, nmp
            c1 = sqrt(tK * BOLTZ_KCAL_MOL / mass(imp))
            do i = 1, 3
               velos(i, imp, irep) = c1 * rnd_boxmuller(irep)
            enddo
         enddo

         ! Set up initial accerelations
         do imp = 1, nmp
            do i = 1, 3
               accels(i, imp, irep) = md_coef_rep(1, imp, irep) * rnd_boxmuller(irep)
            enddo
         end do
      enddo
   endif


   print '(a)', 'Potential and neighbor list cutoffs'
   print '(a,f10.3)', 'bp_cutoff_energy = ', bp_cutoff_energy
   print '(a,f10.3)', 'bp_cutoff_ddist(GC) = ', bp_paras(BPT%GC)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_ddist(AU) = ', bp_paras(BPT%AU)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_ddist(GU) = ', bp_paras(BPT%GU)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_dist(for neighbor list) = ', bp_cutoff_dist
   print '(a,f10.3)', 'wca_sigma = ', wca_sigma
   print '(a,f10.3)', 'stage_sigma = ', stage_sigma
   do irep = 1, nrep_proc
      print '(a,i4,a,f10.3)', 'ele_cutoff(irep=',irep,') = ', ele_cutoff(irep)
   enddo
   print '(a,f10.3)', 'nl_margin = ', nl_margin
   print *

   do irep = 1, nrep_proc
      call neighbor_list(irep)
   enddo
   xyz_move(:,:,:) = 0.0e0_PREC

   if (restarted) then
      call read_rst(RSTBLK%STEP)
   else
      istep = 0_INT64
   endif

   ! Setting up Simulated Annealing
   ianneal = 1
   if (opt_anneal > 0) then

      if (istep >= anneal_step(nanneal)) then
         ianneal = nanneal
         opt_anneal = 0

      else if (restarted) then
         call read_rst(RSTBLK%ANNEAL)

      endif

      tempK = anneal_tempK(ianneal)

      call set_md_coef()
      !call set_bp_map()

      if (ianneal < nanneal) then
         istep_anneal_next = anneal_step(ianneal+1)
      endif
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Initial energies
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do irep = 1, nrep_proc

      flg_bp_energy = .False.
      call energy(irep, energies(0:ENE%MAX, irep))
      call energy_kinetic(irep, Ekinetic(irep))

      ! Open DCD file and write the header
      print '(a,i5,2a)', '# Opening dcd file for irep = ', irep, ', ', trim(cfile_dcd(irep))
      fdcd(irep) = file_dcd(hdl_dcd(irep), cfile_dcd(irep), DCD_OPEN_MODE%WRITE)

      call fdcd(irep)%write_header(nmp)

      ! Open .out file
      if (flg_replica) then
         write(hdl_out(irep), '(a)', advance='no') '#(1)nframe (2)R (3)T   (4)Ekin       (5)Epot       (6)Ebond     '
                                                   !1234567890 1234 123456 1234567890123 1234567890123 1234567890123'
         write(hdl_out(irep), '(a)', advance='no') ' (7)Eangl      (8)Edih      '
                                                   ! 1234567890123 1234567890123'
         write(hdl_out(irep), '(a)', advance='no') ' (9)Ebp        (10)Eexv      (11)Eele'
                                                   ! 1234567890123 1234567890123 1234567890123
         icol = 11
      else
         write(hdl_out(irep), '(a)', advance='no') '#(1)nframe (2)T   (3)Ekin       (4)Epot       (5)Ebond     '
                                                   !1234567890 123456 1234567890123 1234567890123 1234567890123'
         write(hdl_out(irep), '(a)', advance='no') ' (6)Eangl      (7)Edih      '
                                                   ! 1234567890123 1234567890123'
         write(hdl_out(irep), '(a)', advance='no') ' (8)Ebp        (9)Eexv       (10)Eele'
                                                   ! 1234567890123 1234567890123 1234567890123
         icol = 10
      endif
      if (flg_stage) then
         icol = icol + 1
         write(hdl_out(irep), '(a,i2,a)', advance='no') ' (', icol, ')Estage   '
                                                        ! 1   23     4567890123
      endif
      write(hdl_out(irep), *) ''

      if (flg_replica) then
         tK = rep2val(irep2grep(irep), REPT%TEMP)
      else
         tK = tempK
      endif

      ! Output initial structure
      if (restarted) then
         ! Only STDOUT if restarted. No DCD output.
         print '(a)', '##### Energies at the beginning'
         if (flg_replica) then
            print '(a)', '#(1)nframe (2)R (3)T   (4)Ekin       (5)Epot       '
            print '(i10, 1x, i4, 1x, f6.2, 2(1x,g13.6))', istep, irep2grep(irep), tK, Ekinetic(irep), energies(0, irep)
            print '(a)', '(6)Ebond      (7)Eangl      (8)Edih       (9)Ebp        (10)Eexv      (11)Eele      (11)Estage'
            print '(7(1x,g13.6))', (energies(i, irep), i=1, ENE%MAX)
         else
            print '(a)', '#(1)nframe (2)T   (3)Ekin       (4)Epot       '
            print '(i10, 1x, i4, 1x, f6.2, 2(1x,g13.6))', istep, tK, Ekinetic(irep), energies(0, irep)
            print '(a)', '(5)Ebond      (6)Eangl      (7)Edih       (8)Ebp        (9)Eexv       (10)Eele      (11)Estage'
            print '(7(1x,g13.6))', (energies(i, irep), i=1, ENE%MAX)
         endif
         print *

      else
         ! At istep = 0 (not restarted), write both .out and DCD
         if (flg_replica) then
            write(hdl_out(irep), out_fmt, advance='no') istep, irep2grep(irep), tempK, Ekinetic(irep), (energies(i, irep), i=0,ENE%ELE)
         else
            write(hdl_out(irep), out_fmt, advance='no') istep, tempK, Ekinetic(irep), (energies(i, irep), i=0,ENE%ELE)
         endif
         if (flg_stage) then
            write(hdl_out(irep), '(1x, g13.6)', advance='no') energies(ENE%STAGE, irep)
         endif
         write(hdl_out(irep), *) ''

         call fdcd(irep)%write_onestep(nmp, xyz(:,:,irep), fix_com_origin)
      endif
   enddo
   print *
   flush(6)

   if (flg_progress) then
      call progress_init(istep)
   endif

   flg_stop = .False.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Main loop for time integration
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do while (istep < nstep)

      istep = istep + 1

      do irep = 1, nrep_proc

         d2max = 0.0e0_PREC
         d2max_2nd = 0.0e0_PREC
         do imp = 1, nmp
            d2 = dot_product(xyz_move(1:3, imp, irep), xyz_move(1:3, imp, irep))
            if (d2 > d2max) then
               d2max_2nd = d2max
               d2max = d2
            else if (d2 > d2max_2nd) then
               d2max_2nd = d2
            endif
         enddo

         if (sqrt(d2max) + sqrt(d2max_2nd) > nl_margin) then
            call neighbor_list(irep)
            xyz_move(:, :, irep) = 0.0e0_PREC
         endif

         call force(irep, forces(:,:))

         if (maxval(forces(:,:)) > 100.0) then
            print *, 'irep, istep, force ', irep, istep, maxval(forces)
            flush(6)
         endif

         do imp= 1, nmp
            do i = 1, 3
               rnd_bm(i, imp) = rnd_boxmuller(irep)
            enddo
         enddo

         !$omp parallel do private(accels_pre, dxyz)
         do imp = 1, nmp
            !if(fix_mp(imp)) cycle

            !! a = (1 - gamma h / 2m) / (1 + gamma h / 2m)
            !! b = 1 / (1 + gamma h / 2m)

            ! beta(t + h) with the associated coefficient
            accels_pre(1:3) = md_coef_rep(1, imp, irep) * rnd_bm(1:3, imp)
            ! md_coef_rep(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)

            ! v(t + 1/2h) update the half-step velocity
            velos(1:3, imp, irep) =  md_coef(1, imp) * velos(1:3, imp, irep)  &
                             + md_coef(2, imp) * forces(1:3, imp) &
                             + (accels(1:3, imp, irep) + accels_pre(1:3))
            ! md_coef(1) = a
            ! md_coef(2) = sqrt(b) h / m

            ! beta(t) <= beta(t+h) (incluing the coefficient) save for the next iteration
            accels(1:3, imp, irep) = accels_pre(1:3)

            dxyz(1:3) =  md_coef(3, imp) * velos(1:3, imp, irep)
            ! md_coef(3) = sqrt(b) h

            xyz(1:3, imp, irep) = xyz(1:3, imp, irep) + dxyz(1:3)
            xyz_move(1:3, imp, irep) = xyz_move(1:3, imp, irep) + dxyz(1:3)
         end do
         !$omp end parallel do

      enddo ! irep

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Energy calculation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (mod(istep, nstep_save) == 0) then
         do irep = 1, nrep_proc
            if (flg_replica) then
               tK = rep2val(irep2grep(irep), REPT%TEMP)
            else
               tK = tempK
            endif
            call energy(irep, energies(0:ENE%MAX, irep))
            call energy_kinetic(irep, Ekinetic(irep))

            if (flg_replica) then
               write(hdl_out(irep), out_fmt, advance='no') istep, irep2grep(irep), tempK, Ekinetic(irep), (energies(i, irep), i=0,ENE%ELE)
            else
               write(hdl_out(irep), out_fmt, advance='no') istep, tempK, Ekinetic(irep), (energies(i, irep), i=0,ENE%ELE)
            endif
            if (flg_stage) then
               write(hdl_out(irep), '(1x, g13.6)', advance='no') energies(ENE%STAGE, irep)
            endif
            write(hdl_out(irep), *) ''
            call fdcd(irep)%write_onestep(nmp, xyz(:,:,irep), fix_com_origin)
         enddo
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Replica exchange
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flg_replica) then
         if (mod(istep, nstep_rep_exchange) == 0) then

            call replica_exchange(velos, replica_energies, tempk)

            ! Write rep file
            if (myrank == 0 .and. mod(istep, nstep_rep_save) == 0) then
               write(hdl_rep, '(i12,1x)', ADVANCE='no') istep

               do irep = 1, nrep_all - 1
                  write(hdl_rep, '(i5,1x)', ADVANCE='no') rep2lab(irep)
               enddo

               write(hdl_rep, '(i5)') rep2lab(nrep_all)
            endif
         endif
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Update the box size
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flg_variable_box) then
         if (mod(istep, variable_box_step) == 0) then
            call set_pbc_size(pbc_box(:) + variable_box_change(:))
            fdcd(irep)%box(:) = pbc_box(:)

            call neighbor_list(irep)
            xyz_move(:,:,irep) = 0.0e0_PREC

            print '(a,i10,a,f8.3)', 'Box size updated: step = ',istep, ', box size = ', pbc_box(1)
         endif
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Annealing
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (opt_anneal > 0 .and. istep + 1 == istep_anneal_next) then
         ianneal = ianneal + 1
         tempK = anneal_tempK(ianneal)

         call set_md_coef()
         !call set_bp_map()

         if (ianneal < nanneal) then
            istep_anneal_next = anneal_step(ianneal+1)
         else
            opt_anneal = 0
         endif
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Write progress
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (flg_progress) then
         if (mod(istep, step_progress) == 0) then
            call progress_update(istep, nstep)
         endif
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Check wall time
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (stop_wall_time_sec > 0 .and. wall_time_sec() > stop_wall_time_sec) then
         flg_stop = .True.
         print '(a)', 'Wall-clock time limit reached. Stop the job after writing rst file.'
      endif


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! Write restart file (.rst)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (mod(istep, nstep_save_rst) == 0 .or. flg_stop) then
         call write_rst()
      endif

      if (flg_stop) exit

   enddo  ! <--- Main loop for time integration

   do irep = 1, nrep_proc
      call fdcd(irep)%close()
   enddo

contains

   subroutine set_md_coef()

      ! Coefficients that do not depend on replicas
      do imp = 1, nmp
         !fric(imp) = v * mass(imp)
         !if (imp == 1) then
         !   write(*,*) 'fric = ', fric(imp)
         !endif

         !! sqrt(b) = sqrt(1 / (1 + gamma h / 2m))
         c1 = 0.5 * dt * fric(imp) / mass(imp)
         c2 = sqrt(1.0e0_PREC / (1.0_PREC + c1))

         ! md_coef(1) = a
         md_coef(1, imp) = (1.0_PREC - c1) / (1.0_PREC + c1)
         ! md_coef(2) = sqrt(b) h / m
         md_coef(2, imp) = c2 * dt / mass(imp)
         ! md_coef(3) = sqrt(b) h
         md_coef(3, imp) = c2 * dt
      enddo

      ! Coefficients depending on replicas
      do irep = 1, nrep_proc
         do imp = 1, nmp
            c1 = 0.5 * dt * fric(imp) / mass(imp)
            c2 = sqrt(1.0e0_PREC / (1.0_PREC + c1))
            if (flg_replica) then
               tK = rep2val(irep2grep(irep), REPT%TEMP)
            else
               tK = tempK
            endif
            ! md_coef_rep(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)
            md_coef_rep(1, imp, irep) = 0.5_PREC * c2 / mass(imp) * sqrt(2.0_PREC * fric(imp) * BOLTZ_KCAL_MOL * tK * dt)
         enddo
      enddo

   endsubroutine set_md_coef

endsubroutine job_md
