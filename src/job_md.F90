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
   use var_io, only : flg_progress, step_progress, hdl_dcd, hdl_out, cfile_dcd
   use var_potential, only : wca_sigma, bp_paras, bp_cutoff_energy, bp_cutoff_dist
   use var_replica, only : nrep_proc, flg_replica, rep2val, irep2grep
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   integer :: i, irep, imp
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
   logical :: flg_stop

   ! Function
   real(PREC) :: rnd_boxmuller

   allocate(mass(nmp))
   allocate(accels(3, nmp, nrep_proc))
   allocate(velos(3, nmp, nrep_proc))
   allocate(forces(3, nmp))
   allocate(energies(0:ENE%MAX, nrep_proc))
   allocate(Ekinetic(nrep_proc))

   ! set PBC box
   if (flg_pbc) then
      do irep = 1, nrep_proc
         fdcd(irep)%flg_unitcell = .True.
         fdcd(irep)%box(:) = pbc_box(:)
      enddo
   endif

   !! Set up variables for dynamics
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


   !! Set up the initial state
   if (restarted) then
      call read_rst(RSTBLK%VELO)
      call read_rst(RSTBLK%ACCEL)

   else

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
      call set_bp_map()

      if (ianneal < nanneal) then
         istep_anneal_next = anneal_step(ianneal+1)
      endif
   endif

   ! Initial energies
   do irep = 1, nrep_proc

      flg_bp_energy = .False.
      call energy(irep, energies(0:ENE%MAX, irep))
      call energy_kinetic(irep, Ekinetic(irep))

      ! Open DCD file and write the header
      print '(a,i5,2a)', '# Opening dcd file for irep = ', irep, ', ', trim(cfile_dcd(irep))
      fdcd(irep) = file_dcd(hdl_dcd(irep), cfile_dcd(irep), DCD_OPEN_MODE%WRITE)

      call fdcd(irep)%write_header(nmp)

      ! Open .out file
      write(hdl_out(irep), '(a)', advance='no') '#(1)nframe (2)R (3)T   (4)Ekin       (5)Epot       (6)Ebond     '
                                                !1234567890 1234 123456 1234567890123 1234567890123 1234567890123'
      write(hdl_out(irep), '(a)', advance='no') ' (7)Eangl      (8)Edih      '
                                                ! 1234567890123 1234567890123'
      write(hdl_out(irep), '(a)') ' (9)Ebp        (10)Eexv      (11)Eele'
                                  ! 1234567890123 1234567890123 1234567890123

      if (flg_replica) then
         tK = rep2val(irep2grep(irep), REPT%TEMP)
      else
         tK = tempK
      endif

      ! Output initial structure
      if (restarted) then
         ! Only STDOUT if restarted. No DCD output.
         print '(a)', '##### Energies at the beginning'
         print '(a)', '#(1)nframe (2)R (3)T   (4)Ekin       (5)Epot       '
         print '(i10, 1x, i4, 1x, f6.2, 2(1x,g13.6))', istep, irep2grep(irep), tK, Ekinetic, energies(0, irep)
         print '(a)', '(6)Ebond      (7)Eangl      (8)Edih       (9)Ebp        (10)Eexv      (11)Eele'
         print '(6(1x,g13.6))', (energies(i, irep), i=1, ENE%MAX)
         print *

      else
         ! At istep = 0 (not restarted), write both .out and DCD
         write(hdl_out(irep), '(i10, 1x, i4, 1x, f6.2, 8(1x,g13.6))') istep, irep2grep(irep), tK, &
                             Ekinetic(irep), (energies(i, irep), i=0,ENE%MAX)
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
   !!! Time integration
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

      if (mod(istep, nstep_save) == 0) then
         do irep = 1, nrep_proc
            if (flg_replica) then
               tK = rep2val(irep2grep(irep), REPT%TEMP)
            else
               tK = tempK
            endif
            call energy(irep, energies(0:ENE%MAX, irep))
            call energy_kinetic(irep, Ekinetic(irep))
            write(hdl_out(irep), '(i10, 1x, i4, 1x, f6.2, 8(1x,g13.6))') istep, irep2grep(irep), tK, &
                                 Ekinetic(irep), (energies(i, irep), i=0,ENE%MAX)
            call fdcd(irep)%write_onestep(nmp, xyz(:,:,irep), fix_com_origin)
         enddo
      endif

      if (flg_variable_box) then
         if (mod(istep, variable_box_step) == 0) then
            call set_pbc_size(pbc_box(:) + variable_box_change(:))
            fdcd(irep)%box(:) = pbc_box(:)

            call neighbor_list(irep)
            xyz_move(:,:,irep) = 0.0e0_PREC

            print '(a,i10,a,f8.3)', 'Box size updated: step = ',istep, ', box size = ', pbc_box(1)
         endif
      endif

      if (opt_anneal > 0 .and. istep + 1 == istep_anneal_next) then
         ianneal = ianneal + 1
         tempK = anneal_tempK(ianneal)

         call set_md_coef()
         call set_bp_map()

         if (ianneal < nanneal) then
            istep_anneal_next = anneal_step(ianneal+1)
         else
            opt_anneal = 0
         endif
      endif

      if (flg_progress) then
         if (mod(istep, step_progress) == 0) then
            call progress_update(istep, nstep)
         endif
      endif

      if (stop_wall_time_sec > 0 .and. wall_time_sec() > stop_wall_time_sec) then
         flg_stop = .True.
         print '(a)', 'Wall-clock time limit reached. Stop the job after writing rst file.'
      endif

      if (mod(istep, nstep_save_rst) == 0 .or. flg_stop) then
         call write_rst()
      endif

      if (flg_stop) exit

   enddo

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
