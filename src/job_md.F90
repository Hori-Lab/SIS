subroutine job_md()

   use, intrinsic :: iso_fortran_env, Only : iostat_end, INT64
   use const
   use const_phys, only : KCAL2JOUL, N_AVO, PI, BOLTZ_KCAL_MOL
   use const_idx, only : ENE, SEQT, RSTBLK, BPT
   use progress, only : progress_init, progress_update, wall_time_sec
   use pbc, only : pbc_box, set_pbc_size, flg_pbc
   use var_top, only : nmp, seq, mass, lmp_mp, ichain_mp
   use var_state, only : restarted, flg_bp_energy, &
                         viscosity_Pas, xyz,  energies, forces, dt, velos, accels, tempK, &
                         nstep, nstep_save, nstep_save_rst, stop_wall_time_sec, fix_com_origin, &
                         nl_margin, Ekinetic, &
                         flg_variable_box, variable_box_step, variable_box_change, &
                         opt_anneal, nanneal, anneal_tempK, anneal_step, &
                         istep, ianneal, istep_anneal_next
   use var_potential, only : wca_nl_cut2, stage_sigma, wca_sigma, bp_nl_cut2, ele_cutoff, ele_nl_cut2, bp_paras, &
                             bp_cutoff_energy, bp_cutoff_dist
   use var_io, only : flg_progress, step_progress, hdl_dcd, hdl_out, cfile_prefix, cfile_out, cfile_pdb_ini, cfile_xyz_ini
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   integer :: i, imp
   real(PREC) :: dxyz(3)
   real(PREC) :: xyz_move(3, nmp)
   !real(PREC) :: velo(3, nmp), accel1(3, nmp), accel2(3, nmp)
   type(file_dcd) :: fdcd
   real(PREC) :: fric(nmp), radius
   real(PREC) :: v, c1, c2, md_coef(4, nmp)
   real(PREC) :: rnd_bm(3, nmp)
   real(PREC) :: accels_pre(3)
   real(PREC) :: d2, d2max, d2max_2nd
   character(CHAR_FILE_PATH), save :: cfile_dcd_out
   logical :: flg_stop

   ! Function
   real(PREC) :: rnd_boxmuller

   allocate(mass(nmp))
   allocate(xyz(3, nmp))
   allocate(forces(3, nmp))
   allocate(accels(3, nmp))
   allocate(velos(3, nmp))

   ! set PBC box
   if (flg_pbc) then
      fdcd%flg_unitcell = .True.
      fdcd%box(:) = pbc_box(:)
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
            stop (2)
      endselect
   enddo

   call set_md_coef()


   !! Set up the initial state
   if (restarted) then
      call read_rst(RSTBLK%XYZ)
      call read_rst(RSTBLK%VELO)
      call read_rst(RSTBLK%ACCEL)

   else
      ! Load initial coordinates from PDB or XYZ
      if (len(cfile_pdb_ini) > 0) then
         call read_pdb(cfile_pdb_ini, nmp, xyz)

      else if (len(cfile_xyz_ini) > 0) then
         call read_xyz(cfile_xyz_ini, nmp, xyz)

      else
         error stop 'Initial structure not found in job_md'
      endif

      ! Set up initial velocities by Maxwell–Boltzmann distribution
      do imp = 1, nmp
         c1 = sqrt(tempK * BOLTZ_KCAL_MOL / mass(imp))
         do i = 1, 3
            velos(i, imp) = c1 * rnd_boxmuller()
         enddo
      enddo

      ! Set up initial accerelations
      do imp = 1, nmp
         do i = 1, 3
            accels(i, imp) = md_coef(1, imp) * rnd_boxmuller()
         enddo
      end do
   endif

   ! Neighbor list
   wca_nl_cut2 = (wca_sigma + nl_margin) ** 2
   bp_nl_cut2 = (bp_cutoff_dist + nl_margin) ** 2
   ele_nl_cut2 = (ele_cutoff + nl_margin) ** 2

   print '(a)', 'Potential and neighbor list cutoffs'
   print '(a,f10.3)', 'bp_cutoff_energy = ', bp_cutoff_energy
   print '(a,f10.3)', 'bp_cutoff_ddist(GC) = ', bp_paras(BPT%GC)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_ddist(AU) = ', bp_paras(BPT%AU)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_ddist(GU) = ', bp_paras(BPT%GU)%cutoff_ddist
   print '(a,f10.3)', 'bp_cutoff_dist(for neighbor list) = ', bp_cutoff_dist
   print '(a,f10.3)', 'wca_sigma = ', wca_sigma
   print '(a,f10.3)', 'stage_sigma = ', stage_sigma
   print '(a,f10.3)', 'ele_cutoff = ', ele_cutoff
   print '(a,f10.3)', 'nl_margin = ', nl_margin
   print *

   call neighbor_list()
   xyz_move(:,:) = 0.0e0_PREC

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
   flg_bp_energy = .False.
   call energy()
   call energy_kinetic()

   ! Open DCD file and write the header
   cfile_dcd_out = trim(cfile_prefix) // '.dcd'
   print '(2a)', '# Opening dcd file to write: ', trim(cfile_dcd_out)
   fdcd = file_dcd(hdl_dcd, cfile_dcd_out, DCD_OPEN_MODE%WRITE)

   call fdcd%write_header(nmp)

   ! Open .out file
   open(hdl_out, file = cfile_out, status = 'replace', action = 'write', form='formatted')
   write(hdl_out, '(a)', advance='no') '#(1)nframe (2)T   (3)Ekin       (4)Epot       (5)Ebond      (6)Eangl      (7)Edih      '
                                       !1234567890 123456 1234567890123 1234567890123 1234567890123 1234567890123 1234567890123'
   write(hdl_out, '(a)') ' (8)Ebp        (9)Eexv       (10)Eele      (11)Estage'
                         ! 1234567890123 1234567890123 1234567890123 1234567890123

   ! Output initial structure
   if (restarted) then
      ! Only STDOUT if restarted. No DCD output.
      print '(a)', '##### Energies at the beginning'
      print '(a)', '#(1)nframe (2)T   (3)Ekin       (4)Epot       '
      print '(i10, 1x, f6.2, 2(1x,g13.6))', istep, tempK, Ekinetic, energies(0)
      print '(a)', '(5)Ebond      (6)Eangl      (7)Edih       (8)Ebp        (9)Eexv       (10)Eele      (11)Estage'
      print '(7(1x,g13.6))', (energies(i), i=1, ENE%MAX)
      print *

   else
      ! At istep = 0 (not restarted), write both .out and DCD
      write(hdl_out, '(i10, 1x, f6.2, 9(1x,g13.6))') istep, tempK, Ekinetic, (energies(i), i=0,ENE%MAX)
      call fdcd%write_onestep(nmp, xyz, fix_com_origin)
   endif

   if (flg_progress) then
      call progress_init(istep)
   endif

   flg_stop = .False.

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Time integration
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do while (istep < nstep)

      istep = istep + 1

      d2max = 0.0e0_PREC
      d2max_2nd = 0.0e0_PREC
      do imp = 1, nmp
         d2 = dot_product(xyz_move(1:3,imp), xyz_move(1:3,imp))
         if (d2 > d2max) then
            d2max_2nd = d2max
            d2max = d2
         else if (d2 > d2max_2nd) then
            d2max_2nd = d2
         endif
      enddo

      if (sqrt(d2max) + sqrt(d2max_2nd) > nl_margin) then
         call neighbor_list()
         xyz_move(:,:) = 0.0e0_PREC
      endif

      call force()

      if (maxval(forces(:,:)) > 100.0) then
         print *, 'istep, force ', istep, maxval(forces)
         flush(6)
      endif

      do imp= 1, nmp
         do i = 1, 3
            rnd_bm(i, imp) = rnd_boxmuller()
         enddo
      enddo

      !$omp parallel do private(accels_pre, dxyz)
      do imp = 1, nmp
         !if(fix_mp(imp)) cycle

         !! a = (1 - gamma h / 2m) / (1 + gamma h / 2m)
         !! b = 1 / (1 + gamma h / 2m)

         ! beta(t + h) with the associated coefficient
         accels_pre(1:3) = md_coef(1, imp) * rnd_bm(1:3, imp)
         ! md_coef(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)

         ! v(t + 1/2h) update the half-step velocity
         velos(1:3, imp) =  md_coef(2, imp) * velos(1:3, imp)  &
                          + md_coef(3, imp) * forces(1:3, imp) &
                          + (accels(1:3, imp) + accels_pre(1:3))
         ! md_coef(2) = a
         ! md_coef(3) = sqrt(b) h / m

         ! beta(t) <= beta(t+h) (incluing the coefficient) save for the next iteration
         accels(1:3, imp) = accels_pre(1:3)

         dxyz(1:3) =  md_coef(4, imp) * velos(1:3, imp)
         ! md_coef(4) = sqrt(b) h

         xyz(1:3, imp) = xyz(1:3, imp) + dxyz(1:3)
         xyz_move(1:3,imp) = xyz_move(1:3,imp) + dxyz(1:3)
      end do
      !$omp end parallel do

      if (mod(istep, nstep_save) == 0) then
         call energy()
         call energy_kinetic()
         write(hdl_out, '(i10, 1x, f6.2, 8(1x,g13.6))') istep, tempK, Ekinetic, (energies(i), i=0,ENE%MAX)
         call fdcd%write_onestep(nmp, xyz, fix_com_origin)
      endif

      if (flg_variable_box) then
         if (mod(istep, variable_box_step) == 0) then
            call set_pbc_size(pbc_box(:) + variable_box_change(:))
            fdcd%box(:) = pbc_box(:)

            call neighbor_list()
            xyz_move(:,:) = 0.0e0_PREC

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

   call fdcd%close()

   close(hdl_out)

contains

   subroutine set_md_coef()
      do imp = 1, nmp
         !fric(imp) = v * mass(imp)
         !if (imp == 1) then
         !   write(*,*) 'fric = ', fric(imp)
         !endif

         !! sqrt(b) = sqrt(1 / (1 + gamma h / 2m))
         c1 = 0.5 * dt * fric(imp) / mass(imp)
         c2 = sqrt(1.0e0_PREC / (1.0_PREC + c1))

         ! md_coef(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)
         md_coef(1, imp) = 0.5_PREC * c2 / mass(imp) * sqrt(2.0_PREC * fric(imp) * BOLTZ_KCAL_MOL * tempK * dt)
         ! md_coef(2) = a
         md_coef(2, imp) = (1.0_PREC - c1) / (1.0_PREC + c1)
         ! md_coef(3) = sqrt(b) h / m
         md_coef(3, imp) = c2 * dt / mass(imp)
         ! md_coef(4) = sqrt(b) h
         md_coef(4, imp) = c2 * dt
      enddo
   endsubroutine set_md_coef

endsubroutine job_md
