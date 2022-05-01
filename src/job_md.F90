subroutine job_md()

   use, intrinsic :: iso_fortran_env, Only : iostat_end, INT64
   use const
   use const_phys, only : KCAL2JOUL, N_AVO, PI, BOLTZ_KCAL_MOL
   use const_idx, only : ENE, SEQT
   use progress, only : progress_init, progress_update
   use pbc, only : pbc_box, set_pbc_size, flg_pbc
   use var_top, only : nmp, seq, mass, lmp_mp, ichain_mp
   use var_state, only : viscosity_Pas, xyz,  energies, forces, dt, velos, accels, tempK, nstep, nstep_save, &
                         nl_margin, Ekinetic, &
                         flg_variable_box, variable_box_step, variable_box_change, &
                         opt_anneal, nanneal, anneal_tempK, anneal_step
   use var_potential, only : wca_nl_cut2, wca_sigma, bp_nl_cut2, bp_cutoff
   use var_io, only : flg_progress, step_progress, hdl_dcd, hdl_out, cfile_prefix, cfile_out, cfile_pdb_ini, cfile_xyz_ini
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   integer(INT64) :: istep
   integer :: i, imp
   integer :: ianneal
   integer(INT64) :: istep_anneal_next
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

   ! Function
   real(PREC) :: rnd_boxmuller

   allocate(mass(nmp))
   allocate(xyz(3, nmp))
   allocate(forces(3, nmp))

   !! Initial structure from PDB or XYZ
   if (len(cfile_pdb_ini) > 0) then
      call read_pdb(cfile_pdb_ini, nmp, xyz)

   else if (len(cfile_xyz_ini) > 0) then
      call read_xyz(cfile_xyz_ini, nmp, xyz)

   else
      error stop 'Initial structure not found in job_md'
   endif

   !! Open DCD file
   cfile_dcd_out = trim(cfile_prefix) // '.dcd'
   write(*,*) 'Opening dcd file to write: ', trim(cfile_dcd_out)
   fdcd = file_dcd(hdl_dcd, cfile_dcd_out, DCD_OPEN_MODE%WRITE)

   ! set PBC box
   if (flg_pbc) then
      fdcd%flg_unitcell = .True.
      fdcd%box(:) = pbc_box(:)
   endif

   call fdcd%write_header(nmp)

   !! Set up variables for dynamics
   radius = 10.0
   v = viscosity_Pas * sqrt(1.0e3_PREC / KCAL2JOUL) * N_AVO * 1.0e-20_PREC
   write(*,*) 'v =', v
   fric(:) = 6.0e0_PREC * PI * v * radius
   write(*,*) 'fric =', fric(1)

   !v = 0.5 ! ps^(-1)
   !write(*,*) 'v =', v

   write(*,*) 'mass(A) = 328.212'
   write(*,*) 'mass(G) = 344.212'
   write(*,*) 'mass(C) = 304.182'
   write(*,*) 'mass(U) = 305.164'
           
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

   allocate(accels(3, nmp))
   allocate(velos(3, nmp))

   do imp = 1, nmp
      c1 = sqrt(tempK * BOLTZ_KCAL_MOL / mass(imp))
      do i = 1, 3
         velos(i, imp) = c1 * rnd_boxmuller()
      enddo
   enddo
        
   !! Initial accel
   !if (flg_rst) then
   !   call read_rst(RSTBLK%ACCEL)
   !else if (istep_sim == 1 .OR. inmisc%i_reset_struct == 1) then
      do imp = 1, nmp
         do i = 1, 3
            accels(i, imp) = md_coef(1, imp) * rnd_boxmuller()
         enddo
      end do
   !endif

   wca_nl_cut2 = (wca_sigma + nl_margin) ** 2
   bp_nl_cut2 = (bp_cutoff + nl_margin) ** 2
   call neighbor_list()
   xyz_move(:,:) = 0.0e0_PREC

   ! Setting up Simulated Annealing
   if (opt_anneal > 0) then
      tempK = anneal_tempK(1)
      ianneal = 1
      istep_anneal_next = anneal_step(2)
   endif

   open(hdl_out, file = cfile_out, status = 'replace', action = 'write', form='formatted')
   write(hdl_out, '(a)') '#(1)nframe (2)T   (3)Ekin       (4)Epot       (5)Ebond      (6)Eangl      (7)Ebp        (8)Eexv'
                         !1234567890 123456 1234567890123 1234567890123 1234567890123 1234567890123 1234567890123 1234567890123

   ! Write the initial coordinates
   !call fdcd%write_onestep(nmp, xyz)     ! No initial frame for DCD
   call energy()
   call energy_kinetic()
   write(hdl_out, '(i10, 1x, f6.2, 6(1x,g13.6))') 0, tempK, Ekinetic, (energies(i), i=0,ENE%MAX)

   if (flg_progress) then
      call progress_init(0_INT64)
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!! Time integration
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   do istep = 1, nstep

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
         write(hdl_out, '(i10, 1x, f6.2, 6(1x,g13.6))') istep, tempK, Ekinetic, (energies(i), i=0,ENE%MAX)
         call fdcd%write_onestep(nmp, xyz)
      endif

      if (flg_variable_box) then
         if (mod(istep, variable_box_step) == 0) then
            call set_pbc_size(pbc_box(:) + variable_box_change(:))
            fdcd%box(:) = pbc_box(:)

            call neighbor_list()
            xyz_move(:,:) = 0.0e0_PREC

            write(*,'(a,i10,a,f8.3)') 'Box size updated: step = ',istep, ', box size = ', pbc_box(1)
         endif
      endif

      if (flg_progress) then
         if (mod(istep, step_progress) == 0) then
            call progress_update(istep, nstep)
         endif
      endif

      if (opt_anneal > 0 .and. istep + 1 == istep_anneal_next) then
         ianneal = ianneal + 1
         tempK = anneal_tempK(ianneal)

         if (ianneal < nanneal) then
            istep_anneal_next = anneal_step(ianneal+1)
         else
            opt_anneal = 0
         endif
      endif
   enddo

   call fdcd%close()

   close(hdl_out)

endsubroutine job_md
