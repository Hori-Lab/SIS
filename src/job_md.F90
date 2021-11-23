subroutine job_md()

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_phys, only : KCAL2JOUL, N_AVO, PI, BOLTZ_KCAL_MOL
   use const_idx, only : ENE
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, pbc_box, pbc_box_half, flg_pbc, mass
   use var_state, only : viscosity_Pas, xyz,  energies, forces, dt, velos, accels, tempK, nstep, nstep_save, Ekinetic
   use var_io, only : hdl_dcd, hdl_out, cfile_prefix, cfile_out, cfile_pdb_ini
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   integer :: i, istat
   integer :: istep, imp
   real(PREC) :: dxyz(3)
   !real(PREC) :: xyz_move(3, nmp)
   !real(PREC) :: velo(3, nmp), accel1(3, nmp), accel2(3, nmp)
   type(file_dcd) :: fdcd
   real(PREC) :: fric, radius
   real(PREC) :: v, c1, c2, md_coef(4, nmp)
   real(PREC) :: rnd_bm(3, nmp)
   real(PREC) :: accels_pre(3)
   character(CHAR_FILE_PATH), save :: cfile_dcd_out

   ! Function
   real(PREC) :: rnd_boxmuller

   allocate(mass(nmp))
   allocate(xyz(3, nmp))
   allocate(forces(3, nmp))

   !! Initial structure from PDB
   if (len(cfile_pdb_ini) < 1) then
      write(*,*) 'PDB for the initial structure is not specified'
      stop
   endif

   call read_pdb(cfile_pdb_ini, nmp, xyz)

   !! Open DCD file
   cfile_dcd_out = trim(cfile_prefix) // '.dcd'
   write(*,*) 'Opening dcd file to write: ', trim(cfile_dcd_out)
   fdcd = file_dcd(hdl_dcd, cfile_dcd_out, DCD_OPEN_MODE%WRITE)

   call fdcd%write_header(nmp)

   ! set PBC box
   flg_pbc = .True.
   fdcd%box(:) = pbc_box(:)

   call fdcd%write_onestep(nmp, xyz)

   !! Set up variables for dynamics
   radius = 10.0
   mass(:) = 300.0
   v = viscosity_Pas * sqrt(1.0e3_PREC / KCAL2JOUL) * N_AVO * 1.0e-20_PREC
   write(*,*) 'v =', v
   fric = 6.0e0_PREC * PI * v * radius
   write(*,*) 'fric =', fric
           
   do imp = 1, nmp
      !! sqrt(b) = sqrt(1 / (1 + gamma h / 2m))
      c1 = 0.5 * dt * fric / mass(imp)
      c2 = sqrt(1.0e0_PREC / (1.0_PREC + c1))

      ! md_coef(1) = sqrt(b) / 2m * sqrt(2 gamma kT h)
      md_coef(1, imp) = 0.5_PREC * c2 / mass(imp) * sqrt(2.0_PREC * fric * BOLTZ_KCAL_MOL * tempK * dt)
      ! md_coef(2) = a
      md_coef(2, imp) = (1.0_PREC - c1) / (1.0_PREC + c1)
      ! md_coef(3) = sqrt(b) h / m
      md_coef(3, imp) = c2 * dt / mass(imp)
      ! md_coef(4) = sqrt(b) h
      md_coef(4, imp) = c2 * dt
   enddo

   allocate(accels(3, nmp))
   allocate(velos(3, nmp))
        
   !! Initial accel
   !if (flg_rst) then
   !   call read_rst(RSTBLK%ACCEL)
   !else if (istep_sim == 1 .OR. inmisc%i_reset_struct == 1) then
      do imp = 1, nmp
         accels(1:3, imp) = md_coef(1, imp) * rnd_boxmuller()
      end do
   !endif


   open(hdl_out, file = cfile_out, status = 'replace', action = 'write', form='formatted')
   write(hdl_out, '(a)') '#(1)nframe (2)Ekin       (3)Epot       (4)Ebond      (5)Eangl      (6)Ebp        (7)Eexv'
                         !1234567890 1234567890123 1234567890123 1234567890123 1234567890123 1234567890123 1234567890123

   do istep = 1, nstep
      
      do imp= 1, nmp
         do i = 1, 3
            rnd_bm(i, imp) = rnd_boxmuller()
         enddo
      enddo

      call force()

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
         !pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, imp, irep) + dxyz(1:3)
         !xyz_move(1:3,imp) = xyz_move(1:3,imp) + dxyz(1:3)
      end do

      if (mod(istep, nstep_save) == 0) then
          call energy()
          call energy_kinetic()
          write(hdl_out, '(i10, 6(1x,g13.6))') istep, Ekinetic, (energies(i), i=0,ENE%MAX)
          call fdcd%write_onestep(nmp, xyz)
      endif

   enddo

   call fdcd%close()

   close(hdl_out)

endsubroutine job_md
