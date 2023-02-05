subroutine job_check_force()

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_idx, only : ENE
   use pbc, only : flg_pbc, pbc_wrap
   use var_top, only : nmp
   use var_state, only : xyz, energies, forces, flg_bp_energy
   use var_potential, only: flg_stage
   use var_io, only : hdl_out, cfile_out, cfile_pdb_ini, cfile_xyz_ini
   use mt19937_64

   implicit none

   integer :: i, imp, ixyz
   real(PREC) :: xyz_save, e_save
   real(PREC) :: f_energy(3), diff(3)
   real(PREC), parameter :: small = 0.001

   write(6,*) 'Starting job_check_force'

   allocate(xyz(3, nmp))
   allocate(forces(3, nmp))

   if (len(cfile_pdb_ini) > 0) then
      call read_pdb(cfile_pdb_ini, nmp, xyz)

   else if (len(cfile_xyz_ini) > 0) then
      call read_xyz(cfile_xyz_ini, nmp, xyz)

   else
      error stop 'Initial structure not found in job_check_force.'
   endif
   
   if (flg_stage) call check_int_stage()

   if (flg_pbc) then
      call pbc_wrap()
   endif

   open(hdl_out, file = cfile_out, status = 'replace', action = 'write', form='formatted')

   ! Slightly shift the structure to make sure it's not at the energy minimum
   do imp = 1, nmp
      do ixyz = 1, 3
         xyz(ixyz, imp) = xyz(ixyz, imp) + small * genrand64_real3()
      enddo
   enddo

   call force()

   flg_bp_energy = .False.

   do imp = 1, nmp
      do ixyz = 1, 3

         xyz_save = xyz(ixyz, imp)
         xyz(ixyz, imp) = xyz_save + small
         call energy()
         e_save = energies(ENE%TOTAL)

         xyz(ixyz, imp) = xyz_save - small
         call energy()

         f_energy(ixyz) = - (e_save - energies(ENE%TOTAL)) / (2.0 * small)
         xyz(ixyz, imp) = xyz_save
      enddo

      write(hdl_out, "('imp=', i6, 2x, 'f_energy', 3(1x,g11.4), 2x, 'f_force=', 3(1x,g11.4))") &
             imp, (f_energy(i), i = 1, 3), (forces(i, imp), i = 1, 3)

      diff(1:3) = f_energy(1:3) - forces(1:3, imp)
      write(hdl_out, "('f_energy - f_force', i6, 3f11.4)") imp, (diff(i), i = 1, 3)

      if(abs(diff(1)) > small .or. abs(diff(2)) > small .or. abs(diff(3)) > small) then
         write (*, "('f_energy - f_force', 1x, i6, 3(1x,f11.4))") imp, (diff(i), i = 1, 3)
      end if

   enddo

   close(hdl_out)

   write(6,*) 'Done job_check_force'

endsubroutine job_check_force
