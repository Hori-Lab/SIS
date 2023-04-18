subroutine job_check_force()

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_idx, only : ENE
   use pbc, only : flg_pbc, pbc_wrap
   use var_top, only : nmp
   use var_state, only : xyz, energies, flg_bp_energy, mts
   use var_potential, only: flg_stage
   use var_io, only : hdl_out, cfile_pdb_ini, cfile_xyz_ini
   !use mt19937_64
   use mt_stream

   implicit none

   integer :: i, imp, ixyz
   real(PREC) :: xyz_save, e_save
   real(PREC) :: f_energy(3), diff(3)
   real(PREC), allocatable :: forces(:,:)
   integer, parameter :: IREP = 1
   real(PREC), parameter :: small = 0.001

   print '(a)', 'Starting job_check_force'

   allocate(forces(3, nmp))
   allocate(energies(0:ENE%MAX, IREP))

   if (len(cfile_pdb_ini) > 0) then
      call read_pdb(cfile_pdb_ini, nmp, xyz)

   else if (len(cfile_xyz_ini) > 0) then
      call read_xyz(cfile_xyz_ini, nmp, xyz)

   else
      error stop 'Initial structure not found in job_check_force.'
   endif
   
   if (flg_stage) call check_int_stage()

   if (flg_pbc) then
      call pbc_wrap(IREP)
   endif

   ! Slightly shift the structure to make sure it's not at the energy minimum
   do imp = 1, nmp
      do ixyz = 1, 3
         !xyz(ixyz, imp, IREP) = xyz(ixyz, imp, IREP) + small * genrand64_real3()
         xyz(ixyz, imp, IREP) = xyz(ixyz, imp, IREP) + small * genrand_double3(mts(0))
      enddo
   enddo

   call force(IREP, forces)

   flg_bp_energy = .False.

   do imp = 1, nmp
      do ixyz = 1, 3

         xyz_save = xyz(ixyz, imp, IREP)
         xyz(ixyz, imp, IREP) = xyz_save + small
         call energy(IREP, energies)
         e_save = energies(ENE%TOTAL, IREP)

         xyz(ixyz, imp, IREP) = xyz_save - small
         call energy(IREP, energies)

         f_energy(ixyz) = - (e_save - energies(ENE%TOTAL, IREP)) / (2.0 * small)
         xyz(ixyz, imp, IREP) = xyz_save
      enddo

      write(hdl_out(IREP), "('imp=', i6, 2x, 'f_energy', 3(1x,g11.4), 2x, 'f_force=', 3(1x,g11.4))") &
             imp, (f_energy(i), i = 1, 3), (forces(i, imp), i = 1, 3)

      diff(1:3) = f_energy(1:3) - forces(1:3, imp)
      write(hdl_out(IREP), "('f_energy - f_force', i6, 3f11.4)") imp, (diff(i), i = 1, 3)

      if(abs(diff(1)) > small .or. abs(diff(2)) > small .or. abs(diff(3)) > small) then
         write (*, "('f_energy - f_force', 1x, i6, 3(1x,f11.4))") imp, (diff(i), i = 1, 3)
      end if

   enddo

   print '(a)', 'Done job_check_force'

endsubroutine job_check_force
