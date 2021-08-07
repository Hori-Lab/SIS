subroutine job_check_force()

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_idx, only : ENE
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, pbc_box, pbc_box_half, flg_pbc
   use var_state, only : xyz, energies, forces
   use var_io, only : hdl_out, cfile_out, cfile_pdb_ini

   implicit none

   integer :: i, imp, ixyz
   real(PREC) :: r, xyz_save, e_save
   real(PREC) :: f_energy(3), diff(3)
   real(PREC), parameter :: small = 0.0001

   if (len(cfile_pdb_ini) < 1) then
      write(*,*) 'PDB for the initial structure is not specified'
      stop
   endif

   write(*,*) 'Allocating xyz and forces, nmp=', nmp
   allocate(xyz(3, nmp))
   allocate(forces(3, nmp))

   write(*,*) 'calling PDB done'
   call read_pdb(cfile_pdb_ini, nmp, xyz)

   write(*,*) 'Loading PDB done'

   open(hdl_out, file = cfile_out, status = 'replace', action = 'write', form='formatted')
   !write(hdl_out, '(a)') '#(1)nframe  (2)Etotal  (3)Ebond   (4)Eangl   (5)Ebp   (6)Eele'

   do imp = 1, nmp
      do ixyz = 1, 3
         call random_number(r)
         xyz(ixyz, imp) = xyz(ixyz, imp) + small * r
      enddo
   enddo

   call force()

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

      write(hdl_out, "('imp=', i4, 2x, 'f_energy', 3(1x,g10.4), 2x, 'f_force=', 3(1x,g10.4))") &
             imp, (f_energy(i), i = 1, 3), (forces(i, imp), i = 1, 3)

      diff(1:3) = f_energy(1:3) - forces(1:3, imp)
      write(hdl_out, "('f_energy - f_force', i4, 3f10.4)") imp, (diff(i), i = 1, 3)

      if(abs(diff(1)) > small .or. abs(diff(2)) > small .or. abs(diff(3)) > small) then
         write (*, "('f_energy - f_force', 1x, i6, 3(1x,f10.4))") imp, (diff(i), i = 1, 3)
      end if

   enddo

   close(hdl_out)

endsubroutine job_check_force
