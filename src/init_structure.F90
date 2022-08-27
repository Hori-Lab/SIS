subroutine init_structure()

   use const, only : PREC
   use var_io, only : flg_in_pdb, flg_in_xyz, cfile_pdb_ini, cfile_xyz_ini, flg_gen_init_struct
   use var_top, only : nmp
   use var_state, only : xyz
   use var_replica, only : nrep_proc
   use var_parallel

   implicit none

   integer :: ierr, irep
   real(PREC) :: origin(3) = (/0.0_PREC, 0.0_PREC, 0.0_PREC/)

   if (flg_gen_init_struct) then

      do irep = 1, nrep_proc
         call gen_random_coil(irep, nmp, xyz(:,:,irep), origin)
      enddo

   else

      if (myrank == 0) then
         ! Load initial coordinates from PDB or XYZ
         if (flg_in_pdb) then
            call read_pdb(cfile_pdb_ini, nmp, xyz(:,:,1))

         else if (flg_in_xyz) then
            call read_xyz(cfile_xyz_ini, nmp, xyz(:,:,1))

         else
            print *, 'Error: Logical error in init_structure.'
            ! These flags should have been checked in read_input.
            call sis_abort()
         endif
      endif

#ifdef PAR_MPI
      call MPI_BCAST(xyz(:,:,1), 3*nmp, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
#endif

      do irep = 2, nrep_proc
         xyz(:,:,irep) = xyz(:,:,1)
      enddo
   endif

endsubroutine init_structure
