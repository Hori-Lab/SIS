subroutine job_dcd()

   use, intrinsic :: iso_fortran_env, Only : iostat_end, INT64
   use const
   use const_idx, only : ENE
   use pbc, only : flg_pbc, pbc_box, pbc_box_half
   use var_top, only : nmp
   use var_io, only : cfile_dcd_in, iopen_hdl
   use var_state, only : xyz, energies, flg_bp_energy, tempK, nstep_bp_MC, flg_bp_MC, Ekinetic, rg
   use dcd, only : file_dcd, DCD_OPEN_MODE
   use var_replica, only : nrep_proc
#ifdef DUMPFORCE
   use var_state, only : flg_step_dump_force
#endif

   implicit none

   integer(INT64) :: nframe
   integer :: istat, nmp_dcd
   type(file_dcd) :: fdcd
   integer, parameter :: IREP = 1

#ifdef DUMPFORCE
   real(PREC), allocatable :: forces(:, :)
   flg_step_dump_force = .True.
   allocate(forces(3, nmp))
#endif

   allocate(Ekinetic(nrep_proc))
   Ekinetic(:) = 0.0_PREC
   flg_bp_energy = .False.

   iopen_hdl = iopen_hdl + 1
   fdcd = file_dcd(iopen_hdl, cfile_dcd_in, DCD_OPEN_MODE%READ)

   call fdcd%read_header()

   call fdcd%read_nmp(nmp_dcd, istat)

   if (fdcd%box(1) > 0.0) then
      flg_pbc = .True.
      pbc_box(:) = fdcd%box(:)
      pbc_box_half = 0.5 * pbc_box(:)
      write(*,*) '#Box: ', pbc_box
   else
      flg_pbc = .False.
      write(*,*) '#No PBC'
   endif

   if (nmp_dcd /= nmp) then
      write(*,*) "Error: nmp_dcd /= nmp. nmp = ", nmp, ", nmp_dcd = ", nmp_dcd
      stop (2) 
   endif

   allocate(xyz(3, nmp, IREP))
   allocate(energies(0:ENE%MAX, IREP))
   allocate(rg(IREP))

   call write_out_header(IREP)

   flg_bp_MC = .True.
   if (nstep_bp_MC == 0) flg_bp_MC = .False.

   nframe = 0
   do
      call fdcd%read_onestep(nmp, xyz(:,:,IREP), istat)
      if (istat == iostat_end) exit
      nframe = nframe + 1

#ifdef DUMPFORCE
      call force(IREP, forces(:,:))
#endif

      call energy_sumup(IREP, tempK, energies(0:ENE%MAX, IREP))
      
      call write_out(IREP, nframe, IREP, tempK)

      call write_bp(IREP, tempK)
   enddo

   write(*,*) '#nframe:', nframe

   call fdcd%close()

endsubroutine job_dcd
