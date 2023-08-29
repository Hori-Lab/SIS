subroutine job_dcd()

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_idx, only : ENE
   use pbc, only : flg_pbc, pbc_box, pbc_box_half
   use var_top, only : nmp
   use var_state, only : xyz, energies, flg_bp_energy, tempK
   use var_potential, only : flg_stage
   use var_io, only : hdl_out, cfile_dcd_in, iopen_hdl
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   integer :: i, nframe, istat, nmp_dcd, icol
   type(file_dcd) :: fdcd
   integer, parameter :: IREP = 1
   character(len=37) :: out_fmt

   ! Format in .out file
   write(out_fmt, '(a15,i2,a12)') '(i10, 1x, f6.2,', ENE%ELE+2, '(1x, g13.6))'

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

   write(hdl_out(IREP), '(a)', advance='no') '#(1)nframe (2)T   (3)Ekin       (4)Epot       (5)Ebond     '
                                             !1234567890 123456 1234567890123 1234567890123 1234567890123'
   write(hdl_out(IREP), '(a)', advance='no') ' (6)Eangl      (7)Edih      '
                                             ! 1234567890123 1234567890123'
   write(hdl_out(IREP), '(a)', advance='no') ' (8)Ebp        (9)Eexv       (10)Eele'
                                             ! 1234567890123 1234567890123 1234567890123
   icol = 10
   if (flg_stage) then
      icol = icol + 1
      write(hdl_out(IREP), '(a,i2,a)', advance='no') ' (', icol, ')Estage   '
                                                     ! 1   23     4567890123
   endif
   write(hdl_out(IREP), '(a)') ''

   nframe = 0
   do
      call fdcd%read_onestep(nmp, xyz(:,:,IREP), istat)
      if (istat == iostat_end) exit
      nframe = nframe + 1

      call energy_sumup(IREP, tempK, energies(0:ENE%MAX, IREP))
      
      write(hdl_out(IREP), out_fmt, advance='no') nframe, 0.0, 0.0, (energies(i, IREP), i=0,ENE%ELE)
      if (flg_stage) then
         write(hdl_out(IREP), '(1x, g13.6)', advance='no') energies(ENE%STAGE, IREP)
      endif
      write(hdl_out(IREP), '(a)') ''

      call write_bp(IREP, tempK)
   enddo

   write(*,*) '#nframe:', nframe

   call fdcd%close()

endsubroutine job_dcd
