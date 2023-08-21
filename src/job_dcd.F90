subroutine job_dcd()

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_idx, only : ENE
   use pbc, only : flg_pbc, pbc_box, pbc_box_half
   use var_top, only : nmp
   use var_state, only : xyz, energies, flg_bp_energy
   use var_potential, only: flg_stage
   use var_io, only : hdl_dcd, hdl_out, cfile_dcd_in, cfile_out
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   integer :: i, nframe, istat, nmp_dcd, icol
   type(file_dcd) :: fdcd
   character(len=29) :: out_fmt

   flg_bp_energy = .False.

   fdcd = file_dcd(hdl_dcd, cfile_dcd_in, DCD_OPEN_MODE%READ)

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
   
   !Format in .out file
   write(out_fmt, '(a15,i2,a12)') '(i10, 1x, f6.2,', ENE%MAX+1, '(1x, g13.6))'
   allocate(xyz(3, nmp))

   open(hdl_out, file = cfile_out, status = 'replace', action = 'write', form='formatted')
   ! Write header to .out file
   write(hdl_out, '(a)', advance='no') '#(1)nframe (2)T   (3)Ekin       (4)Epot       (5)Ebond      (6)Eangl      (7)Edih      '
                                       !1234567890 123456 1234567890123 1234567890123 1234567890123 1234567890123 1234567890123'
   write(hdl_out, '(a)') ' (8)Ebp        (9)Eexv       (10)Eele'
                         ! 1234567890123 1234567890123 1234567890123
   icol = 10
   if (flg_stage) then
           icol = icol + 1
           write(hdl_out, '(a,i2,a)', advance='no') ' (', icol, ')Estage   '
                                               ! 1   23     4567890123
   endif
   write(hdl_out, *) ''

   nframe = 0
   do
      call fdcd%read_onestep(nmp, xyz, istat)
      if (istat == iostat_end) exit
      nframe = nframe + 1

      call energy()
      
      !write(hdl_out, '(i10, 1x, f6.2, 8(1x,g13.6))') nframe, 0.0, 0.0, (energies(i), i=0,ENE%MAX)
      write(hdl_out, out_fmt, advance='no') nframe, 0.0, 0.0, (energies(i), i=0,ENE%ELE)
      if (flg_stage) then
         write(hdl_out, '(1x, g13.6)', advance='no') energies(ENE%STAGE)
      endif
      write(hdl_out, *) '' 

   enddo

   write(*,*) '#nframe:', nframe

   call fdcd%close()

   close(hdl_out)

endsubroutine job_dcd
