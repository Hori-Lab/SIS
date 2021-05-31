program sn

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_idx, only : ENE
   use var_state, only : nmp, nchains, nmp_chain, xyz, energies
   use var_io, only : hdl_dcd, hdl_out
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   character(CHAR_FILE_PATH) cfile_sis, cfile_dcd, cfile_out

   type(file_dcd) :: fdcd

   integer :: i
   integer :: istat
   integer :: nframe

   integer :: iarg
   integer :: iargc

   iarg = iargc()
   if (iarg /= 3) then
      write(6,*) 'Usage: PROGRAM [sisinfo file] [dcd file] [output prefix]'
      stop (2) 
   end if

   call getarg(1, cfile_sis) 
   call getarg(2, cfile_dcd) 
   call getarg(3, cfile_out) 

   call read_sisinfo(cfile_sis)

   fdcd = file_dcd(hdl_dcd, cfile_dcd, DCD_OPEN_MODE%READ)

   call fdcd%read_header()

   open(hdl_out, file = cfile_out, status = 'unknown', action = 'write')

   call fdcd%read_nmp(nmp, istat)

   if (nmp /= nmp_chain * nchains) then
      write(*,*) "nmp = ", nmp, " is inconsistent with nmp_chain = ", nmp_chain, " and nchains = ", nchains
      stop (2) 
   endif

   allocate(xyz(3, nmp))

   nframe = 0
   do
      call fdcd%read_onestep(nmp, xyz, istat)
      if (istat == iostat_end) exit
      nframe = nframe + 1

      call energy()
      
      write(*,*) (energies(i), i=1,ENE%MAX), energies(ENE%TOTAL)

   enddo

   write(*,*) 'nframe =', nframe

   call fdcd%close()

   close(hdl_out)

   deallocate(xyz)

   stop

end program sn
