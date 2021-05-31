program sn

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_idx, only : ENE, SEQT
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, pbc_box, pbc_box_half
   use var_state, only : xyz, energies
   use var_io, only : hdl_dcd, hdl_out
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   character(CHAR_FILE_PATH) cfile_sis, cfile_dcd, cfile_out

   type(file_dcd) :: fdcd

   integer :: i, j, imp
   integer :: istat
   integer :: nframe
   integer :: nrepeat
   integer :: nmp_dcd

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

   nrepeat = 47
   nchains = 64
   allocate(nmp_chain(nchains))
   nmp_chain(:) = 3 * nrepeat
   nmp = sum(nmp_chain)
   allocate(seq(3*nrepeat, nchains))
   allocate(imp_chain(3*nrepeat, nchains))
   imp = 0
   do i = 1, nchains
      do j = 1, nrepeat
         seq(3*(j-1)+1, i) = SEQT%C
         seq(3*(j-1)+2, i) = SEQT%A
         seq(3*(j-1)+3, i) = SEQT%G
         imp_chain(3*(j-1)+1, i) = imp+1
         imp_chain(3*(j-1)+2, i) = imp+2
         imp_chain(3*(j-1)+3, i) = imp+3
         imp = imp + 3
      enddo
   enddo

   call read_sisinfo(cfile_sis)

   fdcd = file_dcd(hdl_dcd, cfile_dcd, DCD_OPEN_MODE%READ)

   call fdcd%read_header()

   open(hdl_out, file = cfile_out, status = 'unknown', action = 'write')

   call fdcd%read_nmp(nmp_dcd, istat)

   pbc_box(:) = fdcd%box(:)
   pbc_box_half = 0.5 * pbc_box(:)

   if (nmp_dcd /= nmp) then
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
