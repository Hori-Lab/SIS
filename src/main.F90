program sn

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_idx, only : ENE, SEQT
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, pbc_box, pbc_box_half, flg_pbc, ichain_mp
   use var_state, only : xyz, energies
   use var_io, only : hdl_dcd, hdl_out, flg_out_bp, hdl_bp, KIND_OUT_BP
   use dcd, only : file_dcd, DCD_OPEN_MODE

   implicit none

   character(CHAR_FILE_PATH) cfile_dcd, cfile_prefix, cfile_out, cfile_bp

   type(file_dcd) :: fdcd

   integer :: i, j, imp
   integer :: istat
   integer :: nframe
   integer :: nrepeat
   integer :: nmp_dcd

   character(500) :: cline

   if (command_argument_count() /= 4) then
      !write(6,*) 'Usage: PROGRAM [sisinfo file] [dcd file] [output prefix]'
      write(6,*) 'Usage: PROGRAM [nrepeat] [nchain] [dcd file] [output prefix]'
      stop (2) 
   end if

   call get_command_argument(1, cline)
   read(cline, *) nrepeat
   call get_command_argument(2, cline)
   read(cline, *) nchains
   call get_command_argument(3, cfile_dcd)  
   call get_command_argument(4, cfile_prefix)  

   cfile_out = trim(cfile_prefix) // '.out'

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Temporary hard-code the system setup
   flg_out_bp = .True.
   cfile_bp = trim(cfile_prefix) // '.bp'
   open(hdl_bp, file=cfile_bp, status='replace', action='write', form='unformatted',access='stream')
   write(hdl_bp) int(KIND_OUT_BP,kind=4)
   !nrepeat = 47
   !nchains =  1
   allocate(nmp_chain(nchains))
   nmp_chain(:) = 3 * nrepeat
   nmp = sum(nmp_chain)
   allocate(seq(3*nrepeat, nchains))
   allocate(imp_chain(3*nrepeat, nchains))
   allocate(ichain_mp(nmp))
   imp = 0
   do i = 1, nchains
      do j = 1, nrepeat
         seq(3*(j-1)+1, i) = SEQT%C
         seq(3*(j-1)+2, i) = SEQT%A
         seq(3*(j-1)+3, i) = SEQT%G
         imp_chain(3*(j-1)+1, i) = imp+1
         imp_chain(3*(j-1)+2, i) = imp+2
         imp_chain(3*(j-1)+3, i) = imp+3
         ichain_mp(imp+1) = i
         ichain_mp(imp+2) = i
         ichain_mp(imp+3) = i
         imp = imp + 3
      enddo
      write(*,'(a,141(i1))') '# ', seq(:,i)
      write(*,*) '# ', imp_chain(1,i), imp_chain((nrepeat-1)*3+3, i)
   enddo
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call list_local()
   call list_bp()
   !call read_sisinfo(cfile_sis)

   fdcd = file_dcd(hdl_dcd, cfile_dcd, DCD_OPEN_MODE%READ)

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
      write(*,*) "nmp = ", nmp, " is inconsistent with nmp_chain = ", nmp_chain, " and nchains = ", nchains
      stop (2) 
   endif

   allocate(xyz(3, nmp))

   open(hdl_out, file = cfile_out, status = 'replace', action = 'write', form='formatted')
   write(hdl_out, '(a)') '#(1)nframe  (2)Etotal  (3)Ebond   (4)Eangl   (5)Ebp   (6)Eele'

   nframe = 0
   do
      call fdcd%read_onestep(nmp, xyz, istat)
      if (istat == iostat_end) exit
      nframe = nframe + 1

      call energy()
      
      write(hdl_out, *) nframe, (energies(i), i=0,ENE%MAX)

      if (nframe == 1) then
         exit
      endif

   enddo

   write(*,*) '#nframe:', nframe

   call fdcd%close()

   close(hdl_out)
   if (flg_out_bp) then
      close(hdl_bp)
   endif

   deallocate(xyz)

   stop

end program sn
