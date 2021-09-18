program sis

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   use const
   use const_phys, only : BOLTZ_KCAL_MOL
   use const_idx, only : ENE, SEQT, JOBT, seqt2char
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, pbc_box, pbc_box_half, flg_pbc, ichain_mp, nrepeat
   use var_state, only : xyz, tempK, kT, job
   use var_io, only : flg_out_bp, flg_out_bpall, flg_out_bpe, hdl_out, hdl_bp, hdl_bpall, hdl_bpe, KIND_OUT_BP, KIND_OUT_BPE, &
                      cfile_ff, cfile_dcd_in, cfile_prefix, cfile_out, cfile_fasta_in
!$ use omp_lib

   implicit none

   character(CHAR_FILE_PATH) :: cfile_inp, cfile_bp

   integer :: i, j, k, imp
   integer :: nthreads

   character(500) :: cline
   logical :: stat

#ifdef VERGIT
   character(14), parameter :: VERSION_DATE = VERDATE
   character(7),  parameter :: VERSION_BUILD = VERBUILD
   character(30), parameter :: VERSION_BRANCH= VERBRANCH

   write(*, '(13a,7a,9a,30a,14a,14a,6a)') '# Git commit ', VERSION_BUILD, ' (branch:', trim(VERSION_BRANCH), ') compiled on ', VERSION_DATE, ' (UTC)'
#endif

   nthreads = 1
!$  nthreads = omp_get_max_threads()
   if (nthreads > 1) then
      write(*, *) 'OpenMP nthreads = ', nthreads
   endif

   if (command_argument_count() == 1) then
      
      call get_command_argument(1, cfile_inp)  
      call read_input(cfile_inp, stat)

   else if (command_argument_count() == 5) then

      job = JOBT%DCD
      tempK = 273.15 + 22.0
      flg_out_bp = .False.
      flg_out_bpe = .True.

      call get_command_argument(1, cline)
      cfile_ff = trim(cline)
      call get_command_argument(2, cline)
      cfile_dcd_in = trim(cline)
      call get_command_argument(3, cline)
      read(cline, *) nrepeat 
      call get_command_argument(4, cline)
      read(cline, *) nchains
      call get_command_argument(5, cline)
      cfile_prefix = trim(cline)

   else
      write(6,*) 'Usage: PROGRAM input.toml'
      write(6,*) '  or : PROGRAM ff_file dcd_file nrepeat nchains out_prefix'
      stop (2) 
   end if

   !! Load force field
   write(*,*) 'Read force-field file: '//trim(cfile_ff)
   call read_force_field(trim(cfile_ff))

   !! Output files
   cfile_out = trim(cfile_prefix) // '.out'

   if (flg_out_bp) then
      cfile_bp = trim(cfile_prefix) // '.bp'
      open(hdl_bp, file=cfile_bp, status='replace', action='write', form='unformatted',access='stream')
      write(hdl_bp) int(KIND_OUT_BP,kind=4)
      write(hdl_bp) int(KIND_OUT_BPE,kind=4)
   endif

   if (flg_out_bpall) then
      cfile_bp = trim(cfile_prefix) // '.bpall'
      open(hdl_bpall, file=cfile_bp, status='replace', action='write', form='unformatted',access='stream')
      write(hdl_bpall) int(KIND_OUT_BP,kind=4)
      write(hdl_bpall) int(KIND_OUT_BPE,kind=4)
   endif

   if (flg_out_bpe) then
      cfile_bp = trim(cfile_prefix) // '.bpe'
      open(hdl_bpe, file=cfile_bp, status='replace', action='write', form='formatted')
   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (allocated(cfile_fasta_in)) then

      write(*,*) 'reading fasta'
      call read_fasta()

   else if (nrepeat > 0) then

      allocate(nmp_chain(nchains))
      nmp_chain(:) = 3 * nrepeat
      nmp = sum(nmp_chain)
      write(*, '(a,i5)') '#Nrepeat: ', nrepeat
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
         !write(*,'(a,141(i1))') '# ', seq(:,i)
         !write(*,*) '# ', imp_chain(1,i), imp_chain((nrepeat-1)*3+3, i)
      enddo

   else
      write(6,*) 'Error: either FASTA or [repeat] is required.'
      stop (2)

   endif

   write(6, '(a)') '########### System'
   write(6, '(a,i8)') '# Nchain: ', nchains
   write(6, '(a)') '#'
   do i = 1, nchains
      write(6, '(a, i4)') '# Chain ', i
      write(6, '(a, i10)') '# Nnt: ', nmp_chain(i)
      k = 0
      do j = 1, nmp_chain(i)
         write(6, '(a)', advance='no') seqt2char(seq(j,i))
         k = k + 1
         if (mod(k,100) == 0) then
            write(6, *) ''
            k = 0
         endif
      enddo
      if (k /= 0) then
         write(6, *) ''
      endif
      write(6, '(a)') '#'
   enddo

   kT = BOLTZ_KCAL_MOL * tempK
   write(6, '(a,f7.3)') '#T/K: ', tempK
   write(6, '(a,f7.5)') '#kT/kcal/mol: ', kT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call list_local()
   call list_bp()
   call list_exv()
   !call read_sisinfo(cfile_sis)

   if (job == JOBT%DCD) then
      write(6,*) 'Starting job_dcd'
      call job_dcd()

   else if (job == JOBT%CHECK_FORCE) then

      write(6,*) 'Starting job_check_force'
      call job_check_force()

   else if (job == JOBT%MD) then
      continue
      !write(6,*) 'Starting job_md'
      !call job_md()

   endif

   if (flg_out_bp) then
      close(hdl_bp)
   endif
   if (flg_out_bpe) then
      close(hdl_bpe)
   endif

   deallocate(xyz)

   stop

end program sis
