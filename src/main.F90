program sis

   use, intrinsic :: iso_fortran_env, Only : iostat_end, compiler_version, compiler_options
   use const
   use const_phys, only : BOLTZ_KCAL_MOL
   use const_idx, only : ENE, SEQT, JOBT, seqt2char
   use var_top, only : nmp, nchains, nmp_chain, seq, imp_chain, ichain_mp, nrepeat, lmp_mp
   use var_state, only : restarted, xyz, tempK, kT, job, nthreads, rng_seed, opt_anneal
   use var_io, only : flg_out_bp, flg_out_bpall, flg_out_bpe, hdl_out, hdl_bp, hdl_bpall, hdl_bpe, KIND_OUT_BP, KIND_OUT_BPE, &
                      cfile_ff, cfile_dcd_in, cfile_prefix, cfile_out, cfile_fasta_in, hdl_rst
   use mt19937_64, only : init_genrand64
!$ use omp_lib

   implicit none

   character(len=CHAR_FILE_PATH) :: cfile_inp, cfile_bp, cfile_rst

   integer :: nargs
   integer :: i, j, k, imp
   integer :: istat

   character(len=500) :: cline
   logical :: stat

   call write_program_info()

   call init_const()

   nthreads = 1
!$  nthreads = omp_get_max_threads()
   if (nthreads > 1) then
      write(6, *) 'OpenMP nthreads = ', nthreads
   endif

   nargs = command_argument_count()

   if (nargs == 1) then

      call get_command_argument(1, cfile_inp)

      call read_input(cfile_inp, stat)

      if (.not. stat) then
         error stop 'Error in reading input file'
      endif
      
      restarted = .False.

   else if (nargs == 2) then

      call get_command_argument(1, cfile_inp)
      call get_command_argument(2, cfile_rst)

      call read_input(cfile_inp, stat)

      if (.not. stat) then
         error stop 'Error in reading input file'
      endif

      open(hdl_rst, file=cfile_rst, status='old', action='read', iostat=istat, form='unformatted', access='stream')

      if (istat > 0) then
         error stop 'Error: cannot open the restart file ' // trim(cfile_rst)
      endif

      restarted = .True.

   else if (nargs == 5) then

      restarted = .False.
      job = JOBT%DCD
      tempK = 273.15 + 22.0
      kT = BOLTZ_KCAL_MOL * tempK
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
      write(6,*) 'Usage: PROGRAM input.toml [restart fiel]'
      write(6,*) '  or : PROGRAM ff_file dcd_file nrepeat nchains out_prefix'
      stop (2) 
   end if

   !! Set RNG
   call init_genrand64(rng_seed)

   !! Load force field
   call read_force_field(stat)
   if (.not. stat) then
      write(6,*) 'Error in reading force field file'
      stop (2)
   endif

   if (opt_anneal > 0) then
      call read_anneal(stat)
      if (.not. stat) then
         write(6,*) 'Error in reading annealing-schedule file'
         stop (2)
      endif
   endif

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

      call read_fasta()

   else if (nrepeat > 0) then

      allocate(nmp_chain(nchains))
      nmp_chain(:) = 3 * nrepeat
      nmp = sum(nmp_chain)
      write(*, '(a,i5)') '#Nrepeat: ', nrepeat
      allocate(seq(3*nrepeat, nchains))
      allocate(imp_chain(3*nrepeat, nchains))
      allocate(ichain_mp(nmp))
      allocate(lmp_mp(nmp))
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
            lmp_mp(imp+1) = 3*(j-1)+1
            lmp_mp(imp+2) = 3*(j-1)+2
            lmp_mp(imp+3) = 3*(j-1)+3
            imp = imp + 3
         enddo
         !write(*,'(a,141(i1))') '# ', seq(:,i)
         !write(*,*) '# ', imp_chain(1,i), imp_chain((nrepeat-1)*3+3, i)
      enddo

   else
      write(6,*) 'Error: either FASTA or [repeat] is required.'
      stop (2)

   endif

   write(6, '(a)') '############ System ############'
   write(6, '(a,i8)') 'Nchain: ', nchains
   write(6, *) ''
   do i = 1, nchains
      write(6, '(a, i4)') 'Chain ', i
      write(6, '(a, i10)') 'Nnt: ', nmp_chain(i)
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
   enddo
   write(6, '(a)') '################################'
   write(6, *) ''

   write(6, '(a,f7.3)') '#T/K: ', tempK
   write(6, '(a,f7.5)') '#kT/kcal/mol: ', kT

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   call list_local()

   if (job == JOBT%MD) then
      ! Will use neighbor list
      continue

   else
      ! No neighbor list
      call list_bp()
      call list_exv()

   endif
   !call read_sisinfo(cfile_sis)

   if (job == JOBT%DCD) then
      write(6,*) 'Starting job_dcd'
      call job_dcd()

   else if (job == JOBT%CHECK_FORCE) then

      call job_check_force()

   else if (job == JOBT%MD) then

      write(6,*) 'Starting job_md'
      call job_md()

   endif

   if (flg_out_bp) then
      close(hdl_bp)
   endif
   if (flg_out_bpe) then
      close(hdl_bpe)
   endif

   deallocate(xyz)

   stop

contains

   subroutine write_program_info()

      character(len=40) :: githash, git
      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone

      call date_and_time(date, time, zone)
      git = githash()

      write(6, '(a)') '############ Program information ############'
      write(6, '(a)') 'SIS model simulation code by Naoto Hori'
      write(6, '(a)') 'Source: https://github.com/naotohori/sis'
      if (git(1:1) == '?') then
         write(6, '(a)') 'Version: 0.1'
      else
         write(6, '(a)') 'Git commit: ' // git
      endif
      write(6, '(a)') 'Compiler version: ' // compiler_version()
      write(6, '(a)') 'Compiler options: ' // compiler_options()
      write(6, '(a)') 'Executed at ' // date(1:4) // '-' // date(5:6) // '-' // date(7:8) // 'T' &
                      // time(1:2) // ':' // time(3:4) // ':' // time(5:6) // zone(1:3) // ':' // zone(4:5)  ! ISO 8601
      write(6, '(a)') '#############################################'
      write(6, *) ''

   end subroutine write_program_info

end program sis
