program sis

   use, intrinsic :: iso_fortran_env, only : iostat_end, compiler_version, compiler_options, output_unit
   use const, only : CHAR_FILE_PATH, init_const, PREC, FILENAME_DIGIT_REPLICA
   use const_phys, only : BOLTZ_KCAL_MOL
   use const_idx, only : ENE, JOBT, REPT, RSTBLK
   use var_potential, only : flg_ele
   use var_state, only : restarted, xyz, tempK, kT, job, rng_seed, opt_anneal, anneal_tempK
   use var_top, only : nmp
   use var_io, only : flg_out_bp, flg_out_bpall, flg_out_bpe, hdl_in_rst, &
                      hdl_out, hdl_bp, hdl_bpall, hdl_bpe
   use var_parallel, only : init_parallel, end_parallel
   use var_replica, only : nrep_proc
   use mt19937_64, only : init_genrand64

   implicit none

   character(len=CHAR_FILE_PATH) :: cfile_inp, cfile_rst

   integer :: nargs
   integer :: istat
   integer :: irep
   logical :: stat

   call init_const()
   
   call init_parallel()

   call print_program_info()

   nargs = command_argument_count()

   if (nargs < 1 .or. 2 < nargs) then
      print *, 'Usage: PROGRAM input.toml [restart file (.rst)]'
      stop
   endif


   !! Read input file
   call get_command_argument(1, cfile_inp)
   call read_input(cfile_inp)


   !! Initialise replicas
   call init_replica()


   !! Open restart file if given
   restarted = .False.

   if (nargs == 2) then
      call get_command_argument(2, cfile_rst)

      open(hdl_in_rst, file=cfile_rst, status='old', action='read', iostat=istat, form='unformatted', access='stream')

      if (istat > 0) then
         print '(2a)', 'Error: cannot open the restart file ', trim(cfile_rst)
         error stop
      endif

      restarted = .True.
   end if


   !! Set RNG
   call init_genrand64(rng_seed)

   !! Load force field
   call read_force_field()


   !! Read annealing file
   if (opt_anneal > 0) then
      call read_anneal(stat)
      if (.not. stat) then
         print *, 'Error in reading annealing-schedule file'
         stop (2)
      endif

      tempK = anneal_tempK(1)
   endif

   
   !! Prepare output files
   call init_out_files()


   !! Construct the sequences
   call init_sequence()

   !! Read coordinates (xyz)
   if (job /= JOBT%DCD) then

      allocate(xyz(3, nmp, nrep_proc))

      if (restarted) then
         call read_rst(RSTBLK%XYZ)
      else
         call init_structure()
      endif
   endif

   print '(a)', 'Temperature'
   print '(a,f7.3)', '# T(K): ', tempK
   print '(a,f7.5)', '# kT(kcal/mol): ', kT
   print *
   flush(output_unit)

   !! Construct possible combinations of basepairs
   call init_bp()

   !! Allocation and initialisation of electrostatics
   if (flg_ele) call init_ele()

   !! Construct pair lists of local potentials
   call list_local()

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   if (job == JOBT%MD) then
      ! Will use neighbor list
      continue

   else
      ! No neighbor list
      call list_bp()
      call list_exv()
      if (flg_ele) call list_ele()

   endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Main jobs
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (job == JOBT%DCD) then
      print '(a)', 'Starting DCD'
      call job_dcd()

   else if (job == JOBT%CHECK_FORCE) then

      call job_check_force()

   else if (job == JOBT%MD) then

      print '(a)', 'Starting MD'
      call job_md()

   endif

   do irep = 1, nrep_proc
      if (flg_out_bp) close(hdl_bp(irep))
      if (flg_out_bpe) close(hdl_bpe(irep))
      if (flg_out_bpall) close(hdl_bpall(irep))
   enddo

   deallocate(xyz)

   do irep = 1, nrep_proc
      close(hdl_out(irep))
   enddo
   deallocate(hdl_out)

   flush(output_unit)

   call end_parallel()

   stop

contains

   subroutine print_program_info()

      use var_parallel, only : myrank, nprocs, nthreads, ncores

      character(len=40) :: githash, git
      character(len=8) :: date
      character(len=10) :: time
      character(len=5) :: zone
      character(len=500) :: com
      integer :: length

      call date_and_time(date, time, zone)
      call get_command(com, length)
      git = githash()

      print '(a)', '############ Program information ############'
      print '(a)', 'SIS model simulation code by Naoto Hori'
      print '(a)', 'Source: https://github.com/naotohori/sis'
      if (git(1:1) == '?') then
         print '(a)', 'Version: 0.2'
      else
         print '(2a)', 'Git commit: ', git
      endif
      print '(2a)', 'Compiler version: ', compiler_version()
      print '(2a)', 'Compiler options: ', compiler_options()
      print '(2a)', 'Command: ', com(1:length)
      print '(15a)', 'Executed at ', date(1:4), '-', date(5:6), '-', date(7:8), &  ! Date
                     'T', time(1:2), ':', time(3:4), ':', time(5:6), zone(1:3), ':', zone(4:5)  ! ISO 8601 Time
      print '(a,i6)', 'Number of cores: ', ncores
      print '(a,i6)', 'Number of OpenMP threads: ', nthreads
      print '(a,i6)', 'Number of MPI processes: ', nprocs
      if (nprocs > 1) then
         print '(a,i6)', 'MPI rank: ', myrank
      endif
      print '(a)', '#############################################'
      print *

   end subroutine print_program_info

end program sis
