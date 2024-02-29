subroutine init_mts()

   ! Use Multiple stream Mersenne Twister PRNG  (mt_stream)
   ! http://theo.phys.sci.hiroshima-u.ac.jp/~ishikawa/PRNG/mt_stream_en.html
   !use mt_kind_defs
   use mt_stream, only : set_mt19937, new, init, create_stream
   
   use const_idx, only : RSTBLK
   use var_state, only : mts, mts_rep, rng_seed, restarted
   use var_replica, only : irep2grep, nrep_proc, flg_replica

   implicit none

   integer :: irep, grep
   integer :: id_offset
   integer :: rst_status

   print '(a)', 'Initialising MT stream.'
   flush(6)

   allocate(mts(0:nrep_proc))

   call set_mt19937()
   call new(mts(0))
   call init(mts(0), rng_seed)

   id_offset = 0

   ! If REMD, geenerate another stream specifically for replica exchange.
   if (flg_replica) then
      print '(a)', 'Creating an MT stream, mts_rep.'
      call create_stream(mts(0), mts_rep, 1)
      id_offset = 1
   endif     
  
   ! Stream for each process.
   do irep = 1, nrep_proc
      grep = irep2grep(irep)
      if (grep /= 0) then
         print '(a,i4,a,i4,a,i5)', 'Creating an MT stream, mts(', irep, '), grep = ', grep, ', id = ', grep + id_offset
         call create_stream(mts(0), mts(irep), grep + id_offset)
      endif
   enddo
 
   print '(a)', 'Done: Initialising MT stream.'
   print *
   flush(6)

   if (restarted) then
      print '(a)', 'Loading MT streams from the restart file.'

      if (flg_replica) call read_rst(RSTBLK%PRNGREP, rst_status)

      if (rst_status /= 0) then
         print '(a)', '... Failed to load PRNGREP from the restart file. MT state from the given seed will be used.'
      endif

      call read_rst(RSTBLK%PRNG, rst_status)

      if (rst_status /= 0) then
         print '(a)', '... Failed to load PRNG from the restart file. MT states from the given seed will be used.'
      endif

      print '(a)', 'Done: Loading MT streams from the restart file.'
      print *
      flush(6)
   endif

endsubroutine init_mts
