module var_parallel
   
#ifdef PAR_MPI
   use mpi
#endif
   
   use const, only : PREC

!$ use omp_lib

   implicit none

   integer :: myrank
   integer :: nprocs
   integer :: nthreads
   integer :: ncores
   integer :: PREC_MPI

contains

   subroutine init_parallel()

#ifdef PAR_MPI
      integer :: ierr

      if (PREC == 4) then
         PREC_MPI = MPI_REAL
      else if (PREC == 8) then
         PREC_MPI = MPI_DOUBLE_PRECISION
      else if (PREC == 16) then
         PREC_MPI = MPI_REAL16
      else
         print *, 'Error: could not find suitable PREC_MPI.'
         call sis_abort()
      endif

      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)
#else
      myrank = 0
      nprocs = 1
#endif

      nthreads = 1
!$    nthreads = omp_get_max_threads()

      ncores = nprocs * nthreads

   end subroutine init_parallel


   subroutine end_parallel()
#ifdef PAR_MPI
      integer :: ierr
#endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_FINALIZE(ierr)
   endsubroutine end_parallel

end module var_parallel
