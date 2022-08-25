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

      call mpi_init(ierr)
      call mpi_comm_size(MPI_COMM_WORLD, nprocs, ierr)
      call mpi_comm_rank(MPI_COMM_WORLD, myrank, ierr)
#else
      myrank = 0
      nprocs = 1
#endif

      nthreads = 1
!$    nthreads = omp_get_max_threads()

      ncores = nprocs * nthreads

   end subroutine init_parallel


   subroutine sis_abort()
      use,intrinsic :: ISO_FORTRAN_ENV, only: output_unit
      integer :: ierr

      flush(output_unit)
#ifdef PAR_MPI
      call MPI_ABORT(MPI_COMM_WORLD, 2, ierr)
      call MPI_FINALIZE(ierr)
#endif
      error stop
   endsubroutine sis_abort

end module var_parallel
