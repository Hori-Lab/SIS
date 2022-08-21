module var_parallel
   
#ifdef PAR_MPI
   use mpi
#endif
   
!$ use omp_lib

   implicit none

   integer :: myrank
   integer :: nprocs
   integer :: nthreads
   integer :: ncores

contains

   subroutine init_parallel()

#ifdef PAR_MPI
      integer :: ierr
#endif

#ifdef PAR_MPI
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

end module var_parallel
