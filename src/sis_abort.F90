subroutine sis_abort()

   use,intrinsic :: ISO_FORTRAN_ENV, only: output_unit

   use var_parallel

#ifdef PAR_MPI
   integer :: ierr
#endif

   flush(output_unit)

#ifdef PAR_MPI
   call MPI_ABORT(MPI_COMM_WORLD, 2, ierr)
   call MPI_FINALIZE(ierr)
#endif

   error stop

endsubroutine sis_abort
