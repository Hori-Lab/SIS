subroutine init_replica

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const_idx, only : REPT
   use var_replica, only : flg_replica, nrep_all, rep2val, nrep_proc, irep2grep, &
                           n_replica_temp, rep2lab, lab2rep, lab2val, replica_values
   use var_parallel, only : myrank, nprocs

   implicit none

   integer :: irep, grep

   if (flg_replica) then

      print '(a)', 'Initializing replicas'

      nrep_all = n_replica_temp

#ifdef PAR_MPI
      if (nprocs > n_replica_temp) then
         print *, 'Error: the number of MPI processes should be equal to or smaller than the number of replicas.'
         call sis_abort()
      endif

      if (mod(n_replica_temp, nprocs) /= 0) then
         print *, 'Error: the number of replicas has to be a multiple of the number of MPI processes.'
         call sis_abort()
      endif

      nrep_proc = n_replica_temp / nprocs
#else
      nrep_proc = n_replica_temp
#endif

      do grep = 1, nrep_all
         rep2lab(grep) = grep
         lab2rep(grep) = grep
         lab2val(grep, REPT%TEMP) = replica_values(grep, REPT%TEMP)
         print '(a,i5,a,f8.3)', '# Replica', grep, ', temp = ', lab2val(grep, REPT%TEMP)
      enddo
      print '(a)', '#'

      print '(a)', '# This process is responsible for the following replicas.'
      do irep = 1, nrep_proc
         grep = myrank * nrep_proc + irep
         irep2grep(irep) = grep
         print '(a,i5,a,i5)', '# irep = ', irep, ' ==> grep = ', grep
      enddo
      print '(a)', '#'
   else

      nrep_all = 1
      nrep_proc = 1
      irep2grep(:) = 1

   endif

   print *
   flush(output_unit)


endsubroutine init_replica
