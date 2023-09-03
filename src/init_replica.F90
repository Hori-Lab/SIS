subroutine init_replica

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const_idx, only : REPT, MAX_REPLICA, RSTBLK
   use var_state, only : restarted
   use var_replica, only : flg_replica, nrep_all, rep2val, nrep_proc, &
                           irep2grep, grep2irep, grep2rank, &
                           nrep, rep2lab, lab2rep, lab2val, replica_values, &
                           ndim_replica, flg_repvar, &
                           make_replica_exchange_pair_tb, make_replica_type_array
   use var_parallel

   implicit none

   integer :: irep, grep
   integer :: ivar
#ifdef PAR_MPI
   integer :: istat
#endif

   ndim_replica = 0
   nrep_all = 1
   nrep_proc = 1
   irep2grep(:) = 1
   grep2irep(:) = 1
   grep2rank(:) = 0

   if (flg_replica) then

      print '(a)', 'Initializing replicas'

      do ivar = 1, REPT%MAX
         if (flg_repvar(ivar)) then
            ndim_replica = ndim_replica + 1
            nrep_all = nrep_all * nrep(ivar)
         endif
      enddo

      print *, 'nrep_all=',nrep_all
      print *, 'nprocs=',nprocs
      print *, 'nrep(REPT%TEMP)=', nrep(REPT%TEMP)

      if (nprocs > nrep(REPT%TEMP)) then
         print *, 'Error: the number of MPI processes should be equal to or smaller than the number of replicas.'
         call sis_abort()
      endif

#ifdef PAR_MPI
      if (mod(nrep(REPT%TEMP), nprocs) /= 0) then
         print *, 'Error: the number of replicas has to be a multiple of the number of MPI processes.'
         call sis_abort()
      endif

      nrep_proc = nrep_all / nprocs
      print *, 'nrep_proc=',nrep_proc
#else
      nrep_proc = nrep_all
#endif

      do grep = 1, nrep_all
         rep2lab(grep) = grep
         lab2rep(grep) = grep
         lab2val(grep, REPT%TEMP) = replica_values(grep, REPT%TEMP)
         print '(a,i5,a,f8.3)', '# Replica ', grep, ', temp = ', lab2val(grep, REPT%TEMP)
      enddo
      print '(a)', '#'

      print '(a)', '# This process is responsible for the following replicas.'
      do irep = 1, nrep_proc
         grep = myrank * nrep_proc + irep
         irep2grep(irep) = grep
         grep2irep(grep) = irep
         grep2rank(grep) = myrank
         print '(a,i5,a,i5)', '# irep = ', irep, ' ==> grep = ', grep
      enddo
      print '(a)', '#'

#ifdef PAR_MPI
      call MPI_ALLREDUCE(MPI_IN_PLACE, grep2irep, MAX_REPLICA, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, istat)
      call MPI_ALLREDUCE(MPI_IN_PLACE, grep2rank, MAX_REPLICA, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, istat)

      do irep = 1, nrep_all
         print '(a,i5,a,i5)', "# MPI grep2irep(", irep, ") = ", grep2irep(irep)
      enddo
      do irep = 1, nrep_all
         print '(a,i5,a,i5)', "# MPI irep2rank(", irep, ") = ", grep2rank(irep)
      enddo
      print '(a)', '#'
#endif

      if (restarted) then
         call read_rst(RSTBLK%REPLICA)
         do grep = 1, nrep_all
            print '(a,i5,a,i5)', '# Replica ', grep, ', label = ', rep2lab(grep)
         enddo
         print '(a)', '#'
      endif

      call make_replica_type_array()

      call make_replica_exchange_pair_tb()
   endif

   print *
   flush(output_unit)


endsubroutine init_replica
