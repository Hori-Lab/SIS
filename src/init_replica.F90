subroutine init_replica

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const, only : MAX_REPLICA
   use const_idx, only : REPT, RSTBLK
   use const_phys, only : INVALID_VALUE, JOUL2KCAL_MOL
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
   integer :: rst_status
   integer :: l_grep2irep(MAX_REPLICA) !! See comments for MPI_ALLREDUCE below.
   integer :: l_grep2rank(MAX_REPLICA)
#ifdef PAR_MPI
   integer :: istat
#endif

   ndim_replica = 0
   nrep_all = 1
   nrep_proc = 1
   irep2grep(:) = 1
   lab2val(:,:) = INVALID_VALUE

   if (flg_replica) then

      ! The followint two have to be zero in order to use MPI_ALLREDUCE below.
      !grep2irep(:) = 0
      !grep2rank(:) = 0
      l_grep2irep(:) = 0
      l_grep2rank(:) = 0

      print '(a)', 'Initializing replicas'

      do ivar = 1, REPT%MAX
         if (flg_repvar(ivar)) then
            ndim_replica = ndim_replica + 1
            nrep_all = nrep_all * nrep(ivar)
         endif
      enddo

      if (nprocs > nrep_all) then
         print *, 'Error: the number of MPI processes should be equal to or smaller than the number of replicas.'
         call sis_abort()
      endif

#ifdef PAR_MPI
      if (mod(nrep_all, nprocs) /= 0) then
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

         if (flg_repvar(REPT%TEMP)) then 
            lab2val(grep, REPT%TEMP) = replica_values(grep, REPT%TEMP)
            print '(a,i5,a,f8.3)', '# Replica ', grep, ', temp = ', lab2val(grep, REPT%TEMP)
         endif

         if (flg_repvar(REPT%TWZDCF)) then 
            lab2val(grep, REPT%TWZDCF) = replica_values(grep, REPT%TWZDCF)
            print '(a,i5,a,f8.3)', '# Replica ', grep, ', force = ', lab2val(grep, REPT%TWZDCF) / (JOUL2KCAL_MOL * 1.0e-22)
         endif
      enddo
      print '(a)', '#'

      print '(a)', '# This process is responsible for the following replicas.'
      do irep = 1, nrep_proc
         grep = myrank * nrep_proc + irep
         irep2grep(irep) = grep
         !grep2irep(grep) = irep
         !grep2rank(grep) = myrank
         l_grep2irep(grep) = irep
         l_grep2rank(grep) = myrank
         print '(a,i5,a,i5)', '# irep = ', irep, ' ==> grep = ', grep
      enddo
      print '(a)', '#'
      flush(output_unit)

#ifdef PAR_MPI
      !call MPI_ALLREDUCE(MPI_IN_PLACE, grep2irep, MAX_REPLICA, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, istat)
      !call MPI_ALLREDUCE(MPI_IN_PLACE, grep2rank, MAX_REPLICA, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, istat)

      !! Somehow MPI_IN_PLACE does not work properly and all elements become zero after this. (only issue in Mac?). 
      !! For now, the workaround is just use temporal variable l_grep2rank and l_grep2irep.
      call MPI_ALLREDUCE(l_grep2irep, grep2irep, MAX_REPLICA, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, istat)
      call MPI_ALLREDUCE(l_grep2rank, grep2rank, MAX_REPLICA, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, istat)

      do irep = 1, nrep_all
         print '(a,i5,a,i5)', "# MPI grep2irep(", irep, ") = ", grep2irep(irep)
      enddo
      do irep = 1, nrep_all
         print '(a,i5,a,i5)', "# MPI grep2rank(", irep, ") = ", grep2rank(irep)
      enddo
      print '(a)', '#'
#endif

      if (restarted) then
         call read_rst(RSTBLK%REPLICA, rst_status)
         do grep = 1, nrep_all
            print '(a,i5,a,i5)', '# Replica ', grep, ', label = ', rep2lab(grep)
         enddo
         print '(a)', '#'
      endif

      call make_replica_type_array()

      call make_replica_exchange_pair_tb()

   else ! flg_replica

      grep2irep(:) = 1
      grep2rank(:) = 0
   endif

   print *
   flush(output_unit)


endsubroutine init_replica
