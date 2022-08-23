subroutine init_replica

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const_idx, only : REPT
   use var_replica, only : flg_replica, nrep_all, rep2val, nrep_proc, irep2grep, &
                           n_replica_temp, rep2lab, lab2rep, lab2val, replica_values

   implicit none

   integer :: irep

   if (flg_replica) then

      nrep_all = n_replica_temp
      nrep_proc = n_replica_temp

      do irep = 1, nrep_all
         irep2grep(irep) = irep
         rep2lab(irep) = irep
         lab2rep(irep) = irep
         lab2val(irep, REPT%TEMP) = replica_values(irep, REPT%TEMP)
         print *, 'replica', irep, ', temp = ', lab2val(irep, REPT%TEMP)
      enddo

   else

      nrep_all = 1
      nrep_proc = 1
      irep2grep(:) = 1

   endif

   flush(output_unit)


endsubroutine init_replica
