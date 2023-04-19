!#ifdef _DUMP_REP_EX
!#define DUMP
!#else
!#undef DUMP
!#endif

#define DUMP

subroutine replica_exchange(velos, replica_energy, tempk)

   use const
   use const_phys, only : BOLTZ_KCAL_MOL
   use const_idx, only : REPT
   use var_state, only : mts_rep
   use var_replica, only : nrep_all, &
                           rep2lab, lab2rep, lab2val,&
                           flg_repvar, get_pair, set_forward
   use mt_stream
   use var_parallel

   implicit none

   !---------------------------------------------------------------------------
   real(PREC), intent(inout) :: velos(:,:,:)      ! (SDIM, mp, replica)
   real(PREC), intent(in)    :: replica_energy(:,:) ! (2, replica)
   real(PREC), intent(in)    :: tempk

   !---------------------------------------------------------------------------
   ! local
   real(PREC) :: delta  ! energy difference
   real(8)    :: random0(nrep_all) ! for Metropolis
   integer    :: nrand

   ! exchanged targets are represented by "i" and "j"
   integer    :: i
   !integer    :: ivar
   integer    :: rep_i, rep_j
   integer    :: l_i, l_j
   real(PREC) :: temp_m, temp_n
   real(PREC) :: pot_i_m, pot_i_n
   real(PREC) :: pot_j_m, pot_j_n

   !integer :: ntotal_exchange


   ! ###########################################
   !  Generate random numbers
   ! ###########################################
   nrand = 0
   do i = 1, nrep_all, 2
      nrand = nrand + 1
      random0(nrand) = genrand_double1(mts_rep)
   enddo
!#ifdef MPI_REP
!   call MPI_BCAST(random0, nrand, MPI_REAL8, 0, mpi_comm_rep, ierr)
!#endif


   ! ###########################################
   !  Main operation (judgement and exchange)
   ! ###########################################
   !  l_i   : label i
   !  l_j   : label j
   !  rep_i : replica i  (a replica which has label i)
   !  rep_j : replica j  (a replica which has label j)
   ! 
   ! PRE-EXCHANGE
   !           label = temperature & potential
   !   rep_i :  l_i  =   temp_m    &  pot_i_m
   !   rep_j :  l_j  =   temp_n    &  pot_j_n
   !
   ! CANDIDATE
   !           label = temperature & potential
   !   rep_i :  l_j  =   temp_n    &  pot_i_n
   !   rep_j :  l_i  =   temp_m    &  pot_j_m
   !
#ifdef DUMP
   write(6,*) ''
#endif
   nrand = 0
   do l_i = 1, nrep_all

      l_j = get_pair(l_i)

      if (l_j < l_i) then
         cycle
      endif

      nrand = nrand + 1

      rep_i     = lab2rep( l_i)
      rep_j     = lab2rep( l_j)
#ifdef DUMP
      write(6,*) 'exchange label: ', l_i, l_j
      write(6,*) 'exchange replica: ', rep_i, rep_j
#endif

      if (flg_repvar(REPT%TEMP)) then
         temp_m = lab2val( l_i ,REPT%TEMP) 
         temp_n = lab2val( l_j ,REPT%TEMP)
         ! Note that, temp_m < temp_n
         !   i.e.      T(rep_i) < T(rep_j)
      else
         temp_m = tempk
         temp_n = tempk
      endif
#ifdef DUMP
      write(6,*) 'temp_m, temp_n = ', temp_m, temp_n
#endif

      pot_i_m = replica_energy(1, rep_i)
      pot_i_n = replica_energy(2, rep_i)
      pot_j_n = replica_energy(1, rep_j)
      pot_j_m = replica_energy(2, rep_j)

      ! Y. Sugita et al. J Chem phys (2000)  eq.[14]
      delta = ( (pot_j_m - pot_i_m) / temp_m    &
               -(pot_j_n - pot_i_n) / temp_n )  &
             / BOLTZ_KCAL_MOL
#ifdef DUMP
      write(6,*) 'pot_i_m, pot_i_n = ',pot_i_m,pot_i_n
      write(6,*) 'pot_j_m, pot_j_n = ',pot_j_m,pot_j_n
      write(6,*) 'delta = ',delta
#endif

      ! Judgement
      if (delta < 0.0e0_PREC) then
         call exchange_permutation()
#ifdef DUMP
         write(6,*) 'Exchange (delta < 0.0)'
#endif
         if (flg_repvar(REPT%TEMP)) then
            call scale_velo(temp_m, temp_n)
         endif

      else
         if (random0(nrand) <= exp(-delta)) then
            call exchange_permutation()
            if (flg_repvar(REPT%TEMP)) then
               call scale_velo(temp_m, temp_n)
            endif
#ifdef DUMP
         write(6,*) 'Exchange (random(',random0(nrand),') <= delta'
#endif
         else
#ifdef DUMP
         write(6,*) 'NOT Exchange (random(',random0(nrand),') > delta'
#endif
         endif

      endif

   enddo ! l_i 


   ! ###########################################
   !  Aftertreatment
   ! ###########################################
   call set_forward()



!---------------------------------------------------------------------------
contains

   !---------------------------------!
   !-- exchange permutation arrays --!
   !---------------------------------!
   subroutine exchange_permutation()
      rep2lab(rep_i) = l_j
      rep2lab(rep_j) = l_i
   
      lab2rep(l_i) = rep_j
      lab2rep(l_j) = rep_i
   end subroutine exchange_permutation

   !---------------------------------!
   !-- scale velocity ---------------!
   !---------------------------------!
   subroutine scale_velo(temperature_i, temperature_j)
      use var_replica, only : irep2grep, nrep_proc
      implicit none
      real(PREC), intent(in) :: temperature_i, temperature_j
      real(PREC) :: scale_factor
      integer    :: irep, grep

      scale_factor = sqrt( temperature_j / temperature_i )
#ifdef DUMP
      write(6,*) '##scale velo'
      write(6,*) 'Ti = ',temperature_i
      write(6,*) 'Tj = ',temperature_j
      write(6,*) 'Tj / Ti = ', temperature_j / temperature_i 
      write(6,*) 'SQRT(Tj/Ti) = ', sqrt( temperature_j / temperature_i )
#endif

      do irep = 1, nrep_proc
         grep = irep2grep(irep)

         if (rep_i == grep) then
#ifdef DUMP
            write(6,*) 'scale_velo i pre : ',velos(1,1,irep)
#endif
            velos(:,:, irep) = scale_factor * velos(:,:, irep)
#ifdef DUMP
            write(6,*) 'scale_velo i post: ',velos(1,1,irep)
#endif
         endif
         if (rep_j == grep) then
#ifdef DUMP
            write(6,*) 'scale_velo j pre : ',velos(1,1,irep)
#endif
            velos(:,:, irep) = velos(:,:, irep) / scale_factor
#ifdef DUMP
            write(6,*) 'scale_velo j post: ',velos(1,1,irep)
#endif
         endif
      enddo
   end subroutine

end subroutine replica_exchange
#undef DUMP
