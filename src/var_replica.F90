module var_replica

   use,intrinsic :: ISO_FORTRAN_ENV, only: INT64
   use const, only : PREC, MAX_REPLICA, MAX_REP_DIM, MAX_REP_PER_DIM
   use const_idx, only : REPT
   implicit none

   logical :: flg_replica
   logical :: flg_repvar(REPT%MAX)

   integer :: ndim_replica

   integer :: nrep_all
   integer :: nrep_proc

   integer :: nrep(REPT%MAX)

   integer(INT64) :: nstep_rep_exchange
   integer(INT64) :: nstep_rep_save

   logical :: flg_exchange

   real(PREC) :: replica_values(MAX_REP_PER_DIM, REPT%MAX)

   integer, save :: irep2grep(MAX_REPLICA)  !< replica(local) => replica(global)
   integer, save :: grep2irep(MAX_REPLICA)  !< replica(global) => replica(local)
   integer, save :: grep2rank(MAX_REPLICA)  !< replica(global) => local_rank_rep
   !logical, save :: flg_rep(1:REPTYPE%MAX) = .false.


   ! ================================================
   ! Example of parallelization of replica, 
   !    #node=3, #replica=6
   !    -----------------------------------------
   !     Node | replica(global) | replica(local)
   !          |     irep        |    grep
   !    -----------------------------------------
   !      1   |       1         |      1
   !          |       2         |      2
   !      2   |       3         |      1
   !          |       4         |      2
   !      3   |       5         |      1
   !          |       6         |      2
   !    -----------------------------------------
   !     max  |    nrep_all     | n_replica_mpi
   ! ================================================

   ! Permutation function
   integer,    save :: rep2lab(MAX_REPLICA)  !< Permutation function (replica => label)
   integer,    save :: lab2rep(MAX_REPLICA)  !< Permutation function ( label  => replica)
   real(PREC), save :: lab2val(MAX_REPLICA, REPT%MAX) !< Replica variables

!   integer, allocatable, save :: label_order(:,:)  !< (REPTYPE%MAX, n_repica_all) 

   ! ===============================================================
   ! Using these 'Permutation function'
   ! Ex.)
   ! replica => label  (rep2lab)
   ! label => replica  (lab2rep)
   ! --------------------------------
   ! | replica  | 1 | 2 | 3 | 4 | 5 |
   ! | label    | 3 | 1 | 5 | 2 | 4 |
   ! --------------------------------
   ! 
   ! label   => value  (lab2val)
   ! --------------------------------
   ! | label    | 1 | 2 | 3 | 4 | 5 |
   ! | value    |1.0|1.2|1.4|1.6|1.8| (temperature, ionic strength, and so on.)
   ! --------------------------------
   !
   ! As a result, in this situation, 
   ! replica => value  (rep2val (function))
   ! --------------------------------
   ! | replica  | 1 | 2 | 3 | 4 | 5 |
   ! | value    |1.4|1.0|1.8|1.2|1.6| (temperature, ionic strength, and so on.)
   ! --------------------------------
   ! ===============================================================

   !-----------------------------------------------------------
   ! to detect whether current exchange step is odd or even
   ! (These are private variables)
   integer, save, private :: type_array(1:MAX_REP_DIM)
   integer, save, private :: exchange_pair_tb(1:MAX_REPLICA,1:MAX_REP_DIM*2)
   integer, save, private :: idx_type = 1
   integer, save, private :: idx_pair = 1

!! ###########################################################################
contains

   real(PREC) function rep2val (ireplica, vartype)
      implicit none
      integer :: ireplica
      integer :: vartype
      rep2val = lab2val(rep2lab(ireplica), vartype)
   endfunction rep2val


!   integer(L_INT) function rep2step (ireplica)
!      implicit none
!      integer :: ireplica
!      rep2step = lab2step(rep2lab(ireplica))
!   endfunction rep2step
!

   integer function get_pair(i)
      implicit none
      integer, intent(in) :: i
      get_pair = exchange_pair_tb(i, idx_pair)
   endfunction get_pair


   subroutine set_forward()
      implicit none

      if (idx_type == ndim_replica) then
         idx_type = 1
      else
         idx_type = idx_type + 1
      endif

      if (idx_pair == ndim_replica*2) then
         idx_pair = 1
      else
         idx_pair = idx_pair + 1
      endif
   endsubroutine set_forward


   integer function get_type(i)
      implicit none
      integer, intent(in), optional :: i
      if (present(i)) then
         get_type = type_array(i)
      else
         get_type = type_array(idx_type)
      endif
   endfunction get_type


   subroutine make_replica_type_array()
      ! called by replica_settable only once
      use const, only : MAX_REP_DIM
      use const_idx, only : REPT
      implicit none
      integer :: idx, ivar

      type_array(:) = 0
      idx = 0
      do ivar = 1, REPT%MAX
         if (flg_repvar(ivar)) then
            idx = idx + 1
            if (idx > MAX_REP_DIM) then
               print '(a)', 'Erorr: defect in var_replica::make_type_table'
               call sis_abort()
            endif
            type_array(idx) = ivar
         endif
      enddo
   endsubroutine make_replica_type_array


   subroutine make_replica_exchange_pair_tb()
      implicit none
      integer :: icounter, ivar, iset
      integer :: idimn, iparity, icycle, ireplica, icontinue
      integer :: ireplica_start
      integer :: n_division_pre, n_division, n_continue
      logical :: flg_post_zero
#ifdef _DEBUG
      integer :: irep
#endif

      exchange_pair_tb(:,:) = 0

      ! ==========================================================
      ! exchange order  (Ex. 3-dimension case)
      !   -----------------------------------------------------
      !    icounter  |  1   |  2   |  3   |  4   |  5   |  6
      !   -----------------------------------------------------
      !    exchange- | dim1 | dim2 | dim3 | dim1 | dim2 | dim3 
      !      pair    | odd  | odd  | odd  | even | even | even
      !   -----------------------------------------------------
      ! ==========================================================

      icounter = 0
      do iparity = 1, 2   ! odd or even
                          ! 1 odd  : exchange 1-2, 3-4, 5-6 ,.....
                          ! 2 even : exchange 2-3, 4-5, 6-7 ,.....

         n_division     = 1
         do idimn = 1, ndim_replica 

            icounter = icounter + 1
            ivar = type_array(idimn)
            n_division_pre = n_division
            n_division     = n_division * nrep(ivar) 
            n_continue     = nrep_all / n_division

            iset = 0
            do icycle = 1, n_division_pre

               flg_post_zero = .false.
               if (iparity == 1) then
                  if (mod(nrep(ivar),2) == 1) then 
                     ! #replica is odd
                     flg_post_zero = .true.
                  endif
                  ireplica_start = 1
               else
                  ! zeroing 1st replica
                  do icontinue = 1, n_continue
                     iset = iset + 1
                     exchange_pair_tb(iset, icounter) = 0
                  enddo
                  if (mod(nrep(ivar),2) == 0) then
                     ! #replica is even
                     flg_post_zero = .true.
                   endif
                  ireplica_start = 2
               endif

               do ireplica = ireplica_start, (nrep(ivar)-1), 2

                  do icontinue = 1, n_continue*2
                     iset = iset + 1

                     if (icontinue <= n_continue) then
                        exchange_pair_tb(iset, icounter) = iset + n_continue
                     else
                        exchange_pair_tb(iset, icounter) = iset - n_continue
                     endif
                  enddo

               enddo ! ireplica

               ! zeroing last replica
               if (flg_post_zero) then
                  do icontinue = 1, n_continue
                     iset = iset + 1
                     exchange_pair_tb(iset, icounter) = 0
                  enddo
               endif

            enddo ! icycle
         enddo ! iparity
      enddo ! idimn

#ifdef _DEBUG
      do icounter = 1, MAX_REP_DIM*2
         print  '(a)', '#################'
         print '(a,i5)', 'icounter = ',icounter
         do irep = 1, nrep_all
            print *, 'exchange_pair_tb(',irep,',',icounter,')=', exchange_pair_tb(irep, icounter)
         enddo
      enddo
      flush(6)
#endif

   endsubroutine make_replica_exchange_pair_tb

endmodule var_replica
