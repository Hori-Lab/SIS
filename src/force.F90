subroutine force()

!$ use omp_lib
   use const
   use const_idx, only : ENE
   use var_state, only : nthreads, forces
   use var_potential, only : flg_ele, max_bp_per_nt, flg_dih
   use var_top, only : nmp

   implicit none

   integer :: tn
   !logical, save :: flg_first = .True.
   !real(PREC), save, allocatable :: forces_t(:,:,:)
   real(PREC) :: forces_t(3, nmp, 0:nthreads-1)

   !if (flg_first) then
   !   !allocate(forces_t(3, nmp, 0:nthreads-1))
   !   flg_first = .False.
   !endif

!$omp parallel private(tn)
   tn = 0
!$ tn = omp_get_thread_num()
   forces_t(:,:,tn) = 0.0e0_PREC

   call force_bond(forces_t(1,1,tn))

   call force_angl(forces_t(1,1,tn))
   
   if (flg_dih) call force_dihedral(forces_t(1, 1, tn))
      
   
   if (max_bp_per_nt < 1) then
      call force_bp(forces_t(1,1,tn))
   else
      call force_bp_limit(forces_t(1,1,tn))
   endif

   call force_wca(forces_t(1,1,tn))

   if (flg_ele) call force_ele_DH(forces_t(1,1,tn))

!$omp end parallel

   forces(:,:) = forces_t(:,:,0)
   do tn = 1, nthreads-1
      forces(:,:) = forces(:,:) + forces_t(:,:,tn)
   enddo

end subroutine force
