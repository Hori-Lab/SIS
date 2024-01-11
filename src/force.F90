subroutine force(irep, forces)

!$ use omp_lib
   use const
   use const_idx, only : ENE
   use var_parallel, only : nthreads
   use var_potential, only : flg_angl_ReB, flg_ele, flg_dih_cos, flg_dih_exp, bp_model, &
                             flg_stage, flg_twz, flg_bias_rg
   use var_state, only: flg_bp_MC
   use var_top, only : nmp

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(out) :: forces(3, nmp)

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

   call force_bond(irep, forces_t(1,1,tn))

   if (flg_angl_ReB) then
      call force_angl_ReB(irep, forces_t(1, 1, tn))
   else
      call force_angl(irep, forces_t(1,1,tn))
   endif

   if (flg_dih_cos) call force_dih_cos(irep, forces_t(1, 1, tn))

   if (flg_dih_exp) call force_dih_exp(irep, forces_t(1, 1, tn))
   
   if (bp_model == 4 .or. bp_model == 5) then
      if (flg_bp_MC) then
         call force_bp_limit_triplet(irep, forces_t(1,1,tn))
      else
         call force_bp_triplet(irep, forces_t(1,1,tn))
      endif
   else
      if (flg_bp_MC) then
         call force_bp_limit(irep, forces_t(1,1,tn))
      else
         call force_bp(irep, forces_t(1,1,tn))
      endif
   endif

   call force_wca(irep, forces_t(1,1,tn))

   if (flg_ele) call force_ele_DH(irep, forces_t(1,1,tn))

   if (flg_stage) call force_stage(irep, forces_t(1,1,tn))

   if (flg_bias_rg) call force_rg(irep, forces_t(1,1,tn))

!$omp end parallel

   do tn = 1, nthreads-1
      forces_t(:,:,0) = forces_t(:,:,0) + forces_t(:,:,tn)
   enddo
   forces(:,:) = forces_t(:,:,0)

   if (flg_twz) call force_tweezers(irep, forces)

end subroutine force
