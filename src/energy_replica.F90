subroutine energy_replica(energies, replica_energies, flg_replica)

   use const, only : PREC
   use const_idx, only : REPT, ENE
   use var_replica, only : nrep_all, nrep_proc, irep2grep, &
                           lab2val, rep2lab, flg_repvar, get_pair, rep2val
   use var_potential, only : flg_ele
   use var_state, only : flg_bp_energy, tempK, ionic_strength

   implicit none

   real(PREC), intent(out) :: energies(0:ENE%MAX, nrep_proc)
   real(PREC), intent(out) :: replica_energies(2, nrep_all)
   logical,    intent(in)  :: flg_replica

   integer :: irep, grep
   integer :: label_own, label_opp
   real(PREC) :: ionic_strength_rep, tempK_rep
   real(PREC) :: new_energies(0:ENE%MAX)

   interface
   subroutine energy_sumup(irep, tempK_in, energies)
      use const, only : PREC
      use const_idx, only : ENE
      implicit none
      integer,    intent(in)  :: irep
      real(PREC), intent(in)  :: tempK_in
      real(PREC), intent(out) :: energies(0:ENE%MAX)
   endsubroutine energy_sumup
   subroutine set_ele(irep, tempk, ionic_strength, out_lb, out_Zp)
      use const, only : PREC
      integer, intent(in) :: irep
      real(PREC), intent(in) :: tempk
      real(PREC), intent(in) :: ionic_strength
      real(PREC), intent(out), optional :: out_lb
      real(PREC), intent(out), optional :: out_Zp
   endsubroutine set_ele
   endinterface


   energies(:,:)         = 0.0e0_PREC

   do irep = 1, nrep_proc
      if (flg_repvar(REPT%TEMP)) then
         tempK_rep = rep2val(irep2grep(irep), REPT%TEMP)
      else
         tempK_rep = tempK
      endif

      call energy_sumup(irep, tempK_rep, energies(:,irep))
   enddo


   if (.not. flg_replica) then
      return                 ! <<====== If not replica simulation, return here!
   endif

   !#############################################################################

   flg_bp_energy = .False.

   do irep = 1, nrep_proc

      grep = irep2grep(irep)
      replica_energies(1, grep) = energies(ENE%TOTAL, irep)

      label_own = rep2lab(grep)
      label_opp = get_pair(label_own)

      if (label_opp == 0 .OR. label_opp > nrep_all) then
         ! In this case, do not attempt exchange.
         replica_energies(2, grep) = 0.0e0_PREC

      else

         !#########################################
         ! Set replica variables of counterpart
         !#########################################
         if (flg_repvar(REPT%TEMP)) then
            tempK_rep = lab2val(label_opp, REPT%TEMP)
         else
            tempK_rep = tempK
         endif

         !if (flg_repvar(REPT%ION)) then
         !   ionic_strength_rep = lab2val(label_opp, REPT%ION)
         !else
            ionic_strength_rep = ionic_strength
         !endif

         !if (flg_repvar(REPT%PULL)) then
         !   pullforce = lab2val(label_opp, REPT%PULL)
         !   pull_unravel_xyz(:,1:n_pull,grep) = pullforce * pull_direction(:, 1:n_pull)
         !endif

         !#############################################
         ! Re-calculate replica-dependent coefficients
         !#############################################
         if (flg_ele) then
            call set_ele(irep, tempK_rep, ionic_strength_rep)
         endif

         !#############################################
         ! Calculate energy
         !#############################################
         call energy_sumup(irep, tempK_rep, new_energies(:))
         replica_energies(2, grep) = new_energies(ENE%TOTAL)

         !#########################################
         ! Set back replica variables
         !#########################################
         if (flg_repvar(REPT%TEMP)) then
            tempK_rep = lab2val(label_own, REPT%TEMP)
         else
            tempK_rep = tempK
         endif

         !if (flg_repvar(REPT%ION)) then
         !   ionic_strength_rep = lab2val(label_own, REPT%ION)
         !else
            ionic_strength_rep = ionic_strength
         !endif

         !if (flg_repvar(REPT%PULL)) then
         !   pullforce = lab2val(label_own, REPT%PULL)
         !   pull_unravel_xyz(:,1:n_pull,grep) = pullforce * pull_direction(:, 1:n_pull)
         !endif

         !#############################################
         ! Set back original coefficients
         !#############################################
         if (flg_ele) then
            call set_ele(irep, tempK_rep, ionic_strength_rep)
         endif
 
      endif

   enddo  ! irep

end subroutine energy_replica
