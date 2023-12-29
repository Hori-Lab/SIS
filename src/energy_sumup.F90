subroutine energy_sumup(irep, tempK_in, energies)

   use const
   use const_idx, only : ENE
   use var_state, only : flg_bp_MC
   use var_potential, only : flg_angl_ReB, flg_ele, flg_dih_cos, flg_dih_exp, flg_stage, bp_model, &
                             flg_twz

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(in) :: tempK_in
   real(PREC), intent(out) :: energies(0:ENE%MAX)

   energies(:) = 0.0e0_PREC

   call energy_bond(irep, energies(ENE%BOND))

   if (flg_angl_ReB) then
      call energy_angl_ReB(irep, energies(ENE%ANGL))
   else
      call energy_angl(irep, energies(ENE%ANGL))
   endif

   if (flg_dih_cos) call energy_dih_cos(irep, energies(ENE%DIHE))

   if (flg_dih_exp) call energy_dih_exp(irep, energies(ENE%DIHE))

   if (bp_model == 4 .or. bp_model == 5) then
      if (flg_bp_MC) then
         call energy_bp_limit_triplet(irep, tempK_in, energies(ENE%BP))
      else
         call energy_bp_triplet(irep, tempK_in, energies(ENE%BP))
      endif
   else
      if (flg_bp_MC) then
         call energy_bp_limit(irep, energies(ENE%BP))
      else
         call energy_bp(irep, energies(ENE%BP))
      endif
   endif

   call energy_wca(irep, energies(ENE%EXV))

   if (flg_ele) call energy_ele_DH(irep, energies(ENE%ELE))

   if (flg_stage) call energy_stage(irep, energies(ENE%STAGE))
   if (flg_twz) call energy_tweezers(irep, energies(ENE%TWZ))

   energies(ENE%TOTAL) = sum(energies(1:ENE%MAX))

end subroutine energy_sumup
