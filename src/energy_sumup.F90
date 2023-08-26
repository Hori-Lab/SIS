subroutine energy_sumup(irep, tempK_in, energies)

   use const
   use const_idx, only : ENE
   use var_potential, only : flg_angl_ReB, flg_ele, max_bp_per_nt, flg_dih_cos, flg_dih_exp, flg_stage, bp_model
   use var_state, only : tempK

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

   if (max_bp_per_nt < 1) then
      call energy_bp(irep, energies(ENE%BP))
   else
      if (bp_model == 4 .or. bp_model == 5) then
         call energy_bp_limit_triplet(irep, tempK_in, energies(ENE%BP))
      else
         call energy_bp_limit(irep, energies(ENE%BP))
      endif
   endif

   call energy_wca(irep, energies(ENE%EXV))

   if (flg_ele) call energy_ele_DH(irep, energies(ENE%ELE))

   if (flg_stage) call energy_stage(irep, energies(ENE%STAGE))

   energies(ENE%TOTAL) = sum(energies(1:ENE%MAX))

end subroutine energy_sumup
