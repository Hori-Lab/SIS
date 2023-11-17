subroutine energy_kinetic(irep, Ekinetic)
      
   use const
   use var_top, only : nmp, mass, flg_freeze, is_frozen
   use var_state, only : velos
  
   implicit none
  
   integer, intent(in) :: irep
   real(PREC), intent(out) :: Ekinetic

   integer :: imp
   real(PREC) :: e_kin

   e_kin = 0.0e0_PREC

   if (flg_freeze) then
      do imp = 1, nmp
         if (is_frozen(imp)) velos(:,imp,irep) = 0.0_PREC
      enddo
   endif
  
   !$omp parallel do reduction(+:e_kin)
   do imp = 1, nmp
      e_kin = e_kin + mass(imp) * dot_product(velos(:,imp,irep), velos(:,imp,irep))
   enddo
   !$omp end parallel do

   Ekinetic = 0.5e0_PREC * e_kin 
  
end subroutine energy_kinetic
