subroutine energy_kinetic()
      
   use const
   use var_top, only : nmp, mass
   use var_state, only : velos, Ekinetic
  
   implicit none
  
   integer :: imp

   Ekinetic = 0.0e0_PREC
  
   !$omp parallel do reduction(+:Ekinetic)
   do imp = 1, nmp
      Ekinetic = Ekinetic + mass(imp) * dot_product(velos(:,imp), velos(:,imp))
   enddo
   !$omp end parallel do

   Ekinetic = 0.5e0_PREC * Ekinetic
  
end subroutine energy_kinetic
