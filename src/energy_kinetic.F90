subroutine energy_kinetic(Ek)
      
   use const
   use var_top, only : nmp, mass
   use var_state, only : velos
  
   implicit none
  
   real(PREC), intent(out) :: Ek

   integer :: imp

   Ek = 0.0e0_PREC
  
   !$omp parallel do reduction(+:Ek)
   do imp = 1, nmp
      Ek = Ek + mass(imp) * dot_product(velos(:,imp), velos(:,imp))
   enddo
   !$omp end parallel do

   Ek = 0.5e0_PREC * Ek
  
end subroutine energy_kinetic
