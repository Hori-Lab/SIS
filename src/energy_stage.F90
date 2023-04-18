subroutine energy_stage(irep, Estage)

   use const, only : PREC
   use const_phys, only: SMALL_VALUE
   use var_state, only : xyz
   use var_top, only : nmp
   use var_potential, only : stage_eps, stage_sigma

   implicit none
  
   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Estage

   integer :: imp
   real(PREC) :: z
   real(PREC) :: e_stage
   real(PREC) :: stage_eps4

   e_stage = 0.0e0_PREC
   stage_eps4 = 4 * stage_eps

   if (any(xyz(3,:,irep) < 0.85*stage_sigma)) then
      print '(a,f5.1,a)', 'Error: Z-coordinate < ', 0.85*stage_sigma, ' (= 0.85 * stage_sigma) in energy_stage. Abort!'
      flush(6)
      error stop
   endif

   !$omp parallel do private(imp,z) reduction(+:e_stage)
   do imp = 1, nmp

      z = xyz(3, imp, irep)

      z = stage_sigma / z
      e_stage = e_stage + stage_eps4 * (z**12 - z**6)

   enddo
   !$omp end parallel do

   Estage = e_stage

end subroutine energy_stage
