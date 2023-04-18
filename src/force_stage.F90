subroutine force_stage(irep, forces)

   use const, only : PREC
   use const_phys, only : SMALL_VALUE
   use var_state, only : xyz
   use var_potential, only : stage_sigma, stage_eps
   use var_top, only : nmp

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3,nmp)
  
   integer :: imp
   real(PREC) :: z, coef, f_z
   
   coef = 24.0e0_PREC * stage_eps / stage_sigma

   !if (any(xyz(3,:) < 0.85*stage_sigma)) then
   !   print '(a,f5.1,a)', 'Error: Z-coordinate < ', 0.85*stage_sigma, ' (= 0.85 * stage_sigma) in force_stage. Abort!'
   !   flush(6)
   !   error stop
   !endif

   !$omp do private(imp,z,f_z)
   do imp = 1, nmp

      z = xyz(3, imp, irep)

      if (z < SMALL_VALUE) z = SMALL_VALUE

      z = stage_sigma / z

      f_z = coef * (2*z**13 - z**7)

      forces(3, imp) = forces(3, imp) +  f_z
   enddo
   !$omp end do nowait

end subroutine force_stage
