subroutine force_bond(forces)
      
   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_top, only : nmp
   use var_potential, only : nbond, bond_mp, bond_k, bond_r0 !bond_para
  
   implicit none
  
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: ibd, imp1, imp2
   real(PREC) :: v(3), d, delta, f(3)

!$omp do private(imp1,imp2,v,d,delta,f)
   do ibd = 1, nbond
      imp1 = bond_mp(1, ibd)
      imp2 = bond_mp(2, ibd)

      v = pbc_vec_d(xyz(:, imp1), xyz(:, imp2))
      d = norm2(v)
      delta = d - bond_r0
      f(:) = bond_k * delta / d * v(:)

      forces(:, imp1) = forces(:, imp1) - f(:)
      forces(:, imp2) = forces(:, imp2) + f(:)
   enddo
!$omp end do nowait

end subroutine force_bond
