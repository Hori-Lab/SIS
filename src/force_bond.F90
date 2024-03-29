subroutine force_bond(irep, forces)
      
   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_top, only : nmp
   use var_potential, only : nbond, bond_mp, bond_k, bond_r0 !bond_para
#ifdef DUMPFORCE
   use const_idx, only : ENE
   use var_io, only : hdl_force
   use var_state, only: flg_step_dump_force
#endif
  
   implicit none
  
   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: ibd, imp1, imp2
   real(PREC) :: v(3), d, delta, f(3)
#ifdef DUMPFORCE
   integer :: i
   real(PREC) :: force_save(3, 1:nmp)

   if (flg_step_dump_force) then
      !$omp master
      force_save(:,:) = 0.0_PREC
      !$omp end master
      !$omp barrier
   endif
#endif

!$omp do private(imp1,imp2,v,d,delta,f)
   do ibd = 1, nbond
      imp1 = bond_mp(1, ibd)
      imp2 = bond_mp(2, ibd)

      v = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      d = norm2(v)
      delta = d - bond_r0
      f(:) = bond_k * delta / d * v(:)

      forces(:, imp1) = forces(:, imp1) - f(:)
      forces(:, imp2) = forces(:, imp2) + f(:)

#ifdef DUMPFORCE
      if (flg_step_dump_force) then
         do i = 1, 3
            !$omp atomic update
            force_save(i, imp1) = force_save(i, imp1) - f(i)
            !$omp atomic update
            force_save(i, imp2) = force_save(i, imp2) + f(i)
         enddo
      endif
#endif
   enddo
!$omp end do nowait


#ifdef DUMPFORCE
   if (flg_step_dump_force) then
      !$omp barrier
      !$omp master
      do imp1 = 1, nmp
         write(hdl_force(ENE%BOND), '(3(1x,e10.4))', advance='no') force_save(1:3, imp1)
      enddo
      write(hdl_force(ENE%BOND), '(a)') ''
      !$omp end master
   endif
#endif

end subroutine force_bond
