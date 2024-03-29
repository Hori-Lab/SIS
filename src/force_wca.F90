subroutine force_wca(irep, forces)

   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_potential, only : wca_sigma, wca_eps, wca_mp, nwca
   use var_top, only : nmp
#ifdef DUMPFORCE
   use const_idx, only : ENE
   use var_io, only : hdl_force
   use var_state, only: flg_step_dump_force
#endif

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3,nmp)
  
   integer :: iwca, imp1, imp2
   real(PREC) :: v12(3), d2, for(3)
   real(PREC) :: wca_sigma_2, coef
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

   wca_sigma_2 = wca_sigma**2
   coef = 12.0e0_PREC * wca_eps / wca_sigma_2

   !$omp do private(imp1,imp2,v12,d2,for)
   do iwca = 1, nwca(irep)

      imp1 = wca_mp(1, iwca, irep)
      imp2 = wca_mp(2, iwca, irep)

      v12(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      
      d2 = dot_product(v12, v12)

      if (d2 >= wca_sigma_2) cycle

      d2 = wca_sigma_2 / d2

      for(:) = coef * (d2**7 - d2**4) * v12(:)

      forces(:, imp1) = forces(:, imp1) + for(:)
      forces(:, imp2) = forces(:, imp2) - for(:)

#ifdef DUMPFORCE
      if (flg_step_dump_force) then
         do i = 1, 3
            !$omp atomic update
            force_save(i, imp1) = force_save(i, imp1) + for(i)
            !$omp atomic update
            force_save(i, imp2) = force_save(i, imp2) - for(i)
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
         write(hdl_force(ENE%EXV), '(3(1x,e10.4))', advance='no') force_save(1:3, imp1)
      enddo
      write(hdl_force(ENE%EXV), '(a)') ''
      !$omp end master
   endif
#endif

end subroutine force_wca
