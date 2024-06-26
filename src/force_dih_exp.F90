subroutine force_dih_exp(irep, forces)

   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_top, only : nmp
   use var_potential, only : ndih, dih_mp, dih_k, dih_w, dih_p0
#ifdef DUMPFORCE
   use const_idx, only : ENE
   use var_io, only : hdl_force
   use var_state, only: flg_step_dump_force
#endif

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: idih, imp1, imp2, imp3, imp4
   real(PREC) :: dih, akj, akj2, pre, dvijvkj_akj2, dvklvkj_akj2
   real(PREC) :: coef1, coef2
   real(PREC) :: m(3), n(3)
   real(PREC) :: vij(3), vkj(3), vkl(3), fi(3), fj(3), fk(3), fl(3)
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

   coef1 = dih_k * dih_w
   coef2 = - 0.5_PREC * dih_w

   !$omp   do private(imp1, imp2, imp3, imp4, vij, vkj, vkl, m, n, &
   !$omp&             akj2, akj, dih, pre, fi, fj, fk, fl, dvijvkj_akj2, dvklvkj_akj2)
   do idih = 1, ndih
      imp1 = dih_mp(1, idih)
      imp2 = dih_mp(2, idih)
      imp3 = dih_mp(3, idih)
      imp4 = dih_mp(4, idih)

      vij(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      vkj(:) = pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp2, irep))
      vkl(:) = pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp4, irep))

      ! m = rij x rkj
      m(1) = vij(2)*vkj(3) - vij(3)*vkj(2)
      m(2) = vij(3)*vkj(1) - vij(1)*vkj(3)
      m(3) = vij(1)*vkj(2) - vij(2)*vkj(1)

      ! n = rkj x rkl
      n(1) = vkj(2)*vkl(3) - vkj(3)*vkl(2)
      n(2) = vkj(3)*vkl(1) - vkj(1)*vkl(3)
      n(3) = vkj(1)*vkl(2) - vkj(2)*vkl(1)

      akj2 = dot_product(vkj, vkj)
      akj = sqrt(akj2)

      dih = atan2(akj * dot_product(vij, n), dot_product(m, n))

      pre = coef1 * (dih - dih_p0) * exp( coef2 * (dih - dih_p0)**2) * akj
      fi(:) = -pre / dot_product(m,m) * m(:)
      fl(:) = +pre / dot_product(n,n) * n(:)

      dvijvkj_akj2 = dot_product(vij, vkj) / akj2
      dvklvkj_akj2 = dot_product(vkl, vkj) / akj2
      fj(:) = (-1.0_PREC + dvijvkj_akj2) * fi(:) + (          - dvklvkj_akj2) * fl(:)
      fk(:) = (          - dvijvkj_akj2) * fi(:) + (-1.0_PREC + dvklvkj_akj2) * fl(:)

      forces(:, imp1) = forces(:, imp1) + fi(:)
      forces(:, imp2) = forces(:, imp2) + fj(:)
      forces(:, imp3) = forces(:, imp3) + fk(:) 
      forces(:, imp4) = forces(:, imp4) + fl(:)

#ifdef DUMPFORCE
      if (flg_step_dump_force) then
         do i = 1, 3
            !$omp atomic update
            force_save(i, imp1) = force_save(i, imp1) + fi(i)
            !$omp atomic update
            force_save(i, imp2) = force_save(i, imp2) + fj(i)
            !$omp atomic update
            force_save(i, imp3) = force_save(i, imp3) + fk(i)
            !$omp atomic update
            force_save(i, imp4) = force_save(i, imp4) + fl(i)
         enddo
      endif
#endif
   end do 
   !$omp end do nowait

#ifdef DUMPFORCE
   if (flg_step_dump_force) then
      !$omp barrier
      !$omp master
      do imp1 = 1, nmp
         write(hdl_force(ENE%DIHE), '(3(1x,e10.4))', advance='no') force_save(1:3, imp1)
      enddo
      write(hdl_force(ENE%DIHE), '(a)') ''
      !$omp end master
   endif
#endif

end subroutine force_dih_exp
