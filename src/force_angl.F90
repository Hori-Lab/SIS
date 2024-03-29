subroutine force_angl(irep, forces)

   use const
   use pbc, only : pbc_vec_d
   use var_state, only : xyz
   use var_top, only : nmp
   use var_potential, only : nangl, angl_mp, angl_k, angl_t0
#ifdef DUMPFORCE
   use const_idx, only : ENE
   use var_io, only : hdl_force
   use var_state, only: flg_step_dump_force
#endif

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: ibd, imp1, imp2, imp3
   real(PREC) :: t, f, delta, cosine
   real(PREC) :: d21, d32, d2132
   real(PREC) :: v21(3), v32(3), f21(3), f32(3)
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

!$omp do private(imp1,imp2,imp3,v21,v32,d21,d32,d2132,cosine,delta,t,f,f21,f32)
   do ibd = 1, nangl
      imp1 = angl_mp(1, ibd)
      imp2 = angl_mp(2, ibd)
      imp3 = angl_mp(3, ibd)

      v21(:) = pbc_vec_d(xyz(:, imp2, irep), xyz(:, imp1, irep))
      v32(:) = pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp2, irep))

      d21 = dot_product(v21, v21)
      d32 = dot_product(v32, v32)
      d2132 = dot_product(v32, v21)

      cosine = - d2132 / sqrt(d21 * d32)

      if (abs(cosine) > 1.0e0_PREC) then
         cosine = sign(1.0e0_PREC, cosine)
      endif

      delta = acos(cosine) - angl_t0

      t = d21 * d32 - d2132**2
      if (t <= 1.0e0_PREC) then
         t = 1.0e0_PREC
      end if

      f = angl_k * delta / sqrt(t)
      f21(:) = f * (v21(:) * (d2132 / d21) - v32(:))
      f32(:) = f * (v32(:) * (d2132 / d32) - v21(:))

      forces(:, imp1) = forces(:, imp1) - f21(:)
      forces(:, imp2) = forces(:, imp2) + f21(:) - f32(:)
      forces(:, imp3) = forces(:, imp3) + f32(:)

      !Eangl = Eangl + 0.5 * angl_k * (t - angl_t0) ** 2
#ifdef DUMPFORCE
      if (flg_step_dump_force) then
         do i = 1, 3
            !$omp atomic update
            force_save(i, imp1) = force_save(i, imp1) - f21(i)
            !$omp atomic update
            force_save(i, imp2) = force_save(i, imp2) + f21(i) - f32(i)
            !$omp atomic update
            force_save(i, imp3) = force_save(i, imp3) + f32(i)
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
         write(hdl_force(ENE%ANGL), '(3(1x,e10.4))', advance='no') force_save(1:3, imp1)
      enddo
      write(hdl_force(ENE%ANGL), '(a)') ''
      !$omp end master
   endif
#endif

end subroutine force_angl
