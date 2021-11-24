subroutine force_angl()

   use const
   use pbc, only : pbc_vec
   use var_state, only : xyz, forces
   use var_potential, only : nangl, angl_mp, angl_k, angl_t0

   implicit none
  
   integer :: ibd, imp1, imp2, imp3
   real(PREC) :: t, f, delta, cosine
   real(PREC) :: d21, d32, d2132
   real(PREC) :: v21(3), v32(3), f21(3), f32(3)

   do ibd = 1, nangl
      imp1 = angl_mp(1, ibd)
      imp2 = angl_mp(2, ibd)
      imp3 = angl_mp(3, ibd)

      v21(:) = pbc_vec(xyz(:, imp2) - xyz(:, imp1))
      v32(:) = pbc_vec(xyz(:, imp3) - xyz(:, imp2))
     
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
   enddo

end subroutine force_angl
