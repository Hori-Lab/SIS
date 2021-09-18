subroutine force_bp()

   use const
   use const_phys, only : ZERO_JUDGE
   use var_state, only : xyz, kT, forces
   use var_potential
   use var_io, only : flg_out_bp, flg_out_bpe, KIND_OUT_BP

   implicit none
  
   integer :: ibp, imp, jmp
   integer :: imp1, imp2, imp3, imp4, imp5, imp6
   real(PREC) :: u, pre
   real(PREC) :: d, cosine, dih
   real(PREC) :: f(3,6)
   real(PREC) :: f_i(3), f_j(3), f_k(3), f_l(3)
   real(PREC) :: v12(3), v13(3), v42(3), v15(3), v62(3)
   real(PREC) :: a12
   real(PREC) :: d1212, d1313, d4242, d1213, d1242, d1215, d1515, d6262, d1262
   real(PREC) :: d1213over1212, d1242over1212, d1215over1212, d1262over1212
   real(PREC) :: m(3), n(3)
   real(PREC) :: f_bp(3, 6, nbp)

   !#######################################
   ! imp-1 (3) --- imp (1) --- imp+1 (5)
   !                ||
   ! jmp+1 (6) --- jmp (2) --- jmp-1 (4)
   !#######################################

   f_bp(:, :, :) = 0.0e0_PREC

   !$omp parallel do private(imp1, imp2, imp3, imp4, imp5, imp6, d, u, pre, cosine, dih, &
   !$omp&                    f_i, f_j, f_k, f_l, v12, v13, v42, v15, v62, a12, m, n, &
   !$omp&                    d1212, d1313, d4242, d1213, d1242, d1215, d1515, d6262, d1262, &
   !$omp&                    d1213over1212, d1242over1212, d1215over1212, d1262over1212)
   do ibp = 1, nbp

      imp1 = bp_mp(1, ibp)
      imp2 = bp_mp(2, ibp)
      
      v12(:) = pbc_vec(xyz(:,imp1) - xyz(:,imp2))
      d1212 = dot_product(v12,v12)
      a12 = sqrt(d1212)

      if (a12 >= bp_cutoff) cycle

      imp3 = imp1 - 1
      imp4 = imp2 - 1
      imp5 = imp1 + 1
      imp6 = imp2 + 1

      !===== Distance =====
      d = a12 - bp_bond_r

      u = bp_bond_k * d**2
      f_i(:) = (2.0e0_PREC * bp_bond_k * d / a12) * v12(:)
      f_bp(:, 1, ibp) = + f_i(:)
      f_bp(:, 2, ibp) = - f_i(:)

      v13(:) = pbc_vec(xyz(:, imp1) - xyz(:, imp3))
      v15(:) = pbc_vec(xyz(:, imp1) - xyz(:, imp5))
      v42(:) = pbc_vec(xyz(:, imp4) - xyz(:, imp2))
      v62(:) = pbc_vec(xyz(:, imp6) - xyz(:, imp2))
 
      d1313 = dot_product(v13, v13)
      d4242 = dot_product(v42, v42)
      d1213 = dot_product(v12, v13)
      d1242 = dot_product(v12, v42)
      d1215 = dot_product(v12, v15)
      d1515 = dot_product(v15, v15)
      d1262 = dot_product(v12, v62)
      d6262 = dot_product(v62, v62)
      d1213over1212 = d1213 / d1212
      d1242over1212 = d1242 / d1212
      d1215over1212 = d1215 / d1212
      d1262over1212 = d1262 / d1212

      !===== Angle of 3-1=2 (imp-1 -- imp -- jmp) =====
      cosine = d1213 / (sqrt(d1313) * a12)
      d = acos(cosine) - bp_angl_theta1
      pre = 2.0e0_PREC * bp_angl_k * d / sqrt(d1313*d1212 - d1213**2)
      u = u + bp_angl_k * d**2

      f_i(:) = pre * (v12(:) - (d1213 / d1313 * v13(:)))
      f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))

      f_bp(:, 1, ibp) = f_bp(:, 1, ibp) - f_i(:) - f_k(:)
      f_bp(:, 2, ibp) = f_bp(:, 2, ibp) + f_k(:)
      f_bp(:, 3, ibp) = f_i(:)

      !===== Angle of 1=2-4 (imp -- jmp -- jmp-1) =====
      cosine = d1242 / (a12 * sqrt(d4242))
      d = acos(cosine) - bp_angl_theta1
      u = u + bp_angl_k * d**2
      pre = 2.0e0_PREC * bp_angl_k * d / sqrt(d1212*d4242 - d1242**2)

      f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1242 / d4242 * v42(:)))

      f_bp(:, 1, ibp) = f_bp(:, 1, ibp) + f_i(:)
      f_bp(:, 2, ibp) = f_bp(:, 2, ibp) - f_i(:) - f_k(:)
      f_bp(:, 4, ibp) = f_k(:)

      !===== Angle of 5-1=2 (imp+1 -- imp -- jmp) =====
      cosine = d1215 / (sqrt(d1515) * a12)
      d = acos(cosine) - bp_angl_theta2
      u = u + bp_angl_k * d**2
      pre = 2.0e0_PREC * bp_angl_k * d / sqrt(d1515*d1212 - d1215**2)

      f_i(:) = pre * (v12(:) - (d1215 / d1515 * v15(:)))
      f_k(:) = pre * (v15(:) - (d1215over1212 * v12(:)))

      f_bp(:, 1, ibp) = f_bp(:, 1, ibp) - f_i(:) - f_k(:)
      f_bp(:, 2, ibp) = f_bp(:, 2, ibp) + f_k(:)
      f_bp(:, 5, ibp) = f_i(:)

      !===== Angle of 1=2-6 (imp -- jmp -- jmp+1) =====
      cosine = d1262 / (a12 * sqrt(d6262))
      d = acos(cosine) - bp_angl_theta2
      u = u + bp_angl_k * d**2
      pre = 2.0e0_PREC * bp_angl_k * d / sqrt(d1212*d6262 - d1262**2)

      f_i(:) = - pre * (v62(:) - (d1262over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1262 / d6262 * v62(:)))

      f_bp(:, 1, ibp) = f_bp(:, 1, ibp) + f_i(:)
      f_bp(:, 2, ibp) = f_bp(:, 2, ibp) - f_i(:) - f_k(:)
      f_bp(:, 6, ibp) = f_k(:)
 
      !===== Dihedral angle among 4-2=1=3 (jmp-1 -- jmp -- imp -- imp-1) =====
      m(1) = v42(2)*v12(3) - v42(3)*v12(2)
      m(2) = v42(3)*v12(1) - v42(1)*v12(3)
      m(3) = v42(1)*v12(2) - v42(2)*v12(1)
      n(1) = v12(2)*v13(3) - v12(3)*v13(2)
      n(2) = v12(3)*v13(1) - v12(1)*v13(3)
      n(3) = v12(1)*v13(2) - v12(2)*v13(1)

      dih = atan2(dot_product(v42,n)*a12, dot_product(m,n))
      d = dih + bp_dihd_phi1
      u = u + bp_dihd_k * (1.0 + cos(d))

      pre = -bp_dihd_k * sin(d) * a12
      f_i(:) = + pre / dot_product(m, m) * m(:)
      f_l(:) = - pre / dot_product(n, n) * n(:)

      f_bp(:, 4, ibp) = f_bp(:, 4, ibp) + f_i(:)
      f_bp(:, 2, ibp) = f_bp(:, 2, ibp) + (-1.0e0_PREC + d1242over1212) * f_i(:) &
                                        - (              d1213over1212) * f_l(:)
      f_bp(:, 1, ibp) = f_bp(:, 1, ibp) + (-1.0e0_PREC + d1213over1212) * f_l(:) &
                                        - (              d1242over1212) * f_i(:)
      f_bp(:, 3, ibp) = f_bp(:, 3, ibp) + f_l(:)
 
      !===== Dihedral angle among 6-2=1=5 (jmp+1 -- jmp -- imp -- imp+1) =====
      m(1) = v62(2)*v12(3) - v62(3)*v12(2)
      m(2) = v62(3)*v12(1) - v62(1)*v12(3)
      m(3) = v62(1)*v12(2) - v62(2)*v12(1)
      n(1) = v12(2)*v15(3) - v12(3)*v15(2)
      n(2) = v12(3)*v15(1) - v12(1)*v15(3)
      n(3) = v12(1)*v15(2) - v12(2)*v15(1)

      dih = atan2(dot_product(v62,n)*a12, dot_product(m,n))
      d = dih + bp_dihd_phi2
      u = u + bp_dihd_k * (1.0 + cos(d))

      pre = -bp_dihd_k * sin(d) * a12
      f_i(:) = + pre / dot_product(m, m) * m(:)
      f_l(:) = - pre / dot_product(n, n) * n(:)

      f_bp(:, 6, ibp) = f_bp(:, 6, ibp) + f_i(:)
      f_bp(:, 2, ibp) = f_bp(:, 2, ibp) + (-1.0e0_PREC + d1262over1212) * f_i(:) &
                                        - (              d1215over1212) * f_l(:)
      f_bp(:, 1, ibp) = f_bp(:, 1, ibp) + (-1.0e0_PREC + d1215over1212) * f_l(:) &
                                        - (              d1262over1212) * f_i(:)
      f_bp(:, 5, ibp) = f_bp(:, 5, ibp) + f_l(:)

      !===== Total =====
      f_bp(:, :, ibp) = bp_mp(3, ibp) * bp_U0 * exp(-u) * f_bp(:, :, ibp)

   enddo
   !$omp end parallel do

   do ibp = 1, nbp
      imp1 = bp_mp(1, ibp)
      imp2 = bp_mp(2, ibp)
      imp3 = imp1 - 1
      imp4 = imp2 - 1
      imp5 = imp1 + 1
      imp6 = imp2 + 1

      forces(:, imp1) = forces(:, imp1) + f_bp(:, 1, ibp)
      forces(:, imp2) = forces(:, imp2) + f_bp(:, 2, ibp)
      forces(:, imp3) = forces(:, imp3) + f_bp(:, 3, ibp)
      forces(:, imp4) = forces(:, imp4) + f_bp(:, 4, ibp)
      forces(:, imp5) = forces(:, imp5) + f_bp(:, 5, ibp)
      forces(:, imp6) = forces(:, imp6) + f_bp(:, 6, ibp)
   enddo

contains

   function pbc_vec(v) result (new_vec)
      
      use var_top, only : flg_pbc, pbc_box, pbc_box_half

      real(PREC) :: new_vec(3)
      real(PREC), intent(in) :: v(3)

      integer :: i

      if (.not. flg_pbc) then
         new_vec(:) = v(:)
         return
      endif

      do i = 1, 3
         if(v(i) > pbc_box_half(i)) then
            new_vec(i) = v(i) - pbc_box(i)
   
         else if(v(i) < -pbc_box_half(i)) then
            new_vec(i) = v(i) + pbc_box(i)
   
         else
            new_vec(i) = v(i)
   
         end if
      end do

   end function pbc_vec

end subroutine force_bp
