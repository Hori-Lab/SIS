subroutine force_bp_triplet(irep, forces)

   !use mt19937_64, only : genrand64_real1, genrand64_real3
   use mt_stream
   use const, only : PREC
   use const_idx, only : REPT
   use const_phys, only : BOLTZ_KCAL_MOL
   use pbc, only : pbc_vec_d
   use var_top, only : nmp
   use var_state, only : xyz, bp_status, ene_bp, kT, flg_bp_energy, tempK, &
                         nstep_bp_MC, bp_status_MC
   use var_potential, only : nbp, bp_cutoff_energy, bp_mp, bp_paras, bp_coef, &
                             basepair_parameters, flg_bias_ss, bias_ss_force
   use var_replica, only : flg_repvar, rep2val, irep2grep

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: ibp
   integer :: imp1, imp2, imp3, imp4, imp5, imp6
   type(basepair_parameters) :: bpp
   real(PREC) :: tK, dG
   real(PREC) :: u, pre, beta, ene
   real(PREC) :: d, cosine, dih
   real(PREC) :: f_i(3), f_j(3), f_k(3), f_l(3)
   real(PREC) :: v12(3), v13(3), v42(3), v15(3), v62(3)
   real(PREC) :: a12
   real(PREC) :: d1212, d1313, d4242, d1213, d1242, d1215, d1515, d6262, d1262
   real(PREC) :: d1213over1212, d1242over1212, d1215over1212, d1262over1212
   real(PREC) :: m(3), n(3)
   real(PREC) :: f_bp(3, 6)  ! This can be private in each thread. Do not have to use for_bp in var_state.

   !#######################################
   ! imp-1 (3) --- imp (1) --- imp+1 (5)
   !                ||
   ! jmp+1 (6) --- jmp (2) --- jmp-1 (4)
   !#######################################

   if (flg_repvar(REPT%TEMP)) then
      tK = rep2val(irep2grep(irep), REPT%TEMP)
      beta = 1.0_PREC / (BOLTZ_KCAL_MOL * tK)
   else
      tK = tempK
      beta = 1.0_PREC / kT
   endif

   !$omp master
   bp_status(1:nbp(irep), irep) = .False.
   ene_bp(1:nbp(irep), irep) = 0.0_PREC
   !$omp end master

   ! Wait until the master initializes the arrays
   !$omp barrier

   !$omp do private(bpp, imp1, imp2, imp3, imp4, imp5, imp6, d, u, ene, pre, cosine, dih, &
   !$omp&           f_i, f_j, f_k, f_l, v12, v13, v42, v15, v62, a12, m, n, &
   !$omp&           d1212, d1313, d4242, d1213, d1242, d1215, d1515, d6262, d1262, &
   !$omp&           d1213over1212, d1242over1212, d1215over1212, d1262over1212, &
   !$omp&           dG, f_bp)
   do ibp = 1, nbp(irep)

      if (nstep_bp_MC > 0) then
         if (.not. bp_status_MC(ibp, irep)) cycle
      endif

      imp1 = bp_mp(1, ibp, irep)
      imp2 = bp_mp(2, ibp, irep)
      imp3 = bp_mp(3, ibp, irep)
      imp4 = bp_mp(4, ibp, irep)
      imp5 = bp_mp(5, ibp, irep)
      imp6 = bp_mp(6, ibp, irep)
      bpp = bp_paras(bp_mp(7, ibp, irep))

      ! dG = dH - T * dS
      dG = bp_coef(1, ibp, irep) - tK * bp_coef(2, ibp, irep)
      if (dG >= 0.0_PREC) cycle

      v12(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      d1212 = dot_product(v12,v12)
      a12 = sqrt(d1212)
      d = a12 - bpp%bond_r

      if (flg_bias_ss) then
         if (abs(d) > bpp%cutoff_ddist) then
            v12(:) = bias_ss_force / a12 * v12(:)
            !$omp critical
            forces(:, imp1) = forces(:, imp1) - v12(:)
            forces(:, imp2) = forces(:, imp2) + v12(:)
            !$omp end critical
            cycle
         endif
      else
         if (abs(d) > bpp%cutoff_ddist) cycle
      endif

      !===== Distance =====
      !d = a12 - bpp%bond_r

      u = bpp%bond_k * d**2
      f_i(:) = (2.0e0_PREC * bpp%bond_k * d / a12) * v12(:)
      f_bp(:, 1) = + f_i(:)
      f_bp(:, 2) = - f_i(:)

      v13(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp3, irep))
      v15(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp5, irep))
      v42(:) = pbc_vec_d(xyz(:, imp4, irep), xyz(:, imp2, irep))
      v62(:) = pbc_vec_d(xyz(:, imp6, irep), xyz(:, imp2, irep))
 
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
      d = acos(cosine) - bpp%angl_theta2
      pre = 2.0e0_PREC * bpp%angl_k2 * d / sqrt(d1313*d1212 - d1213**2)
      u = u + bpp%angl_k2 * d**2

      f_i(:) = pre * (v12(:) - (d1213 / d1313 * v13(:)))
      f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))

      f_bp(:, 1) = f_bp(:, 1) - f_i(:) - f_k(:)
      f_bp(:, 2) = f_bp(:, 2) + f_k(:)
      f_bp(:, 3) = f_i(:)

      !===== Angle of 1=2-4 (imp -- jmp -- jmp-1) =====
      cosine = d1242 / (a12 * sqrt(d4242))
      d = acos(cosine) - bpp%angl_theta1
      u = u + bpp%angl_k1 * d**2
      pre = 2.0e0_PREC * bpp%angl_k1 * d / sqrt(d1212*d4242 - d1242**2)

      f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1242 / d4242 * v42(:)))

      f_bp(:, 1) = f_bp(:, 1) + f_i(:)
      f_bp(:, 2) = f_bp(:, 2) - f_i(:) - f_k(:)
      f_bp(:, 4) = f_k(:)

      !===== Angle of 5-1=2 (imp+1 -- imp -- jmp) =====
      cosine = d1215 / (sqrt(d1515) * a12)
      d = acos(cosine) - bpp%angl_theta4
      u = u + bpp%angl_k4 * d**2
      pre = 2.0e0_PREC * bpp%angl_k4 * d / sqrt(d1515*d1212 - d1215**2)

      f_i(:) = pre * (v12(:) - (d1215 / d1515 * v15(:)))
      f_k(:) = pre * (v15(:) - (d1215over1212 * v12(:)))

      f_bp(:, 1) = f_bp(:, 1) - f_i(:) - f_k(:)
      f_bp(:, 2) = f_bp(:, 2) + f_k(:)
      f_bp(:, 5) = f_i(:)

      !===== Angle of 1=2-6 (imp -- jmp -- jmp+1) =====
      cosine = d1262 / (a12 * sqrt(d6262))
      d = acos(cosine) - bpp%angl_theta3
      u = u + bpp%angl_k3 * d**2
      pre = 2.0e0_PREC * bpp%angl_k3 * d / sqrt(d1212*d6262 - d1262**2)

      f_i(:) = - pre * (v62(:) - (d1262over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1262 / d6262 * v62(:)))

      f_bp(:, 1) = f_bp(:, 1) + f_i(:)
      f_bp(:, 2) = f_bp(:, 2) - f_i(:) - f_k(:)
      f_bp(:, 6) = f_k(:)
 
      !===== Dihedral angle among 4-2=1=3 (jmp-1 -- jmp -- imp -- imp-1) =====
      m(1) = v42(2)*v12(3) - v42(3)*v12(2)
      m(2) = v42(3)*v12(1) - v42(1)*v12(3)
      m(3) = v42(1)*v12(2) - v42(2)*v12(1)
      n(1) = v12(2)*v13(3) - v12(3)*v13(2)
      n(2) = v12(3)*v13(1) - v12(1)*v13(3)
      n(3) = v12(1)*v13(2) - v12(2)*v13(1)

      dih = atan2(dot_product(v42,n)*a12, dot_product(m,n))
      d = dih + bpp%dihd_phi1
      u = u + bpp%dihd_k1 * (1.0 + cos(d))

      pre = -bpp%dihd_k1 * sin(d) * a12
      f_i(:) = + pre / dot_product(m, m) * m(:)
      f_l(:) = - pre / dot_product(n, n) * n(:)

      f_bp(:, 4) = f_bp(:, 4) + f_i(:)
      f_bp(:, 2) = f_bp(:, 2) + (-1.0e0_PREC + d1242over1212) * f_i(:) &
                              - (              d1213over1212) * f_l(:)
      f_bp(:, 1) = f_bp(:, 1) + (-1.0e0_PREC + d1213over1212) * f_l(:) &
                              - (              d1242over1212) * f_i(:)
      f_bp(:, 3) = f_bp(:, 3) + f_l(:)
 
      !===== Dihedral angle among 6-2=1=5 (jmp+1 -- jmp -- imp -- imp+1) =====
      m(1) = v62(2)*v12(3) - v62(3)*v12(2)
      m(2) = v62(3)*v12(1) - v62(1)*v12(3)
      m(3) = v62(1)*v12(2) - v62(2)*v12(1)
      n(1) = v12(2)*v15(3) - v12(3)*v15(2)
      n(2) = v12(3)*v15(1) - v12(1)*v15(3)
      n(3) = v12(1)*v15(2) - v12(2)*v15(1)

      dih = atan2(dot_product(v62,n)*a12, dot_product(m,n))
      d = dih + bpp%dihd_phi2
      u = u + bpp%dihd_k2 * (1.0 + cos(d))

      pre = -bpp%dihd_k2 * sin(d) * a12
      f_i(:) = + pre / dot_product(m, m) * m(:)
      f_l(:) = - pre / dot_product(n, n) * n(:)

      f_bp(:, 6) = f_bp(:, 6) + f_i(:)
      f_bp(:, 2) = f_bp(:, 2) + (-1.0e0_PREC + d1262over1212) * f_i(:) &
                              - (              d1215over1212) * f_l(:)
      f_bp(:, 1) = f_bp(:, 1) + (-1.0e0_PREC + d1215over1212) * f_l(:) &
                              - (              d1262over1212) * f_i(:)
      f_bp(:, 5) = f_bp(:, 5) + f_l(:)

      !===== Total =====
      !ene = bpp%U0 * bp_map_dG(imp1, imp2, irep) * exp(-u)
      ene = bpp%U0 * dG * exp(-u)

      if (ene <= bp_cutoff_energy) then
         bp_status(ibp, irep) = .True.
         ene_bp(ibp, irep) = ene

         f_bp(:, :) = ene * f_bp(:, :)

         forces(:, imp1) = forces(:, imp1) + f_bp(:, 1)
         forces(:, imp2) = forces(:, imp2) + f_bp(:, 2)
         forces(:, imp3) = forces(:, imp3) + f_bp(:, 3)
         forces(:, imp4) = forces(:, imp4) + f_bp(:, 4)
         forces(:, imp5) = forces(:, imp5) + f_bp(:, 5)
         forces(:, imp6) = forces(:, imp6) + f_bp(:, 6)
      endif

   enddo
   !$omp end do nowait

   !$omp master
   flg_bp_energy = .True.
   !$omp end master

end subroutine force_bp_triplet
