subroutine force_bp_limit_triplet(forces)

   use mt19937_64, only : genrand64_real1, genrand64_real3
   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_top, only : nmp
   use var_state, only : xyz, bp_status, ene_bp, for_bp, kT, flg_bp_energy, nt_bp_excess
   use var_potential, only : max_bp_per_nt, nbp, bp_cutoff_energy, bp_mp, bp_paras, &
                             basepair_parameters, bp_map_dG

   implicit none

   real(PREC), intent(inout) :: forces(3, nmp)

   integer :: i, ibp, jbp, ibp_delete, nt_delete
   integer :: imp1, imp2, imp3, imp4, imp5, imp6
   integer :: i_save, i_swap
   integer :: nbp_seq
   integer :: bp_seq(nbp)
   integer :: nnt_bp_excess
   integer :: ntlist_excess(nmp)
   type(basepair_parameters) :: bpp
   real(PREC) :: u, pre, beta, ene, ratio, rnd
   real(PREC) :: d, cosine, dih
   real(PREC) :: f_i(3), f_j(3), f_k(3), f_l(3)
   real(PREC) :: v12(3), v13(3), v42(3), v15(3), v62(3)
   real(PREC) :: a12
   real(PREC) :: d1212, d1313, d4242, d1213, d1242, d1215, d1515, d6262, d1262
   real(PREC) :: d1213over1212, d1242over1212, d1215over1212, d1262over1212
   real(PREC) :: m(3), n(3)

   !real(PREC) :: for_bp(3, 6, nbp)
   ! In this subroutine, for_bp will be used in more than one loop.
   ! Therefore for_bp cannot be private and it has to be stored in var_state modele.

   !#######################################
   ! imp-1 (3) --- imp (1) --- imp+1 (5)
   !                ||
   ! jmp+1 (6) --- jmp (2) --- jmp-1 (4)
   !#######################################

   beta = 1.0_PREC / kT

   !$omp master
   bp_status(1:nbp) = .False.
   ene_bp(1:nbp) = 0.0_PREC
   for_bp(1:3,1:6,1:nbp) = 0.0_PREC
   nt_bp_excess(1:nmp) = -max_bp_per_nt
   !$omp end master

   ! Wait until the master initializes the arrays
   !$omp barrier

   !$omp do private(bpp, imp1, imp2, imp3, imp4, imp5, imp6, d, u, ene, pre, cosine, dih, &
   !$omp&           f_i, f_j, f_k, f_l, v12, v13, v42, v15, v62, a12, m, n, &
   !$omp&           d1212, d1313, d4242, d1213, d1242, d1215, d1515, d6262, d1262, &
   !$omp&           d1213over1212, d1242over1212, d1215over1212, d1262over1212)
   do ibp = 1, nbp

      imp1 = bp_mp(1, ibp)
      imp2 = bp_mp(2, ibp)
      bpp = bp_paras(bp_mp(3, ibp))

      v12(:) = pbc_vec_d(xyz(:,imp1), xyz(:,imp2))
      d1212 = dot_product(v12,v12)
      a12 = sqrt(d1212)
      d = a12 - bpp%bond_r

      if (abs(d) > bpp%cutoff_ddist) cycle

      imp3 = imp1 - 1
      imp4 = imp2 - 1
      imp5 = imp1 + 1
      imp6 = imp2 + 1

      !===== Distance =====
      !d = a12 - bpp%bond_r

      u = bpp%bond_k * d**2
      f_i(:) = (2.0e0_PREC * bpp%bond_k * d / a12) * v12(:)
      for_bp(:, 1, ibp) = + f_i(:)
      for_bp(:, 2, ibp) = - f_i(:)

      v13(:) = pbc_vec_d(xyz(:, imp1), xyz(:, imp3))
      v15(:) = pbc_vec_d(xyz(:, imp1), xyz(:, imp5))
      v42(:) = pbc_vec_d(xyz(:, imp4), xyz(:, imp2))
      v62(:) = pbc_vec_d(xyz(:, imp6), xyz(:, imp2))
 
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
      d = acos(cosine) - bpp%angl_theta1
      pre = 2.0e0_PREC * bpp%angl_k1 * d / sqrt(d1313*d1212 - d1213**2)
      u = u + bpp%angl_k1 * d**2

      f_i(:) = pre * (v12(:) - (d1213 / d1313 * v13(:)))
      f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))

      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) - f_i(:) - f_k(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) + f_k(:)
      for_bp(:, 3, ibp) = f_i(:)

      !===== Angle of 1=2-4 (imp -- jmp -- jmp-1) =====
      cosine = d1242 / (a12 * sqrt(d4242))
      d = acos(cosine) - bpp%angl_theta1
      u = u + bpp%angl_k1 * d**2
      pre = 2.0e0_PREC * bpp%angl_k1 * d / sqrt(d1212*d4242 - d1242**2)

      f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1242 / d4242 * v42(:)))

      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) + f_i(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) - f_i(:) - f_k(:)
      for_bp(:, 4, ibp) = f_k(:)

      !===== Angle of 5-1=2 (imp+1 -- imp -- jmp) =====
      cosine = d1215 / (sqrt(d1515) * a12)
      d = acos(cosine) - bpp%angl_theta2
      u = u + bpp%angl_k2 * d**2
      pre = 2.0e0_PREC * bpp%angl_k2 * d / sqrt(d1515*d1212 - d1215**2)

      f_i(:) = pre * (v12(:) - (d1215 / d1515 * v15(:)))
      f_k(:) = pre * (v15(:) - (d1215over1212 * v12(:)))

      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) - f_i(:) - f_k(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) + f_k(:)
      for_bp(:, 5, ibp) = f_i(:)

      !===== Angle of 1=2-6 (imp -- jmp -- jmp+1) =====
      cosine = d1262 / (a12 * sqrt(d6262))
      d = acos(cosine) - bpp%angl_theta2
      u = u + bpp%angl_k2 * d**2
      pre = 2.0e0_PREC * bpp%angl_k2 * d / sqrt(d1212*d6262 - d1262**2)

      f_i(:) = - pre * (v62(:) - (d1262over1212 * v12(:)))
      f_k(:) = - pre * (v12(:) - (d1262 / d6262 * v62(:)))

      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) + f_i(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) - f_i(:) - f_k(:)
      for_bp(:, 6, ibp) = f_k(:)
 
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

      for_bp(:, 4, ibp) = for_bp(:, 4, ibp) + f_i(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) + (-1.0e0_PREC + d1242over1212) * f_i(:) &
                                        - (              d1213over1212) * f_l(:)
      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) + (-1.0e0_PREC + d1213over1212) * f_l(:) &
                                        - (              d1242over1212) * f_i(:)
      for_bp(:, 3, ibp) = for_bp(:, 3, ibp) + f_l(:)
 
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

      for_bp(:, 6, ibp) = for_bp(:, 6, ibp) + f_i(:)
      for_bp(:, 2, ibp) = for_bp(:, 2, ibp) + (-1.0e0_PREC + d1262over1212) * f_i(:) &
                                        - (              d1215over1212) * f_l(:)
      for_bp(:, 1, ibp) = for_bp(:, 1, ibp) + (-1.0e0_PREC + d1215over1212) * f_l(:) &
                                        - (              d1262over1212) * f_i(:)
      for_bp(:, 5, ibp) = for_bp(:, 5, ibp) + f_l(:)

      !===== Total =====
      ene = bpp%U0 * bp_map_dG(imp1, imp2) * exp(-u)

      if (ene <= bp_cutoff_energy) then
         bp_status(ibp) = .True.
         ene_bp(ibp) = ene

         for_bp(:, :, ibp) = ene * for_bp(:, :, ibp)

         !$omp atomic
         nt_bp_excess(imp1) = nt_bp_excess(imp1) + 1
         !$omp atomic
         nt_bp_excess(imp2) = nt_bp_excess(imp2) + 1
      endif
   enddo
   !$omp end do

   !$omp master
   do     ! loop while (nnt_bp_excess > 0)

      nnt_bp_excess = 0  ! Number of nucleotides that have more than one base pair
      ntlist_excess(:) = 0
      do imp1 = 1, nmp
         if (nt_bp_excess(imp1) > 0) then
            nnt_bp_excess = nnt_bp_excess + 1
            ntlist_excess(nnt_bp_excess) = imp1
         endif
      enddo

      if (nnt_bp_excess == 0) exit

      ! Randomely choose one nucleotide (nt) that has more than one base pair
      rnd = genrand64_real3()    ! (0,1)-real-interval

      nt_delete = ntlist_excess( ceiling(rnd * nnt_bp_excess) )
      !  1 <= nt_delete <= nnt_bp_excess

      ! Generate a sequence of nucleotides that form basepairs involving nt_delete
      nbp_seq = 0
      bp_seq(:) = 0
      do ibp = 1, nbp
         if (bp_status(ibp)) then
            imp1 = bp_mp(1, ibp)
            imp2 = bp_mp(2, ibp)
            if (imp1 == nt_delete .or. imp2 == nt_delete) then
               nbp_seq = nbp_seq + 1
               bp_seq(nbp_seq) = ibp
            endif
         endif
      enddo

      ! Shuffle
      do i = 1, nbp_seq
         rnd = genrand64_real3()   ! (0,1)-real-interval
         i_swap = ceiling(rnd*nbp_seq)
         i_save = bp_seq(i)
         bp_seq(i) = bp_seq(i_swap)
         bp_seq(i_swap) = i_save
      enddo

      ! Randomely choose one "ibp" that will be deleted, depending on the energies
      ibp_delete = bp_seq(1)
      do i = 2, nbp_seq
         jbp = bp_seq(i)

         ratio = exp( (ene_bp(jbp) - ene_bp(ibp_delete)) * beta )
         rnd = genrand64_real1()  ! [0,1]-real-interval

         if (rnd < ratio) then
            ibp_delete = jbp
         endif
      enddo

      ! Delete
      bp_status(ibp_delete) = .False.
      !ene_bp(ibp_delete) = 0.0_PREC
      !! This line is commented out because ene_bp has to be kept for CHECK_FORCE.
      !! In normal run, ene_bp will never be refered when bp_status is False,
      !! thus it does not have to be zero cleared.

      ! Update nt_bp_excess
      nt_bp_excess(bp_mp(1, ibp_delete)) = nt_bp_excess(bp_mp(1, ibp_delete)) - 1
      nt_bp_excess(bp_mp(2, ibp_delete)) = nt_bp_excess(bp_mp(2, ibp_delete)) - 1
   enddo

   !$omp end master

   ! Wait until the master finishes deletions
   !$omp barrier

   !$omp do private(imp1, imp2, imp3, imp4, imp5, imp6)
   do ibp = 1, nbp

      if (.not. bp_status(ibp)) cycle

      imp1 = bp_mp(1, ibp)
      imp2 = bp_mp(2, ibp)
      imp3 = imp1 - 1
      imp4 = imp2 - 1
      imp5 = imp1 + 1
      imp6 = imp2 + 1

      forces(:, imp1) = forces(:, imp1) + for_bp(:, 1, ibp)
      forces(:, imp2) = forces(:, imp2) + for_bp(:, 2, ibp)
      forces(:, imp3) = forces(:, imp3) + for_bp(:, 3, ibp)
      forces(:, imp4) = forces(:, imp4) + for_bp(:, 4, ibp)
      forces(:, imp5) = forces(:, imp5) + for_bp(:, 5, ibp)
      forces(:, imp6) = forces(:, imp6) + for_bp(:, 6, ibp)
   end do
   !$omp end do nowait

   !$omp master
   flg_bp_energy = .True.
   !$omp end master

end subroutine force_bp_limit_triplet
