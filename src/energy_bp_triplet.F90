subroutine energy_bp_triplet(irep, tempK_in, Ebp)

   !use mt19937_64, only : genrand64_real1, genrand64_real3
   use mt_stream
   use const, only : PREC
   use pbc, only : pbc_vec_d
   use var_state, only : xyz, bp_status, ene_bp, flg_bp_energy, temp_independent
   use var_potential, only : bp_cutoff_energy, nbp, bp_mp, bp_paras, bp_coef, &
                             basepair_parameters

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(in) :: tempK_in
   real(PREC), intent(inout) :: Ebp

   integer :: ibp
   integer :: imp, jmp
   type(basepair_parameters) :: bpp
   real(PREC) :: tK, dG
   real(PREC) :: u
   real(PREC) :: d, theta, phi
   real(PREC) :: ene

   if (.not. flg_bp_energy) then

      if (temp_independent == 0) then
         tK = tempK_in
      else
         tK = 0.0_PREC
      endif

      bp_status(1:nbp(irep), irep) = .False.
      ene_bp(1:nbp(irep), irep) = 0.0_PREC

      !$omp barrier

      !$omp parallel do private(imp, jmp, d, u, theta, phi, ene, bpp, dG)
      do ibp = 1, nbp(irep)

         imp = bp_mp(1, ibp, irep)
         jmp = bp_mp(2, ibp, irep)
         bpp = bp_paras(bp_mp(3, ibp, irep))

         ! dG = dH - T * dS
         dG = bp_coef(1, ibp, irep) - tK * bp_coef(2, ibp, irep)
         if (dG >= 0.0_PREC) cycle

         d = norm2(pbc_vec_d(xyz(:,imp,irep), xyz(:, jmp,irep))) - bpp%bond_r

         if (abs(d) > bpp%cutoff_ddist) cycle

         u = bpp%bond_k * d**2

         theta = mp_angle(imp, jmp, jmp-1)
         u = u + bpp%angl_k1 * (theta - bpp%angl_theta1)**2

         theta = mp_angle(imp-1, imp, jmp)
         u = u + bpp%angl_k2 * (theta - bpp%angl_theta2)**2

         theta = mp_angle(imp, jmp, jmp+1)
         u = u + bpp%angl_k3 * (theta - bpp%angl_theta3)**2

         theta = mp_angle(imp+1, imp, jmp)
         u = u + bpp%angl_k4 * (theta - bpp%angl_theta4)**2

         phi = mp_dihedral(imp-1, imp, jmp, jmp-1)
         u = u + bpp%dihd_k1 * (1.0_PREC + cos(phi + bpp%dihd_phi1))

         phi = mp_dihedral(imp+1, imp, jmp, jmp+1)
         u = u + bpp%dihd_k2 * (1.0_PREC + cos(phi + bpp%dihd_phi2))

         ene = bpp%U0 * dG * exp(-u)

         if (ene <= bp_cutoff_energy) then
            ene_bp(ibp, irep) = ene
            bp_status(ibp, irep) = .True.
         endif

      enddo
      !$omp end parallel do

   endif

   Ebp = sum(ene_bp(1:nbp(irep), irep))

contains

   pure function mp_angle(imp1, imp2, imp3) result(theta)

      real(PREC) :: theta
      integer, intent(in) :: imp1, imp2, imp3
      real(PREC) :: v12(3), v32(3)
      real(PREC) :: co

      v12(:) = pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      v32(:) = pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp2, irep))

      co = dot_product(v32, v12) / sqrt(dot_product(v12,v12) * dot_product(v32,v32))

      if(co > 1.0e0_PREC) then
         co = 1.0e0_PREC
      else if(co < -1.0e0_PREC) then
         co = -1.0e0_PREC
      end if

      theta = acos(co)

   endfunction mp_angle

   pure function mp_dihedral(imp1, imp2, imp3, imp4) result(phi)

      real(PREC) :: phi
      integer, intent(in) :: imp1, imp2, imp3, imp4
      real(PREC) :: v12(3), v32(3), v34(3), m(3), n(3)

      v12(:) =  pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))
      v32(:) =  pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp2, irep))
      v34(:) =  pbc_vec_d(xyz(:, imp3, irep), xyz(:, imp4, irep))

      m(1) = v12(2)*v32(3) - v12(3)*v32(2)
      m(2) = v12(3)*v32(1) - v12(1)*v32(3)
      m(3) = v12(1)*v32(2) - v12(2)*v32(1)
      n(1) = v32(2)*v34(3) - v32(3)*v34(2)
      n(2) = v32(3)*v34(1) - v32(1)*v34(3)
      n(3) = v32(1)*v34(2) - v32(2)*v34(1)

      phi = atan2(dot_product(v12,n)*norm2(v32), dot_product(m,n))

   endfunction mp_dihedral

end subroutine energy_bp_triplet
