subroutine energy_bp(irep, Ebp)

   use const
   use const_phys, only : ZERO_JUDGE
   use pbc, only : pbc_vec_d
   use var_state, only : xyz, bp_status, bp_status_MC, ene_bp, nstep_bp_MC, flg_bp_energy
   use var_potential, only : bp_paras, nbp, bp_mp, basepair_parameters, bp_cutoff_energy

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(inout) :: Ebp

   integer :: ibp, imp1, imp2, imp3, imp4, imp5, imp6
   type(basepair_parameters) :: bpp
   real(PREC) :: u
   real(PREC) :: d, theta, phi

   if (.not. flg_bp_energy) then

      ene_bp(1:nbp(irep), irep) = 0.0e0_PREC
      bp_status(1:nbp(irep), irep) = .False.

      !$omp parallel do private(imp1, imp2, imp3, imp4, imp5, imp6, d, u, theta, phi, bpp)
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

         d = norm2(pbc_vec_d(xyz(:, imp1, irep), xyz(:, imp2, irep))) - bpp%bond_r

         if (abs(d) > bpp%cutoff_ddist) cycle

         u = bpp%bond_k * d**2

         theta = mp_angle(imp1, imp2, imp4)
         u = u + bpp%angl_k1 * (theta - bpp%angl_theta1)**2

         theta = mp_angle(imp3, imp1, imp2)
         u = u + bpp%angl_k2 * (theta - bpp%angl_theta2)**2

         theta = mp_angle(imp1, imp2, imp6)
         u = u + bpp%angl_k3 * (theta - bpp%angl_theta3)**2

         theta = mp_angle(imp5, imp1, imp2)
         u = u + bpp%angl_k4 * (theta - bpp%angl_theta4)**2

         phi = mp_dihedral(imp3, imp1, imp2, imp4)
         u = u + bpp%dihd_k1 * (1.0_PREC + cos(phi + bpp%dihd_phi1))

         phi = mp_dihedral(imp5, imp1, imp2, imp6)
         u = u + bpp%dihd_k2 * (1.0_PREC + cos(phi + bpp%dihd_phi2))

         u = bpp%U0 * exp(-u)
         if (u <= bp_cutoff_energy) then
            ene_bp(ibp, irep) = u
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

end subroutine energy_bp
