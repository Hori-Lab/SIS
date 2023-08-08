subroutine energy_bp(Ebp)

   use :: ieee_exceptions, only : IEEE_GET_HALTING_MODE, IEEE_SET_HALTING_MODE, IEEE_UNDERFLOW

   use const
   use const_phys, only : ZERO_JUDGE
   use pbc, only : pbc_vec_d
   use var_state, only : xyz, kT
   use var_potential, only : bp_paras, nbp, bp_mp, bp_type2nhb, basepair_parameters
   use var_io, only : flg_out_bp, flg_out_bpall, flg_out_bpe, hdl_bp, hdl_bpall, hdl_bpe, KIND_OUT_BP, KIND_OUT_BPE

   implicit none
  
   real(PREC), intent(inout) :: Ebp

   integer :: ibp, imp, jmp, nhb
   type(basepair_parameters) :: bpp
   real(PREC) :: u
   real(PREC) :: d, theta, phi
   real(PREC) :: e_bp(nbp)
   logical :: halt_mode

   e_bp(:) = 0.0e0_PREC

   !$omp parallel do private(imp, jmp, d, u, theta, phi, bpp)
   do ibp = 1, nbp

      imp = bp_mp(1, ibp)
      jmp = bp_mp(2, ibp)
      bpp = bp_paras(bp_mp(3, ibp))
      
      d = norm2(pbc_vec_d(xyz(:,imp), xyz(:, jmp))) - bpp%bond_r

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

      e_bp(ibp) = bpp%U0 * exp(-u)

   enddo
   !$omp end parallel do

   Ebp = sum(e_bp)

   if (flg_out_bp) then

      call ieee_get_halting_mode(IEEE_UNDERFLOW, halt_mode)
      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=.false. )

      do ibp = 1, nbp

         nhb = bp_type2nhb(bp_mp(3, ibp))

         if (e_bp(ibp) < - nhb * kT) then
            imp = bp_mp(1, ibp)
            jmp = bp_mp(2, ibp)
            write(hdl_bp) int(imp,kind=KIND_OUT_BP), int(jmp,kind=KIND_OUT_BP), real(e_bp(ibp), kind=KIND_OUT_BPE)
         endif
      enddo

      write(hdl_bp) int(0,kind=KIND_OUT_BP), int(0,kind=KIND_OUT_BP), real(0.0, kind=KIND_OUT_BPE)

      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=halt_mode)
   endif

   if (flg_out_bpall) then

      call ieee_get_halting_mode(IEEE_UNDERFLOW, halt_mode)
      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=.false. )

      do ibp = 1, nbp

         if (e_bp(ibp) < -ZERO_JUDGE) then  ! To output all
            imp = bp_mp(1, ibp)
            jmp = bp_mp(2, ibp)
            write(hdl_bpall) int(imp,kind=KIND_OUT_BP), int(jmp,kind=KIND_OUT_BP), real(e_bp(ibp), kind=KIND_OUT_BPE)
         endif
      enddo

      write(hdl_bpall) int(0,kind=KIND_OUT_BP), int(0,kind=KIND_OUT_BP), real(0.0, kind=KIND_OUT_BPE)

      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=halt_mode)
   endif


   if (flg_out_bpe) then

      do ibp = 1, nbp

         !nhb = bp_type2nhb(bp_mp(3, ibp))

         if (e_bp(ibp) < -ZERO_JUDGE) then  ! To output all
         !if (e_bp(ibp) < - nhb * kT) then
            imp = bp_mp(1, ibp)
            jmp = bp_mp(2, ibp)
            write(hdl_bpe, '(1x,i5,1x,i5,1x,f5.2)', advance='no') imp, jmp, e_bp(ibp)
         endif
      enddo

      write(hdl_bpe, '(a)') ''
   endif

contains

   pure function mp_angle(imp1, imp2, imp3) result(theta)

      real(PREC) :: theta
      integer, intent(in) :: imp1, imp2, imp3
      real(PREC) :: v12(3), v32(3)
      real(PREC) :: co

      v12(:) = pbc_vec_d(xyz(:, imp1), xyz(:, imp2))
      v32(:) = pbc_vec_d(xyz(:, imp3), xyz(:, imp2))

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

      v12(:) =  pbc_vec_d(xyz(:, imp1), xyz(:, imp2))
      v32(:) =  pbc_vec_d(xyz(:, imp3), xyz(:, imp2))
      v34(:) =  pbc_vec_d(xyz(:, imp3), xyz(:, imp4))

      m(1) = v12(2)*v32(3) - v12(3)*v32(2)
      m(2) = v12(3)*v32(1) - v12(1)*v32(3)
      m(3) = v12(1)*v32(2) - v12(2)*v32(1)
      n(1) = v32(2)*v34(3) - v32(3)*v34(2)
      n(2) = v32(3)*v34(1) - v32(1)*v34(3)
      n(3) = v32(1)*v34(2) - v32(2)*v34(1)

      phi = atan2(dot_product(v12,n)*norm2(v32), dot_product(m,n))

   endfunction mp_dihedral

end subroutine energy_bp
