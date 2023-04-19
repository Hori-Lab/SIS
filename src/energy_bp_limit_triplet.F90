subroutine energy_bp_limit_triplet(irep, tempK_in, Ebp)

   use :: ieee_exceptions, only : IEEE_GET_HALTING_MODE, IEEE_SET_HALTING_MODE, IEEE_UNDERFLOW

   !use mt19937_64, only : genrand64_real1, genrand64_real3
   use mt_stream
   use const, only : PREC
   use const_idx, only : REPT
   use const_phys, only : BOLTZ_KCAL_MOL
   use pbc, only : pbc_vec_d
   use var_top, only : nmp
   use var_state, only : xyz, bp_status, ene_bp, flg_bp_energy, mts, temp_independent
   use var_potential, only : max_bp_per_nt, bp_cutoff_energy, nbp, bp_mp, bp_paras, bp_coef, &
                             basepair_parameters
   use var_io, only : flg_out_bp, flg_out_bpall, flg_out_bpe, hdl_bp, hdl_bpall, hdl_bpe, KIND_OUT_BP, KIND_OUT_BPE
   use var_replica, only : flg_repvar, rep2val, irep2grep

   implicit none
  
   integer, intent(in) :: irep
   real(PREC), intent(in) :: tempK_in
   real(PREC), intent(inout) :: Ebp

   integer :: i, ibp, jbp
   integer :: nt_delete, ibp_delete
   integer :: imp, jmp
   integer :: i_save, i_swap
   type(basepair_parameters) :: bpp
   real(PREC) :: tK, dG
   real(PREC) :: u, beta, ratio
   real(PREC) :: d, theta, phi
   real(PREC) :: ene
   real(PREC) :: rnd
   integer :: nt_bp_excess(nmp)
   integer :: nbp_seq
   integer :: bp_seq(nbp(irep))
   integer :: nnt_bp_excess
   integer :: ntlist_excess(nmp)
   logical :: halt_mode

   if (flg_repvar(REPT%TEMP)) then
      tK = rep2val(irep2grep(irep), REPT%TEMP)
      beta = 1.0_PREC / (BOLTZ_KCAL_MOL * tK)
   else
      tK = tempK_in
      beta = 1.0_PREC / (BOLTZ_KCAL_MOL * tK)
   endif

   if (temp_independent == 0) then
      tK = tempK_in
      beta = 1.0_PREC / (BOLTZ_KCAL_MOL * tK)
   else
      tK = 0.0_PREC
      beta = HUGE(beta)   ! so that the highest energy BP will be deleted.
   endif

   if (.not. flg_bp_energy) then

      bp_status(1:nbp(irep), irep) = .False.
      nt_bp_excess(1:nmp) = -max_bp_per_nt
      ene_bp(1:nbp(irep)) = 0.0_PREC

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
         u = u + bpp%angl_k1 * (theta - bpp%angl_theta1)**2

         theta = mp_angle(imp, jmp, jmp+1)
         u = u + bpp%angl_k2 * (theta - bpp%angl_theta2)**2

         theta = mp_angle(imp+1, imp, jmp)
         u = u + bpp%angl_k2 * (theta - bpp%angl_theta2)**2

         phi = mp_dihedral(imp-1, imp, jmp, jmp-1)
         u = u + bpp%dihd_k1 * (1.0_PREC + cos(phi + bpp%dihd_phi1))

         phi = mp_dihedral(imp+1, imp, jmp, jmp+1)
         u = u + bpp%dihd_k2 * (1.0_PREC + cos(phi + bpp%dihd_phi2))

         ene = bpp%U0 * dG * exp(-u)

         if (ene <= bp_cutoff_energy) then
            ene_bp(ibp) = ene
            bp_status(ibp, irep) = .True.

            !$omp atomic
            nt_bp_excess(imp) = nt_bp_excess(imp) + 1
            !$omp atomic
            nt_bp_excess(jmp) = nt_bp_excess(jmp) + 1
         endif

      enddo
      !$omp end parallel do


      do   ! loop while (nnt_bp_excess > 0)

         nnt_bp_excess = 0  ! Number of nucleotides that have more than one base pair
         ntlist_excess(:) = 0
         do imp = 1, nmp
            if (nt_bp_excess(imp) > 0) then
               nnt_bp_excess = nnt_bp_excess + 1
               ntlist_excess(nnt_bp_excess) = imp
            endif
         enddo

         if (nnt_bp_excess == 0) exit

         ! Randomely choose one nucleotide (nt) that has more than one base pair
         !rnd = genrand64_real3()    ! (0,1)-real-interval
         rnd = genrand_double3(mts(irep))    ! (0,1)-real-interval

         nt_delete = ntlist_excess( ceiling(rnd * nnt_bp_excess) )
         !  1 <= nt_delete <= nnt_bp_excess

         ! Generate a sequence of nucleotides that form basepairs involving nt_delete
         nbp_seq = 0
         bp_seq(:) = 0
         do ibp = 1, nbp(irep)
            if (bp_status(ibp, irep)) then
               imp = bp_mp(1, ibp, irep)
               jmp = bp_mp(2, ibp, irep)
               if (imp == nt_delete .or. jmp == nt_delete) then
                  nbp_seq = nbp_seq + 1
                  bp_seq(nbp_seq) = ibp
               endif
            endif
         enddo

         ! Shuffle
         do i = 1, nbp_seq
            !rnd = genrand64_real3()   ! (0,1)-real-interval
            rnd = genrand_double3(mts(irep))   ! (0,1)-real-interval
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
            !rnd = genrand64_real1()  ! [0,1]-real-interval
            rnd = genrand_double1(mts(irep))  ! [0,1]-real-interval

            if (rnd < ratio) then
               ibp_delete = jbp
            endif
         enddo

         ! Delete
         bp_status(ibp_delete, irep) = .False.
         !ene_bp(ibp_delete) = 0.0_PREC
         !! This line is commented out because ene_bp has to be kept for CHECK_FORCE.
         !! In normal run, ene_bp will never be refered when bp_status is False,
         !! thus it does not have to be zero cleared.

         ! Update nt_bp_excess
         nt_bp_excess(bp_mp(1, ibp_delete, irep)) = nt_bp_excess(bp_mp(1, ibp_delete, irep)) - 1
         nt_bp_excess(bp_mp(2, ibp_delete, irep)) = nt_bp_excess(bp_mp(2, ibp_delete, irep)) - 1

      enddo

   endif

   ! Sum of ene_bp masked by bp_status (Note: bp_status(1:nbp_max))
   Ebp = sum(ene_bp(1:nbp(irep)), bp_status(1:nbp(irep), irep))

   if (flg_out_bp) then

      call ieee_get_halting_mode(IEEE_UNDERFLOW, halt_mode)
      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=.false. )

      do ibp = 1, nbp(irep)

         !nhb = bp_type2nhb(bp_mp(3, ibp, irep))
         !if (ene_bp(ibp) < - nhb * kT) then
         if (bp_status(ibp, irep)) then
            imp = bp_mp(1, ibp, irep)
            jmp = bp_mp(2, ibp, irep)
            write(hdl_bp(irep)) int(imp,kind=KIND_OUT_BP), int(jmp,kind=KIND_OUT_BP), &
                                real(ene_bp(ibp), kind=KIND_OUT_BPE)
         endif
      enddo

      write(hdl_bp(irep)) int(0,kind=KIND_OUT_BP), int(0,kind=KIND_OUT_BP), &
                          real(0.0, kind=KIND_OUT_BPE)

      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=halt_mode)
   endif

   if (flg_out_bpall) then

      call ieee_get_halting_mode(IEEE_UNDERFLOW, halt_mode)
      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=.false. )

      do ibp = 1, nbp(irep)

         !if (ene_bp(ibp) < -ZERO_JUDGE) then  ! To output all
         if (bp_status(ibp, irep)) then
            imp = bp_mp(1, ibp, irep)
            jmp = bp_mp(2, ibp, irep)
            write(hdl_bpall(irep)) int(imp,kind=KIND_OUT_BP), int(jmp,kind=KIND_OUT_BP), &
                                   real(ene_bp(ibp), kind=KIND_OUT_BPE)
         endif
      enddo

      write(hdl_bpall(irep)) int(0,kind=KIND_OUT_BP), int(0,kind=KIND_OUT_BP), &
                             real(0.0, kind=KIND_OUT_BPE)

      call ieee_set_halting_mode(IEEE_UNDERFLOW, halting=halt_mode)
   endif

   if (flg_out_bpe) then

      do ibp = 1, nbp(irep)

         !nhb = bp_type2nhb(bp_mp(3, ibp, irep))

         !if (ene_bp(ibp) < - nhb * kT) then
         !if (ene_bp(ibp) < -ZERO_JUDGE) then  ! To output all
         if (bp_status(ibp, irep)) then
            imp = bp_mp(1, ibp, irep)
            jmp = bp_mp(2, ibp, irep)
            write(hdl_bpe(irep), '(1x,i5,1x,i5,1x,f5.2)', advance='no') imp, jmp, ene_bp(ibp)
         endif
      enddo

      write(hdl_bpe(irep), '(a)') ''
   endif

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

end subroutine energy_bp_limit_triplet
