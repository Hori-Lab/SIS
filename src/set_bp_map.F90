subroutine set_bp_map()

   use, intrinsic :: iso_fortran_env, Only : output_unit
   use const
   use const_idx, only : SEQT, seqt2char, seqt2nnt, seqt2bpt
   use var_io, only : iopen_hdl, cfile_prefix
   use var_top, only : nmp, seq, lmp_mp, ichain_mp, nmp_chain
   use var_potential, only : bp_paras, basepair_parameters, bp_model, bp_map, bp_map_0, bp_map_paras, &
                             NN_dG, NN_dH, NN_dS, dH0, dS0, coef_dG
   use var_state, only : tempK, temp_independent

   implicit none

   integer :: imp, jmp
   integer :: i, j, ichain, jchain
   integer :: hdl, istat
   type(basepair_parameters) :: bpp, new_bpp
   real(PREC) :: dG, dH, dS, multi, ln, a, b, c, delta

   if (.not. (bp_model == 4 .or. bp_model == 5)) then
      return
   endif

   bp_map(:,:) = bp_map_0(:,:)

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl
   open(hdl, file=trim(cfile_prefix)//'.bpcoef', status='unknown', action='write', position='append', iostat=istat)
   if (istat /= 0) then
      print '(2a)', 'Error: failed to open bpcoef file. ', trim(cfile_prefix)//'.bpcoef'
      flush(output_unit)
      error stop
   endif

   if (bp_model == 5) then

      if (temp_independent == 0) then
         write(hdl, '(a, f8.3)') '########## tempK = ', tempK
      else
         write(hdl, '(a)') '########## temperature independent'
      endif
   endif

   !bp_map_dG(:,:,:) = 0.0_PREC

   do imp = 1, nmp-1
      i = lmp_mp(imp)
      ichain = ichain_mp(imp)

      ! Either 5' or 3' end
      if (i == 1 .or. i == nmp_chain(ichain)) cycle

      do jmp = imp+1, nmp

         if (bp_map_0(imp, jmp) == 0) cycle

         j = lmp_mp(jmp)
         jchain = ichain_mp(jmp)

         if (bp_model == 4) then

            dG = 0.0_PREC
            if (bp_map_0(imp-1, jmp+1) > 0) then
               dG = 0.5 * NN_dG(seqt2nnt(seq(i-1, ichain), seq(i, ichain), seq(j+1, jchain), seq(j, jchain)))
            endif

            if (bp_map_0(imp+1, jmp-1) > 0) then
               dG = dG + 0.5 * NN_dG(seqt2nnt(seq(i, ichain), seq(i+1, ichain), seq(j, jchain), seq(j-1, jchain)))
            endif

         else if (bp_model == 5) then

            dH = 0.0_PREC
            dS = 0.0_PREC

            if (bp_map_0(imp-1, jmp+1) > 0) then
               dH = 0.5_PREC * (NN_dH(seqt2nnt(seq(i-1, ichain), seq(i, ichain), seq(j+1, jchain), seq(j, jchain))) - dH0)
               dS = 0.5_PREC * (NN_dS(seqt2nnt(seq(i-1, ichain), seq(i, ichain), seq(j+1, jchain), seq(j, jchain))) - dS0)
            endif

            if (bp_map_0(imp+1, jmp-1) > 0) then
               dH = dH + 0.5_PREC * (NN_dH(seqt2nnt(seq(i, ichain), seq(i+1, ichain), seq(j, jchain), seq(j-1, jchain))) - dH0) 
               dS = dS + 0.5_PREC * (NN_dS(seqt2nnt(seq(i, ichain), seq(i+1, ichain), seq(j, jchain), seq(j-1, jchain))) - dS0)
            endif

            ! Default
            if (temp_independent == 0) then
               dG = coef_dG * (dH - tempK * 1.0e-3_PREC * dS)
            else
               dG = coef_dG * dH
            endif

         endif

         if (dG < 0.0_PREC) then
            bpp = bp_paras(seqt2bpt(seq(i,ichain), seq(j,jchain)))
            new_bpp = bpp
            new_bpp%U0 = dG
            if (dG < bpp%U0) then
                !new_k = para_fitted[1]/(2*np.log(new_U0*np.sqrt(math.e)/para_fitted[0]))
                !1/np.sqrt(math.e)=0.6
                multi = 1/(2*log(dG/bpp%U0/0.6_PREC))
                new_bpp%bond_k = bpp%bond_k*multi 
                new_bpp%angl_k1 = bpp%angl_k1*multi
                new_bpp%angl_k2 = bpp%angl_k2*multi
                new_bpp%angl_k3 = bpp%angl_k3*multi
                new_bpp%angl_k4 = bpp%angl_k4*multi 
                !ln=np.log(para_fitted[0]/new_U0)
                !a=1-2*ln
                !b=(2*para_fitted[1]+ln)*(ln-1)
                !c=para_fitted[1]*(para_fitted[1]+ln)
                !delta=b**2-4*a*c
                !new_k = (-b-np.sqrt(delta))/(2*a)
            
                ln=log(bpp%U0/dG)
                a=1.0_PREC-2.0_PREC*ln
                
                b=(2.0_PREC*bpp%dihd_k1+ln)*(ln-1.0_PREC)
                c=bpp%dihd_k1*(bpp%dihd_k1+ln)
                delta=b**2-4.0_PREC*a*c
                new_bpp%dihd_k1 = (-b-sqrt(delta))/(2.0_PREC*a)
                
                b=(2.0_PREC*bpp%dihd_k2+ln)*(ln-1.0_PREC)
                c=bpp%dihd_k2*(bpp%dihd_k2+ln)
                delta=b**2-4.0_PREC*a*c
                new_bpp%dihd_k2 = (-b-sqrt(delta))/(2.0_PREC*a)
            else
                !new_k = 2*para_fitted[1]*np.log(np.sqrt(math.e)*para_fitted[0]/new_U0)
                multi = 2.0_PREC*log(bpp%U0/dG/0.6_PREC)
                new_bpp%bond_k = bpp%bond_k*multi
                new_bpp%angl_k1 = bpp%angl_k1*multi
                new_bpp%angl_k2 = bpp%angl_k2*multi
                new_bpp%angl_k3 = bpp%angl_k3*multi
                new_bpp%angl_k4 = bpp%angl_k4*multi
                !new_k = para_fitted[1]-np.log(new_U0/para_fitted[0])/(1+1/(2*para_fitted[1])-np.sqrt(1+1/((2*para_fitted[1])**2)))
                new_bpp%dihd_k1 = bpp%dihd_k1 - log(dG/bpp%U0)/(1.0_PREC+1.0_PREC/(2.0_PREC*bpp%dihd_k1)-sqrt(1.0_PREC+1.0_PREC/((2.0_PREC*bpp%dihd_k1)**2)))
                new_bpp%dihd_k2 = bpp%dihd_k2 - log(dG/bpp%U0)/(1.0_PREC+1.0_PREC/(2.0_PREC*bpp%dihd_k2)-sqrt(1.0_PREC+1.0_PREC/((2.0_PREC*bpp%dihd_k2)**2)))
            endif

            bp_map_paras(imp, jmp) = new_bpp
            bp_map_paras(jmp, imp) = new_bpp

            write(hdl, '(i5,1x,i5,3x,7a1,3x,f8.3,3x,f6.2,3x,f6.2,3x,f6.2,3x,f6.2,3x,f6.2,3x,f6.2,3x,f6.2)') imp, jmp, &
                      seqt2char(seq(i-1,ichain)), seqt2char(seq(i,ichain)), seqt2char(seq(i+1,ichain)), '/', &
                      seqt2char(seq(j+1,jchain)), seqt2char(seq(j,jchain)), seqt2char(seq(j-1,jchain)), new_bpp%U0, &
                      new_bpp%bond_k, new_bpp%angl_k1, new_bpp%angl_k2, new_bpp%angl_k3, new_bpp%angl_k4, &
                      new_bpp%dihd_k1, new_bpp%dihd_k2

         else
            bp_map(imp, jmp) = 0
            bp_map(jmp, imp) = 0

         endif

      enddo
   enddo

   close(hdl)
   iopen_hdl = iopen_hdl - 1

endsubroutine set_bp_map
