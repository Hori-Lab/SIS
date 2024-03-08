subroutine write_rep_table()

   use const
   use const_idx, only : REPT
   use const_phys,only : JOUL2KCAL_MOL
   use var_io,     only : hdl_rep
   use var_replica, only : flg_repvar, rep2val, lab2val, nrep_all
   use var_parallel, only : myrank

   implicit none
   integer     :: irep, ivar

   if (myrank == 0) then

      write(hdl_rep,'(a)') '# Table of replica valiable and label'
      write(hdl_rep,'(a)',ADVANCE="NO") '#label '
      do ivar = 1, REPT%MAX
         if (flg_repvar(ivar)) then
            select case (ivar)
            case (REPT%TEMP)
               write(hdl_rep, '(a)', ADVANCE = "NO") 'Temperature   '
            !case (REPT%ION)
            !   write(hdl_rep, '(a)', ADVANCE = "NO") 'IonicStrength '
            case (REPT%TWZDCF)
               write(hdl_rep, '(a)', ADVANCE = "NO") 'External-force '
            endselect
         endif
      enddo
      write(hdl_rep,*) ''

      do irep = 1, nrep_all
         write(hdl_rep, '(i5,1x)', ADVANCE="NO") irep
         do ivar = 1, REPT%MAX
            if (flg_repvar(ivar)) then
               if (ivar == REPT%TWZDCF) then
                  write(hdl_rep, '(f12.4,1x)', ADVANCE="NO") lab2val(irep, ivar) / (JOUL2KCAL_MOL * 1.0e-22)
               else
                  write(hdl_rep, '(f12.4,1x)', ADVANCE="NO") lab2val(irep, ivar)
               endif
            endif
         enddo
         write(hdl_rep, *) ''
      enddo
      write(hdl_rep, '(a)') ''
      write(hdl_rep, '(a)') ''

      write(hdl_rep, '(a)') '# History of replica labels'

      flush(hdl_rep)
   endif

endsubroutine write_rep_table
