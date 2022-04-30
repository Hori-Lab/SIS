subroutine list_exv()

   use const
   use const_idx, only : SEQT
   use var_top, only : nmp_chain, imp_chain, nchains
   use var_potential, only : nwca, wca_mp

   implicit none
  
   integer :: ichain, jchain
   integer :: n, i, j, j_start, imp, jmp
   integer :: iwca

   if (allocated(wca_mp)) then
      deallocate(wca_mp)
   endif

   do n = 1, 2
   
      iwca = 0

      do ichain = 1, nchains 

         do jchain = ichain, nchains 

            do i = 1, nmp_chain(ichain)

               imp = imp_chain(i, ichain)
    
               if (ichain == jchain) then
                  j_start = i + 3
               else
                  j_start = 1
               endif
    
               do j = j_start, nmp_chain(jchain)
                  
                  jmp = imp_chain(j, jchain)
    
                  iwca = iwca + 1

                  if (n == 2) then
                     wca_mp(1,iwca) = imp
                     wca_mp(2,iwca) = jmp
                  endif
    
               enddo
            enddo
    
         enddo
      enddo

      if (n == 1) then
         nwca = iwca
         allocate(wca_mp(2, nwca))
      endif
   enddo

   write(*,*) '#nwca: ', nwca

end subroutine list_exv
