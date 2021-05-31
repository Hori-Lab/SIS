subroutine read_sisinfo(cfilepath)
      
   use const
   use var_potential, only : nbond, bond_mp, bond_para
   use var_io, only : iopen_hdl
  
   implicit none

   character(len=CHAR_FILE_PATH), intent(in) :: cfilepath
  
   character(len=CHAR_SISINFO_LINE) :: cline
   integer :: istat

   integer :: hdl
   integer :: s, imp1, imp2, iunit1, iunit2, imp1_unit, imp2_unit, ibd
   real(PREC) :: k, l

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   ! Count the nubmer of interactions
   ibd = 0
   do
      read (hdl, '(a)', iostat = istat) cline
      if(istat < 0) then
         exit
      else if(istat > 0) then
         write(*,*) 'Error: cannot read sisinfo file'
         stop (2)
      end if
 
      ! ------------------------------------------------------------------
      ! bond
      if(cline(1:4) == 'bond') then
         ibd = ibd + 1
      endif
   enddo

   nbond = ibd

   ! Allocate
   allocate(bond_mp(2, nbond))
   allocate(bond_para(2, nbond))

   rewind(hdl)

   ibd = 0
   do
      read (hdl, '(a)', iostat = istat) cline
      if(istat < 0) then
         exit
      else if(istat > 0) then
         write(*,*) 'Error: cannot read sisinfo file'
         stop (2)
      end if
 
      ! ------------------------------------------------------------------
      ! bond
      if(cline(1:4) == 'bond') then
         read(cline(5:), *) s, imp1, imp2, iunit1, iunit2, imp1_unit, imp2_unit, k, l
          
         ibd = ibd + 1

         if (imp1 < imp2) then
            bond_mp(1, ibd) = imp1
            bond_mp(2, ibd) = imp2
         else
            bond_mp(1, ibd) = imp2
            bond_mp(2, ibd) = imp1
         endif

         bond_para(1, ibd) = k
         bond_para(2, ibd) = l

         !write(*,*) bond_mp(1, ibd), bond_mp(2, ibd), bond_para(1, ibd), bond_para(2, ibd)

      endif

   enddo

   close(hdl)
   iopen_hdl = iopen_hdl - 1

end subroutine read_sisinfo
