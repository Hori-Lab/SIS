subroutine read_sisinfo(cfilepath)
      
   use const
   use var_potential, only : nbond, bond_mp, bond_para, nangl, angl_mp, angl_para, ndihedral, dihedral_mp, dihedral_para
   use var_io, only : iopen_hdl
  
   implicit none

   character(len=*), intent(in) :: cfilepath
  
   character(len=CHAR_FILE_LINE) :: cline
   integer :: istat

   integer :: hdl
   integer :: s, imp1, imp2, imp3, imp4, iunit1, iunit2, iunit3, imp1_unit, imp2_unit, imp3_unit, iunit4, imp4_unit
   integer :: ibd, iangl, idihedral
   real(PREC) :: k, l

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   ! Count the nubmer of interactions
   ibd = 0
   iangl = 0
   idihedral = 0
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

      ! angl
      if(cline(1:4) == 'angl') then
         iangl = iangl + 1
      endif
      
      !Dihedral
      if(cline(1:4) == 'Dihe') then
         idihedral = idihedral + 1
      endif
   enddo

   nbond = ibd
   nangl = iangl
   ndihedral = idihedral
   ! Allocate
   allocate(bond_mp(2, nbond))
   allocate(bond_para(2, nbond))
   allocate(angl_mp(3, nangl))
   allocate(angl_para(2, nangl))
   allocate(dihedral_mp(3, ndihedral))
   allocate(dihedral_para(2, ndihedral))

   rewind(hdl)

   ibd = 0
   iangl = 0
   idihedral = 0
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
         read(cline(5:), *) s, imp1, imp2, iunit1, imp1_unit, iunit2, imp2_unit, k, l
          
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

      ! ------------------------------------------------------------------
      ! angl
      else if(cline(1:4) == 'angl') then
         read(cline(5:), *) s, imp1, imp2, imp3, iunit1, imp1_unit, iunit2, imp2_unit, iunit3, imp3_unit, k, l
          
         iangl = iangl + 1

         if (imp1 < imp3) then
            angl_mp(1, iangl) = imp1
            angl_mp(2, iangl) = imp2
            angl_mp(3, iangl) = imp3
         else
            angl_mp(1, iangl) = imp3
            angl_mp(2, iangl) = imp2
            angl_mp(3, iangl) = imp1
         endif

         angl_para(1, iangl) = k
         angl_para(2, iangl) = l

      endif
      !--------------------------------------------------------------------------------
      ! dihedral
      else if(cline(1:4) == 'Dihe') then
         read(cline(5:), *) s, imp1, imp2, imp3, imp4, iunit1, imp1_unit, iunit2, imp2_unit, iunit3, imp3_unit, iunit4, imp4_unit, k, l
          
         idihedral = idihedral + 1

         if (imp1 < imp4) then
            dihedral_mp(1, idihedral) = imp1
            dihedral_mp(2, idihedral) = imp2
            dihedral_mp(3, idihedral) = imp3
            dihedral_mp(4, idihedral) = imp4
         else
            dihedral_mp(1, idihedral) = imp4
            dihedral_mp(2, idihedral) = imp3
            dihedral_mp(3, idihedral) = imp2
            dihedral_mp(4, idihedral) = imp1
         endif

         dihedral_para(1, idihedral) = k
         dihedral_para(2, idihedral) = l

      endif

   enddo

   close(hdl)
   iopen_hdl = iopen_hdl - 1

end subroutine read_sisinfo
