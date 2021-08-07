subroutine read_pdb(cfilepath, nmp, xyz)

   use const
   !use const_idx
   use var_io, only : iopen_hdl

   implicit none

   character(len=CHAR_FILE_PATH), intent(in) :: cfilepath
   integer, intent(in) :: nmp
   real(PREC), intent(out) :: xyz(3, nmp)

   ! ---------------------------------------------------------------------
   integer :: imp, iatom!, iunit, ires
   integer :: iresnum!, iresnum_save
   !logical :: flg_reading
   real(PREC) :: x, y, z
   real(PREC) :: tempfactor, occupancy  ! These variables are just read, not used.
   character(1) :: multistruct
   character(4) :: nameofatom
   character(6) :: nameid
   character(3) :: nameofmp
   character(2) :: chainid, chainid_save

   integer :: istat, hdl
   character(72) :: char72

   imp = 0
   !iunit = 0
   !ires = 0

   iopen_hdl = iopen_hdl + 1
   hdl = iopen_hdl

   open(hdl, file=cfilepath, status='old', action='read', iostat=istat)

   !flg_reading = .False.

   do
      read (hdl, '(a72)', iostat = istat) char72

      if (istat < 0) then
         exit

      else if (istat > 0) then
         write(*,*) 'Error: cannot read pdb file'
         stop (2)

      end if

      if (char72(1:4) == 'ATOM' .or. char72(1:6) == 'HETATM') then
         read (char72, '(a6, i5, 1x, a4, a1, a3, a2, i4, 4x, 3f8.3, 2f6.2)', iostat = istat) &
              nameid, iatom, nameofatom, multistruct, nameofmp, &
              chainid, iresnum, x, y, z, occupancy, tempfactor
         if (istat > 0) then
            write(*,*) 'Error: cannot read an ATOM line in pdb file'
            stop (2)
         end if

!         if (.not. flg_reading) then
!            iunit = iunit + 1
!            !ires = ires + 1
!            chainid_save = chainid
!            !iresnum_save = iresnum
!            !flg_reading = .true.
!
!         else
!            if (chainid /= chainid_save) then
!               write(*,*) 'Error: chainid /= chainid_save in read_pdb'
!               write(*,*) 'The line causing the error is:'
!               write(*,*) char72
!               stop (2)
!            endif
!            !if (iresnum /= iresnum_save) then
!               !ires = ires + 1
!               !iresnum_save = iresnum
!            !endif
!         end if

         imp = imp + 1
         xyz(1, imp) = x
         xyz(2, imp) = y
         xyz(3, imp) = z

      end if
   end do

   !nunit = iunit
   !nmp   = imp
   !nres  = ires

end subroutine read_pdb
