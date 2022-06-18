program sisbp_Q

   use, intrinsic :: iso_fortran_env, Only : iostat_end
   implicit none

   integer, parameter :: PREC = 8
   integer(4) :: kind_int, kind_real
   integer :: i, j, n
   integer :: nnt
   integer :: istat
   integer :: idummy
   integer :: hdl_bp, hdl_ref, hdl_out
   integer :: nbp_ref
   real(PREC) :: e, q
   real(PREC) :: e_cutoff
   logical :: flg_ct
   logical, allocatable :: ref(:,:), map(:,:)
   character(1) :: cdummy
   character(64) :: cpara
   character(len=1000) :: cfile_bp, cfile_ref, cfile_out

   if (command_argument_count() == 4) then
      
      call get_command_argument(1, cfile_bp)  
      call get_command_argument(2, cfile_ref)  
      call get_command_argument(3, cpara)  
      read(cpara, *) e_cutoff
      call get_command_argument(4, cfile_out)  

   else
      print '(a)', 'Usage: PROGRAM (input bp file) (basepairs .ct or .bpseq file) (Ebp cutoff) (output Q file)'
      print '(a)', '  (e.g.,  a.out md.bp reference.ct -0.614 q.out)'
      stop
   endif

   n = len(trim(cfile_ref))
   if (cfile_ref(n-2:n) == '.ct') then
      flg_ct = .True.
   else if (cfile_ref(n-5:n) == '.bpseq') then
      flg_ct = .False.
   else
      print *, 'Error: reference file has to be either .ct or .bpseq format'
      stop
   endif

   open(newunit=hdl_ref, file=cfile_ref, status='old', action='read')

   if (flg_ct) then
      read(hdl_ref, *) nnt

      allocate(ref(nnt, nnt))
      ref(:,:) = .false.
      allocate(map(nnt, nnt))
      map(:,:) = .false.

      do
         read(hdl_ref, *, iostat=istat) i, cdummy, idummy, idummy, j, idummy
         if (istat == iostat_end) exit

         if (i > j) then
            idummy = i
            i = j
            j = idummy
         endif

         ref(i,j) = .true.
      enddo

   else
      nnt = 0
      do
         read(hdl_ref, *, iostat=istat) i, cdummy, j
         if (istat == iostat_end) exit
         nnt = nnt + 1
      enddo

      allocate(ref(nnt, nnt))
      ref(:,:) = .false.
      allocate(map(nnt, nnt))
      map(:,:) = .false.

      rewind(hdl_ref)

      do
         read(hdl_ref, *, iostat=istat) i, cdummy, j
         if (istat == iostat_end) exit

         if (i > j) then
            idummy = i
            i = j
            j = idummy
         endif

         ref(i,j) = .true.
      enddo

   endif

   close(hdl_ref)

   nbp_ref = count(ref)

   open(newunit=hdl_bp, file=cfile_bp, status='old', action='read', form='unformatted',access='stream')
   open(newunit=hdl_out, file=cfile_out, status='unknown', action='write')

   read(hdl_bp) kind_int
   read(hdl_bp) kind_real

   loop_read: do while (.True.)

      map(:,:) = ref(:,:)

      loop_frame: do

         call read_bp(hdl_bp, i, j, e, istat)
         if (istat == iostat_end) then
            exit loop_read
         endif

         if (i == 0 .and. j == 0) then
            exit loop_frame
         endif

         if (e < e_cutoff) then
            map(i,j) = .false.
         endif

      enddo loop_frame

      q = (nbp_ref - count(map)) / real(nbp_ref, kind=PREC)

      write(hdl_out, '(f8.6)') q
   enddo loop_read

   close(hdl_bp)
   close(hdl_out)


contains
   subroutine read_bp(hdl, i, j, e, istat)

      integer, intent(in) :: hdl
      integer, intent(out) :: i, j
      real(PREC), intent(out) :: e
      integer(2) :: i2, j2
      real(4) :: e4
      integer, intent(out) :: istat

      if (kind_int == 2 .and. kind_real == 4) then
         read(hdl, iostat=istat) i2, j2, e4
         i = int(i2)
         j = int(j2)
         e = real(e4, kind=PREC)
      endif

   endsubroutine read_bp

endprogram sisbp_Q
