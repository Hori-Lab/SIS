subroutine init_out_files

   use const, only : CHAR_FILE_PATH, FILENAME_DIGIT_REPLICA
   use const_idx, only : JOBT
   use var_state, only : job
   use var_io, only : iopen_hdl, cfile_prefix, KIND_OUT_BP, KIND_OUT_BPE, &
                      flg_out_dcd, flg_out_rst, flg_out_bp, flg_out_bpall, flg_out_bpe, &
                      hdl_out, hdl_dcd, hdl_rst, hdl_bp, hdl_bpall, hdl_bpe, hdl_rep, &
                      cfile_rst, cfile_dcd
   use var_replica, only : nrep_proc, flg_replica, irep2grep
   use var_parallel, only : myrank

   implicit none

   integer :: irep
   character(20) :: crep = ''
   character(len=:), allocatable :: cfilebase
   !character(len=CHAR_FILE_PATH) :: cfilename
   character(len=:), allocatable  :: cfilename

   if (job == JOBT%CHECK_FORCE) then
      flg_out_bp = .False.
      flg_out_bpall = .False.
      flg_out_bpe = .False.
   endif

   if (job == JOBT%MD) then
      flg_out_dcd = .True.
      flg_out_rst = .True.
   endif


   !! Allocate
   allocate(hdl_out(nrep_proc))

   !! rst
   if (flg_out_rst) then
      allocate(hdl_rst(nrep_proc))
      allocate(cfile_rst(nrep_proc))
   endif

   !! dcd
   if (flg_out_dcd) then
      allocate(hdl_dcd(nrep_proc))
      allocate(cfile_dcd(nrep_proc))
   endif

   if (flg_out_bp) allocate(hdl_bp(nrep_proc))
   if (flg_out_bpall) allocate(hdl_bpall(nrep_proc))
   if (flg_out_bpe) allocate(hdl_bpe(nrep_proc))


   !! Generate filenames and open
   do irep = 1, nrep_proc

      if (flg_replica) then
         write(crep, '(i20)') irep2grep(irep) + 10**FILENAME_DIGIT_REPLICA 
         cfilebase = trim(cfile_prefix) // '_' // crep(21-FILENAME_DIGIT_REPLICA:20)
      else
         cfilebase = trim(cfile_prefix)
      endif

      !! out
      cfilename = trim(cfilebase) // '.out'
      iopen_hdl = iopen_hdl + 1
      hdl_out(irep) = iopen_hdl
      open(hdl_out(irep), file = cfilename, status = 'replace', action = 'write', form='formatted')

      !! dcd  (The file will be opened in job_md using dcd module.)
      if (flg_out_dcd) then
         iopen_hdl = iopen_hdl + 1
         hdl_dcd(irep) = iopen_hdl
         cfile_dcd(irep) = cfilebase // '.dcd'
      endif

      !! rst  (The file will be opened each time of writing.)
      if (flg_out_rst) then
         iopen_hdl = iopen_hdl + 1
         hdl_rst(irep) = iopen_hdl
         cfile_rst(irep) = trim(cfilebase) // '.rst'
      endif

      !! bp
      if (flg_out_bp) then
         cfilename = trim(cfilebase) // '.bp'
         iopen_hdl = iopen_hdl + 1
         hdl_bp(irep) = iopen_hdl
         open(hdl_bp(irep), file=cfilename, status='replace', action='write', form='unformatted',access='stream')
         write(hdl_bp(irep)) int(KIND_OUT_BP,kind=4)
         write(hdl_bp(irep)) int(KIND_OUT_BPE,kind=4)
      endif

      !! bpall
      if (flg_out_bpall) then
         cfilename = trim(cfilebase) // '.bpall'
         iopen_hdl = iopen_hdl + 1
         hdl_bpall(irep) = iopen_hdl
         open(hdl_bpall(irep), file=cfilename, status='replace', action='write', form='unformatted',access='stream')
         write(hdl_bpall(irep)) int(KIND_OUT_BP,kind=4)
         write(hdl_bpall(irep)) int(KIND_OUT_BPE,kind=4)
      endif

      !! bpe
      if (flg_out_bpe) then
         cfilename = trim(cfilebase) // '.bpe'
         iopen_hdl = iopen_hdl + 1
         hdl_bpe(irep) = iopen_hdl
         open(hdl_bpe(irep), file=cfilename, status='replace', action='write', form='formatted')
      endif
   enddo

   if (flg_replica .and. myrank == 0) then
      cfilename = trim(cfile_prefix) // '.rep'
      iopen_hdl = iopen_hdl + 1
      hdl_rep = iopen_hdl
      open(hdl_rep, file=cfilename, status='replace', action='write', form='formatted')
   endif

endsubroutine init_out_files
