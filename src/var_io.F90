module var_io
   
   use const
   use const_idx, only : ENE

   logical :: flg_out_bpcoef = .False.
   logical :: flg_out_bp = .False.
   logical :: flg_out_bpall = .False.
   logical :: flg_out_bpe = .False.
   logical :: flg_out_T = .False.
   logical :: flg_out_dcd = .False.
   logical :: flg_out_rst = .False.
   logical :: flg_out_twz = .False.

   logical :: flg_in_bpseq = .False.
   logical :: flg_in_ct = .False.
   logical :: flg_in_fasta = .False.
   logical :: flg_in_xyz = .False.
   logical :: flg_in_pdb = .False.

   logical :: flg_gen_init_struct = .False.

   integer, parameter :: KIND_OUT_BP  = 2  ! Defines the format of bp output file.
   integer, parameter :: KIND_OUT_BPE = 4  ! Defines the format of bp output file.

   integer, parameter :: hdl_in_rst = 10  ! input rst

   integer, allocatable :: hdl_out(:)
   integer, allocatable :: hdl_bpcoef(:)
   integer, allocatable :: hdl_dcd(:)
   integer, allocatable :: hdl_rst(:)
   integer, allocatable :: hdl_bp(:)
   integer, allocatable :: hdl_bpall(:)
   integer, allocatable :: hdl_bpe(:)
   integer  :: hdl_rep
   integer  :: hdl_twz

   integer :: iopen_hdl = 15

#ifdef DUMPFORCE
   integer :: hdl_force(1:ENE%MAX)
#endif

   logical :: flg_progress
   integer :: step_progress

   !character(len=CHAR_FILE_PATH) :: cfile_out
   character(len=CHAR_FILE_PATH), allocatable :: cfile_rst(:)  ! (nrep_proc)
   character(len=CHAR_FILE_PATH), allocatable :: cfile_dcd(:)  ! (nrep_proc)
   character(len=:), allocatable :: cfile_ff, cfile_prefix, cfile_dcd_in, cfile_pdb_ini, &
                                    cfile_fasta_in, cfile_ct_in, cfile_bpseq_in, &
                                    cfile_anneal_in, cfile_xyz_ini

end module var_io
