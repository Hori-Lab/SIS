Files
================

Input file should be in the [TOML](https://toml.io/) format.

Input files
^^^^^^^^^^^

Minimal input files for a simple RNA simulation

::

    [Files]
    [Files.In]
    ff = "./forcefield.ff"
    fasta = "./sequence.fasta"
    
    pdb_ini = "initial.pdb"
    # xyz_ini = "initial.xyz"


Force field (.ff) file
-----------------------

The force field file is a text file that contains the force field parameters for the UNISIS model. 


FASTA (.fasta) file
-------------------

The FASTA file is a text file that contains the sequence of the RNA strands. This is required and only the source of information for the sequence, i.e. the sequence informaiton is not read from PDB files, for example. All the sequences in the simulation system must be specified in a single FASTA file in the order of the sequence appearance in coordinate files.

Multiple sequences can be specified in the FASTA file. Each sequence block should start from a single line with `>` followed by the sequence consists of A, U, G, C, and D, where D represents a dummy bead.

Below is an example of a valid FASTA file containing three sequences. Note that any empty lines are simply ignored. 
::

    >sequence1
    AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

    >sequence2
    AUGCGAUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGC
    UAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCUAGCU
    >sequence3
    DUAGGCUACGGAUCGUACGGAUCCGAUUGCCGAUCGGA
    UCCGAUUGCCGAUCGGAUCCGAUUGCCGAUCGGAUCCG
    AUUGCCGAUCGGAUCCGAUUGCCGAUCGGAUCCGAUUGCCGAUCGGAUCCGAUUGCCGAUD

Initial coordinate file (.pdb)
------------------------------

The initial coordinate file is a PDB file that contains the initial coordinates of the simulation system. The format of the PDB file is as follows:



Input files to specify secondary structures
-------------------------------------------

::

    [Files]
        [Files.In]
            ...

        ct = "./test.ct"
        # CT (Connectivity Table) file
        bpseq = "./test.bpseq"
        # BPSEQ file (no comment line allowed)
        bpl = "./test.bpl"
        # Basepair list file (no comment line allowed)



Other optional files for specific purposes
------------------------------------------

::

    [Files]
        [Files.In]

        ...

        dcd = "./trajectory.dcd"

        anneal = "./annealing_schedule.txt"

        bias_Rg = "./bias_Rg_schedule.txt"

        restraint = "./restraint.txt"
