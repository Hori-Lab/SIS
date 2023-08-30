# How to run replica exchange simulations.

## Input file

```
[Replica]
nrep_temp = 14
nstep_exchange = 100
nstep_save = 10000

   [Replica.Temperature]
   1 = 273.15
   2 = 283.15
   3 = 293.15
   4 = 303.15
   5 = 313.15
   6 = 323.15
   7 = 328.15
   8 = 333.15
   9 = 338.15
   10 = 343.15
   11 = 348.15
   12 = 353.15
   13 = 363.15
   14 = 373.15
```

`nrep_temp` is the number of replicas. In this example, there will be 14 replcias. Therfore, you need to specify eight temperatures in `[Replica.Temperature]`.

`nstep_exchange` is the frequency of attempting replica exchanges.

`nstep_save` is the frequency of saving replica states into `.rep` file. I recommend that you set `nstep_save` to be the same as `[MD] nstep_save`.

## Submitting a job

### Serial job
If you run a serial job, all the replicas will be simulated in a single CPU. For this, you can just submit the job as usual. Each step of simulations will be ran for all replicas sequentially, so it will be significantly slower (e.g. at least 14 times slower than a normal simulation if there are 14 replicas).

### MPI job

First, you need to compile the code with MPI. In case of Sulis, ` module load iimpi` in addition to the `intel` module. This can be written in `~/.bash_profile`.

You can find a typical compile command in `build_mpi.sh`, or simply execute the script. It will compile the code for you and generate `sis.mpi`. 

To submit an MPI job in Sulis, add the following into your submission batch file.

```
   ....
#SBATCH --nodes=1
#SBATCH --ntasks=14
#SBATCH --mem-per-cpu=500
   ...
   
module purge
module load GCCcore/11.3.0 intel iimpi

srun mpirun.sh 1> out 2> err
```

Separately, you need `mpirun.sh` with the following.

``` 
#!/bin/sh

# Uncomment if OpenMPI
#rank=$OMPI_COMM_WORLD_RANK

# IntelMPI (Sulis)
rank=$PMI_RANK

suffix=`printf "%3.3d" $rank`

../sis.mpi input.toml 1> ./out.$suffix 2> ./err.$suffix

## In case restarting simulation
#../sis.mpi input.toml run1.rst 1> ./out.$suffix 2> ./err.$suffix
```

`../sis.mpi` should be replaced by the path of your execute file.

## Output files

Output files of a replica simulation will be something like `md_0001.out`, `md_0002.out`, `md_0003.out` and so on, where the prefix is the one specified in the input file, `[Files.Out] prefix = md`.

In each of these output files, a trajectory is recoreded for a replica specified by the ID. For example, `md_0001.out` contains the energy trajectory of Replica 1, and `md_0001.dcd` is for the coordinates of Replica 1. 

Here is an example of such output file.

```
$ head -n 5 md_0001.out
#(1)nframe (2)R (3)T   (4)Ekin       (5)Epot       (6)Ebond      (7)Eangl      (8)Edih       (9)Ebp        (10)Eexv      (11)Eele
         0    1 273.15   33.4404      -2.85342       1.51499      0.185399      -10.4430       0.00000       0.00000       5.88923
    100000    3 293.15   28.6941      -5.66499       5.66644       4.35206      -22.2032     -0.129127E-07   0.00000       6.51973
    200000    9 338.15   30.1093       10.1126       9.88961       9.22375      -16.0226     -0.860057E-11   0.00000       7.02182
    300000   12 353.15   29.5965       7.05127       9.49548       8.64525      -18.5929     -0.134088E-10   0.00000       7.50339
```

The second column (R) indicates the label of replica variable(s) (temperature in this example) at each time step. The labels corresponds to the ones specified in the input file.

```
   1 = 273.15
   2 = 283.15
   3 = 293.15
   4 = 303.15
   5 = 313.15
   6 = 323.15
   ....
```

In the example above, Replica 1 was at T = 273.15 K (label 1) at time step = 0, at T = 293.15 K (label 3) at time step = 100000, and then T = 338.15 K (label 9) at time step 200000.

## Analysis

It is often convenient to convert these trajectory files into files correspond to each replica variables (labels). Use a python script, `utils/replica_to_label.py`, to do that.

To show the usage, run with `-h`.

```
$ ./replica_to_label.py -h
usage: replica_to_label.py [-h] indir name nrep outdir

Convert REMD simulation result from files for each replica to files for each label.

positional arguments:
  indir       Input directory path
  name        Filename prefix, e.g. "md" for md_0001.out, md_0002.out ...
  nrep        Number of replicas
  outdir      Output directory path

optional arguments:
  -h, --help  show this help message and exit
```

Here is an example.

```
$ util/replica_to_label.py ./REMD_trajs/ md 14 ./REMD_analysis/
```

The first argument `./REMD_trajs/` is the directory where the replica simulation restuls are stored. The second argument `md` is the prefix of the simulation files. The next, `14`, is the number of replicas. The last argument `./REMD_analysis/` is the directory where the new files will be stored.

This will creat a new set of files under `REMD_analysis/`.

```
$ head -n 5 ./REMD_analysis/md_0001.out
#(1)nframe (2)R (3)T   (4)Ekin       (5)Epot       (6)Ebond      (7)Eangl      (8)Edih       (9)Ebp        (10)Eexv      (11)Eele
         0    1 273.15   33.4404      -2.85342       1.51499      0.185399      -10.4430       0.00000       0.00000       5.88923
    100000    1 273.15   23.4387       3.72253       8.91529       5.74424      -17.1527     -0.252869E-12   0.00000       6.21570
    200000    1 273.15   28.8156       4.49741       10.0573       5.35580      -17.5358     -0.845919E-17   0.00000       6.62015
    300000    1 273.15   20.3474       1.77313       11.6629       5.16187      -21.3556     -0.496145E-38   0.00000       6.30401
```

As seen above, each file now contains a trajecotry for a specifoc replica label.

These files can be used for e.g. WHAM analysis to calculate the melting curve, in the same ways as normal constant-temperature simulations.

