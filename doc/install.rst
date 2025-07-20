Getting started
================

.. _installation-1:

Installation
------------

The code is written in Fortran and the build environment uses CMake,
therefore a Fortran compiler and ``cmake`` are needed to compile.

The build command should be something like,

::

   $ mkdir build
   $ FC=gfortran cmake ..
   $ cd ../

A convenient build script is provided as ``build.sh``. In many cases,
the user can just execute the script from the base directory.

Running a simulation
--------------------

`SIS` is a command-line program (sorry, no GUI). 

The simplest way to run the program is,

::

   $ ./sis input.toml

`input.toml` is an input file.


It is **highly recommended to redirect the standard output and error to files**,

::

   $ ./sis input.toml 1> out 2> err



Restarting a simulation
~~~~~~~~~~~~~~~~~~~~~~~

End of the simulation
~~~~~~~~~~~~~~~~~~~~~

There are three wasy of a simulation stops.

1. Reaching the number of simulation steps specified in the input file.

::

   [MD]
   nstep = 500000

2. Reaching the wall time specified in the input file.

(optional, default: -1) Wall-clock time limit for the job. The default
is -1 thus no time limit. If a positive value is specified, the job will
be terminated once the wall time reaches the limit value. A restart file
(.rst) will be generated at the time. The unit is hour and can be any
positive real number, e.g. 0.25 for 15 minutes. If the specified value
is negative, no limit will be set. It is recommended to set a slightly
smaller value than the resource time allowed, Since it may take some
time for the program to exit normally, such as when saving files, it is
recommended that you specify a slightly smaller value, e.g. 47.9 for
48-hour job.

::

   [MD]
   stop_wall_time_hour = 47.9

   nstep_check_stop = 10000
       # (optional, default: 10000) Frequency to check the wall time to stop. 
       # The value has to be a positive integer. If a value less than 1 is specified, 
       # then nstep_check_stop = 1 will be used (check evey step).
       # The simulation can also be stopped by placing an empty file named STOP_SIS at the same directory
       # specified in [Files.In] prefix.
       # e.g. if prefix = "./output/md", either ./output/stop_sis or ./output/STOP_SIS will terminate the job.

1. The program detects a STOP file.
