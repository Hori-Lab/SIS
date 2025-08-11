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
   $ cd build/
   $ FC=gfortran cmake ..
   $ cd ../

A convenient build script is provided as ``build.sh``. In many cases,
the user can just execute the script from the base directory.
By default, ``gfortran`` will be used as a compiler.

::

   $ ./build.sh

Running a simulation
--------------------

`SIS` is a command-line program. The simplest way to run the program is,

::

   $ ./sis input.toml

`input.toml` is an input file.


It is, however, **highly recommended to redirect the standard output and error to files**,

::

   $ ./sis input.toml 1> out 2> err



Restarting a simulation
~~~~~~~~~~~~~~~~~~~~~~~

End of the simulation
~~~~~~~~~~~~~~~~~~~~~

There are three ways of a simulation stops.

1. Reaching the number of simulation steps specified in the input file.

::

   [MD]
   nstep = 500000

2. Reaching the wall time specified in the input file.
   
User can specify the wall-clock time limit for the job by :code:`stop_wall_time_hour` in the :code:`[MD]` section. The default is -1,
which means no time limit. If a positive value is specified, the job will
be terminated once the wall time reaches the limit value. A restart file
(.rst) will be generated at the time. The unit is hour and can be any
positive real number, e.g. 0.25 for 15 minutes. Since it may take some time for the program to exit normally saving files, it is recommended that the user specify a slightly shorter time, e.g. 47.9 for 48-hour job.

::

   [MD]
   stop_wall_time_hour = 47.9
   nstep_check_stop = 10000

:code:`nstep_check_stop` is the frequency to check the wall time to stop. The value has to be a positive integer. If this is not specified in the input file, then the default value of 10,000 is used. If a value less than 1 is specified, then :code:`nstep_check_stop = 1` will be used (check evey step).

3. The program detects a STOP_SIS file.

The simulation can also be stopped by placing an empty file named STOP_SIS at the same directory specified in :code:`[Files.In] prefix`. 
For example, if :code:`prefix = "./output/md"`, either :code:`./output/stop_sis` or :code:`./output/STOP_SIS` will terminate the job. This is a way to stop the simulation safely from the outside, ensuring that the output files are saved properly. 
