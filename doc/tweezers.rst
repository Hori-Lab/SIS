Tweezers option
===============

``[Tweezers]`` option is to mimic the experiments using optical tweezers
`(wikipedia) <https://en.wikipedia.org/wiki/Optical_tweezers>`__, in
which highly focused lasers are used to grasp and manipulate
micron-scale beads typically attached to biomolecules of interest. See
e.g. Bustamante et al (2021) *Nat. Rev. Methods Primers*
`10.1038/s43586-021-00021-6 <https://doi.org/10.1038/s43586-021-00021-6>`__
for a comprehensive review.

There are two different modes of operation implemented as follows.

Constant Force mode
-------------------

In this mode, it is assumed that one or more pairs of coarse-grained
particles (:math:`i` and :math:`j`) are gripped with optical tweezers
and a constant force :math:`\vec{F}` is applied to each pair in opposite
directions.

.. figure:: images/Tweezers_DCF.png
   :alt: Example set up of Dual_Constant_Force

   Example set up of Dual_Constant_Force

Input file
~~~~~~~~~~

::

   [Tweezers]
       [Tweezers.Dual_Constant_Force]
       id_pairs = [[1, 5],]
       forces_pN = [[5.0, 0.0, 0.0],]

- ``id_pairs`` specfies the particle pairs that are subject of constant
  force. In the example above, particles 1 and 5 will be the subject of
  external force.
- ``forces_pN`` specifies the force vector applied to each pair. In the
  above example, :math:`F_i = (f_x, f_y, f_z) = (5.0, 0, 0)` pN is
  applied to the first particle. The exact opposite force,
  :math:`F_j = (-5.0, 0, 0)` pN, is applied to the particle 5.

Force Ramp mode
---------------

In the ``Force_Ramp`` mode, a pair of particles :math:`i` and :math:`j`
are each trapped by tweezers with harmonic potentials with the spring
constants :math:`k_i` and :math:`k_j`. The positions of the two traps
are initially located at the same position as the particles :math:`i`
and :math:`j` at time = 0 (:math:`\vec{R_i^0}` and :math:`\vec{R_j^0}`).
The trap positions move along the simulation with velocity
:math:`\vec{v_i}` and :math:`\vec{v_j}`, which are defined by the
initial-position vector :math:`\left(\vec{R_j^0} - \vec{R_i^0}\right)`
and given absolute speeds :math:`s_i` and :math:`s_j`.

.. figure:: images/Tweezers_FR.png
   :alt: Example set up of Force_Ramp

   Example set up of Force_Ramp

.. _input-file-1:

Input file
~~~~~~~~~~

::

   [Files]
       [Files.Out]
       types = [...., "twz"]

   [Tweezers]
       [Tweezers.Force_Ramp]
       id_pairs = [[1, 27]]
       spring_const = [[5.0, 5.0],]
       trap_speed = [[0.0, 1.0e-6]]

- ``id_pairs`` specfies the particle pairs (:math:`i` and :math:`j`)
  being trapped.

- ``spring_const`` specifies the spring constant :math:`k` of harmonic
  potential applied to each pair, in the unit of
  :math:`\rm{kcal/mol/Å^2}`.

- ``trap_speed`` specifies the absolute speed of the trap movement. The
  unit is Å per simulation time step.

- In the above example, the pair of particles 1 and 27 are each trapped
  with spring constants :math:`k_1 = k_{27} = 5`
  :math:`\rm{kcal/mol/Å^2}`. The positions of the two traps are
  initially located at the same position as the particles at time =
  :math:`0` (:math:`\vec{R_{1}^0}` and :math:`\vec{R_{27}^0}`). The trap
  positions move along the simulation at the speed given by
  :math:`s_i = 0` and :math:`s_j = 10^{-6}` Å per simulation step. In
  practice, as in this example, the first trap is fixed in the space by
  specifying zero speed.

**DO NOT FORGET to add ``"twz"`` in ``[Files.Out]`` ``types`` to get
information of extension and forces.**

Output files
~~~~~~~~~~~~
