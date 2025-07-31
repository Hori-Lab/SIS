Theory
======

Single-Interaction-Site RNA model
---------------------------------

[Nguyen2022]_

The UNISIS model is a coarse-grained model for RNA strands. It is a single-interaction-site model, where each nucleotide is represented by a single site.

.. figure:: images/unisi_model.png
   :width: 50%
   :align: center
   :alt: UNISIS model

   The UNISIS model.

The interactions between the sites are given by the following potential energy terms.

.. math::

   U_{\textrm{SIS}} = U_{\textrm{Bond}} + U_{\textrm{Angle}} + U_{\textrm{Dihedral}} +U_{EV} + U_{\textrm{BP}} + U_{\textrm{Ele}}

Local energy terms
^^^^^^^^^^^^^^^^^^^

.. math::

   U_{\textrm{Bond}}=\frac{1}{2}\sum_{\textrm{bond}}{k_{\textrm{b}} (b - b_0)^2} .

where :math:`k_{\textrm{b}}` is the force constant, :math:`b_0` is the equilibrium bond length, and :math:`b` is the bond length.

For consecutive three nucleotides, the angle potential is defined as follows. This *restricted bending* potential prevents the angle approach to 180 degrees [Bulacu2013]_.

.. math::

   U_{\textrm{Angle}} = \frac{1}{2}\sum_{\textrm{angle}}k_{\theta} \frac{(\cos \theta - \cos \theta_0)^2}{\sin^2 \theta}

where :math:`k_{\theta}` is the force constant, :math:`\theta_0` is the equilibrium angle, and :math:`\theta` is the angle between the three sites.

To reproduce the intrinsic helical nature of RNA strands, the dihedral potential is defined as,

.. math::

   U_{\textrm{Dihedral}} = -\sum_{\textrm{dih}} k_{\phi} \exp{\left[-\frac{1}{2}w (\phi - \phi_0)^2\right]}

where :math:`\phi_0` is set to the dihedral angle of ideal A-type helix structure, and :math:`k_{\phi}` and :math:`w` are parameters to adjust the helical propensity.

Excluced-volume energy term
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. math::

   U_{\textrm{EV}} = \sum_{i,j>i+1}{\Theta\left(\sigma-r_{ij}\right) \varepsilon \left[ \left(\frac{\sigma}{r_{ij}} \right)^{12} - 2 \left(\frac{\sigma}{r_{ij}} \right)^{6} + 1\right]}.

Base-pair energy term
^^^^^^^^^^^^^^^^^^^^^

.. math::
   :label: eq_U_BP

    U_{\textrm{BP}} = \sum U_{\textrm{BP}}^0 \exp{\left[ U_{\textrm{BP,bond}} + U_{\textrm{BP,angle}} + U_{\textrm{BP,dihedral}} \right]}

where :math:`U_{\textrm{BP}}^0` is base pair stability and the three local geometry dependent terms are hamonic potentials defined as follows.

.. math::
   :label: eq_U_BP_terms

    \begin{split}
        U_{\textrm{BP,bond}} =& -k_r \left( r_{ij} - r_{0} \right) ^ 2, \\
        U_{\textrm{BP,angle}} =&-k_{\theta_1} \left( \theta_{i,j,j-1} - \theta_{1} \right) ^ 2
        -k_{\theta_2} \left( \theta_{i-1,i,j} - \theta_{2} \right) ^ 2 \\
        & -k_{\theta_3} \left( \theta_{i,j,j+1} - \theta_{3} \right) ^ 2
        -k_{\theta_4} \left( \theta_{i+1,i,j} - \theta_{4} \right) ^ 2, \\
        U_{\textrm{BP,dihedral}} =&
        -k_{\phi_1} \left[ 1 + \cos \left( \phi_{i-1,i,j,j-1} + \phi_{1} \right) \right] \\
        & -k_{\phi_2} \left[ 1 + \cos \left( \phi_{i+1,i,j,j+1} + \phi_{2} \right) \right].
    \end{split}

:math:`U_{\textrm{BP}}^0` 

Electrostatic energy term
^^^^^^^^^^^^^^^^^^^^^^^^^

The electrostatic potential was introduced into the previous SIS model in \citep{Maity23PNAS} to investigate salt-concentration dependent folding of trinucleotide RNA sequences. We use the same Debye-H\"uckel approximation with consideration of Oosawa-Manning ion condensation,

.. math::
   :label: eq_U_EL

    U_{\textrm{EL}} = \frac{Q_p^{2}e^{2}}{4\pi\varepsilon_{0}\varepsilon}\sum_{i<j}\frac{1}{r_{ij}}\exp\left(-\frac{r_{ij}}{\lambda_{D}}\right)

where :math:`\varepsilon_0` is the electric constant, :math:`e` is the elementary charge. 

Because we treat the water implicitly, :math:`\varepsilon` is the dielectric constant of water and we use a formula determined experimentally, which depends on temperature,

.. math::
   :label: eq_diele

   \varepsilon = 87.740-0.40008T_c+9.398\times10^{-4}T_c^{2}-1.410\times10^{-6}T_c^{3}

where :math:`T_c` is the temperature in degree Celsius [Malmberg1956]_.

The Debye length :math:`\lambda_D` is calculated as :math:`\lambda_{D}^{-2}=\frac{8\pi I}{\varepsilon(T)k_{B}T}` where :math:`I` is the ionic strength and :math:`k_B` is the Boltzmann constant. 
Finally, the magnitude of phosphate charge :math:`Q_p` is derived from the polyelectrolyte theory by Oosawa and Manning [Manning1969]_,

.. math::
   :label: eq_Qp

   Q_p=\frac{b\varepsilon(T)k_{B}T}{e^{2}}

where :math:`b` is the length per unit charge. For RNA molecules, the optimized value of :math:`b` was determined to be 0.44 nm in previous studies to match known thermodynamics of small RNAs [Denesyuk2013]_ [Nguyen2022]_.

Simulation methods
------------------


----

.. [Nguyen2022]
   H.T. Nguyen, N. Hori, D. Thirumalai (2022) Condensates in RNA Repeat
   Sequences are Heterogeneously Organized and Exhibit Reptation Dynamics.
   *Nat. Chem.* 14: 775–785.
   `10.1038/s41557-022-00934-z <https://doi.org/10.1038/s41557-022-00934-z>`__

.. [Bulacu2013]
   M. Bulacu, N. Goga, W. Zhao, G. Rossi, L. Monticelli, X. Periole, D.P. Tieleman, S.J. Marrink (2013) Improved Angle Potentials for Coarse-Grained Molecular Dynamics Simulations. *J. Chem. Theory Comput.* 9: 3282–3292
   `10.1021/ct400219n <https://doi.org/10.1021/ct400219n>`_

.. [Denesyuk2013]
   N.A. Denesyuk, and D. Thirumalai (2013) Coarse-Grained Model for Predicting RNA Folding Thermodynamics. *J Phys Chem B* 117: 4901–4911.
   `10.1021/jp401087x <https://doi.org/10.1021/jp401087x>`_

.. [Malmberg1956]
   C.G. Malmberg, and A.A. Maryott (1956) Dielectric constant of water from 0 to 100 C. *J Res Natl Inst Stan* 56: 1.
   `10.6028/jres.056.001 <https://doi.org/10.6028/jres.056.001>`_

.. [Manning1969]
    G.S. Manning (1969) Limiting Laws and Counterion Condensation in Polyelectrolyte Solutions I. Colligative Properties. *J Chem Phys* 51: 924–933
   `10.1063/1.1672157 <https://doi.org/10.1063/1.1672157>`_
  