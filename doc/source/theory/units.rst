Units
-----

The wave equation is solved using dimensionless variables. Physical quantities are normalised following standard Particle-in-Cell (PIC) code conventions [Derouillat2018]_.


----

Notation
^^^^^^^^

Standard notation is adopted for physical constants: :math:`c` for the speed of light, :math:`m_e` for the electron mass, 
:math:`e` for the elementary charge, :math:`\varepsilon_0` for the vacuum permittivity and :math:`\hbar` for the reduced Planck constant. 
The laser wavelength is denoted as :math:`\lambda`, its angular frequency as :math:`\omega=2\pi c/\lambda` and its vacuum wavevector :math:`k=\omega/c`.

----

Reference quantities
^^^^^^^^^^^^^^^^^^^^

To simplify the expression of the wave equation, all quantities used in  :program:`HarP` are expressed in normalised units based on the laser angular frequency :math:`\omega`. The corresponding reference units are:

        * Velocity :math:`V_r = c`
        * Time :math:`T_r = 1/\omega`
        * Length :math:`L_r = c/\omega` 
        * Energy :math:`\mathcal{E}_r = \hbar \omega`
        * Electric field :math:`E_r = m_e c \omega / e`
        * Density :math:`\rho_r = \varepsilon_0 m_e \omega^2 /e^2`
        * Polarisability :math:`\alpha_r = e^2/(m_e \omega^2)`
        * Second hyperpolarisability :math:`\gamma_r = e^4/(m_e^3 c^2 \omega^4)`
        * Dipole :math:`D_r = ec/\omega`

Note that, although the reference time is defined by :math:`T_r`, the temporal grid does not necessarily need to resolve this 
timescale. Only the envelope of the field needs to be resolved.
Similarly, the spatial grid resolution is determined by the laser waist and not by the laser wavelength.

.. warning::
     All quantities defined in the input file must be provided in normalised units unless other units are explicitly requested.

Output quantities can be automatically converted and saved in SI units upon request, see :doc:`namelist <../practice/namelist>`.

----

.. _newsystemofeq:

Dimensionless equations
^^^^^^^^^^^^^^^^^^^^^^^^

In this system of units, the propagation equations introduced in :doc:`HarP Equations <harp_equations>` section take the following dimensionless form:

.. math::
     \begin{gathered}
          \partial_r^2 E_q + \frac{1}{r}\partial_r E_q+ 2ik_{q,0}(z)\partial_z E_q +2k_{q,0} [\delta k_q(r,z)+ \Delta k_{q,0}(z)]E_q= - 4 q^2 \rho_a(r,z)d_q(r,z)e^{\mathbf{i}q\phi_1(r,z)}, \\
          \,k_{1,0}(z)= 1+\rho_g(z)\frac{\alpha_0}{2}, \quad \delta k_{1}(r,z)= -\frac{\rho_e(r,z)}{2}(1+\alpha_0)+\rho_a(r,z)\frac{3 }{8}\gamma E_1^2(r,z), \\
          \,k_{q,0}(z)= q-\frac{\rho_g(z)}{2q}[f_1(q)-\mathbf{i}f_2(q)], \quad \delta k_{q}(r,z)= -\frac{\rho_e(r,z)}{2q}\left[1-f_1(q)+\mathbf{i}f_2(q)\right].
     \end{gathered}

The partial derivatives are normalised to the reference length :math:`L_r` or time :math:`T_r`. The envelope profile :math:`E_q` is expressed in units of :math:`E_r`, the dipole :math:`d_q` in units of :math:`D_r`, and the polarisation and second hyperpolarisability are normalised to :math:`\alpha_r` and :math:`\gamma_r`, respectively. 
Likewise, the wavevectors :math:`k_{q,0}` and their perturbations :math:`\delta k_q` are normalised to :math:`L_r`. 
Consequently, the vacuum wavevector of the q-th harmonic simply equals :math:`q`.
The system is closed by the equations governing electron and atomic densities, both expressed in units of :math:`\rho_r`.

This set of equations is solved using the numerical procedure described in :doc:`HarP Algorithm <harp_algorithm>` section.
