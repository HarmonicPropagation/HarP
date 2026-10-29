HarP Equations
--------------

Several models have been developed in order to solve the macroscopic propagation of the harmonics. 
One of the earliest theoretical models capable of reproducing the main trends of harmonic generation and their propagation was introduced in [LHuillier1992]_.

The :program:`HarP` code aims to solve the equations introduced in this pioneering work that are re-derived in the current section.

In the following, we consider a **cylindrically symmetric** geometry using radial and longitudinal coordinates :math:`(r,z)`.
We defined the real electric field and nonlinear polarisation through their :math:`q` th Fourier components :math:`\mathcal{E}_q` and :math:`\mathcal{P}_q`:

.. math::
    \mathcal{E}(r,z,t)&= \frac{1}{2} \sum_q \mathcal{E}_q(r,z)e^{-\mathbf{i}q \omega t} +c.c. \\
    \mathcal{P}(r,z,t)&= \frac{1}{2} \sum_q \mathcal{P}_q(z,t)e^{-\mathbf{i}q \omega t} +c.c.

----

Wave Equation
^^^^^^^^^^^^^

The derivation begins from the wave equation in Fourier space for the electric field of the harmonic of order :math:`q`:

.. math::
    \nabla^2 \mathcal{E}_q + k_q^2 \mathcal{E}_q=-\frac{1}{\epsilon_0} \left(\frac{q\omega}{c}\right)^2 \mathcal{P}_q

The **polarisation** is assumed to be **induced solely by the fundamental field** (:math:`q=1`).
The wave vector of the :math:`q` th harmonic is defined as :math:`k_q=nq\omega/c` where :math:`n` is the refractive index of the medium evaluated at energy :math:`q \hbar \omega` and position :math:`(r,z)`. 
**Neglecting the ionic susceptibility**, the index of refraction is determined by the atomic and electronic, denoted :math:`\chi_a` and :math:`\chi_e`, respectively:

.. math::
    n(q\omega,r,z) &= \sqrt{1+\chi_a(q\omega,r,z)+\chi_e(q\omega,r,z)}, \\
    \chi_a(q \omega,r,z) &= -\rho_a(r,z)\frac{ e^2}{\epsilon_0 m_e q^2 \omega^2}\left[f_1(q\omega) -\mathbf{i} f_2(q\omega)\right], \\
    \chi_e(q \omega,r,z) &= -\rho_e(r,z)\frac{ e^2}{\epsilon_0 m_e q^2 \omega^2},

Here :math:`f_1`, :math:`f_2`, represent the atomic form factors. The atomic form factor :math:`f_2` is responsible for the absorption and should be understood as a positive quantity in agreement with [Henke1993]_ table.
Here, :math:`\rho_a` and :math:`\rho_e` correspond to the atomic and electronic densities.

----

Wavevector
^^^^^^^^^^

The wavevector is now expressed using a perturbative expansion: :math:`k_q(r,z)=k_{q,0}(z)+\delta k_q(r,z)`. 
The zeroth-order term :math:`k_{q,0}` is assumed to be independent of the transverse coordinate,
while the first-order perturbation :math:`\delta k_q` satisfies :math:`k_{q,0}(z) \gg \delta k_q(r,z)` everywhere in space. 
We assume a **low-density target** such that :math:`\chi_{a}, \chi_{e} \ll 1` and consider only the first ionisation relevant to the interaction such that the ionic density equals the electronic density. Thus, the initial target density is :math:`\rho_g(z)=\rho_e(r,z)+\rho_a(r,z)` (assumed radially uniform at this stage). We have :math:`\forall q>1`:

.. math::
    :label: EqWavevectors

    \begin{gathered}
        \frac{c}{q\omega}\,k_{q,0}(z)= 1-\rho_g(z)\frac{e^2}{2\epsilon_0m_eq^2\omega^2}[f_1(q \omega)-\mathbf{i}f_2(q\omega)] \\ 
        \frac{c}{q\omega} \delta k_{q}(r,z)= -\rho_e(r,z)\frac{e^2}{2\epsilon_0m_eq^2\omega^2}\left[1-f_1(q\omega)+\mathbf{i}f_2(q\omega)\right].
    \end{gathered}

This decomposition follows the approach introduced in Refs. [LHuillier1992]_, [Antoine1996]_. 
As discussed in Ref. [Antoine1996]_, an alternative formulation can be adopted by setting :math:`k_{q,0}=q\omega/c` 
and including the spatial dependencies into the perturbation term: :math:`\delta k_q(r,z)=\chi_a(q\omega,r,z)/2+\chi_e(q\omega,r,z)/2`. 
This latter approach is implemented for initial target profiles exhibiting transverse variations (see :doc:`namelist <../practice/namelist>`). 
For radially uniform initial target profiles, the system defined in Eqs. :eq:`EqWavevectors` provides higher accuracy.

The wavevector of the fundamental field is defined directly through the polarisability :math:`\alpha_0` rather than the form factors.
Furthermore, at high laser intensities (of the order of or exceeding :math:`10^{14}` W/cm\ :sup:`2`), the Kerr effect [Kerr1875]_ must be accounted for in the pulse propagation. 
Adding the Kerr nonlinearity to the refractive index [Boyd2008]_ for :math:`q=1`, yields the wavevector components :math:`k_{1,0}(z), \delta k_1(r,z)`:

.. math::
    
    \begin{gathered}
        \frac{c}{\omega} \,k_{1,0}(z)= 1+\rho_g(z)\frac{\alpha_0}{2\epsilon_0}, \\
        \frac{c}{\omega} \delta k_{1}(r,z)= -\frac{\rho_e(r,z)}{2}\left(\frac{\alpha_0}{\epsilon_0}+\frac{e^2}{\epsilon_0m_eq^2\omega^2}\right)+\rho_a(r,z)\frac{3\gamma }{8\epsilon_0}E_1^2(r,z).
    \end{gathered}

where :math:`\gamma` represents the second hyperpolarisability. The definition of the second hyperpolarisability may differ by a factor of 6 across different communities. 
In :program:`HarP`, the definition of :math:`\gamma` is 6 times smaller than the value reported in Ref. [Soldán2001]_.

----

Envelope Equation
^^^^^^^^^^^^^^^^^

We define the modified envelope function of the electric field :math:`\hat{E}_q` and of the polarisability :math:`P_q` through the following relations:

.. math::

    \mathcal{E}_q&=\hat{E}_q e^{i\int^zk_{q,0}(z') dz'}, \\ 
    \mathcal{P}_q&= P_q e^{i\int^z qk_{1,0}(z') dz'}.

Applying the **paraxial approximation** :math:`\partial_z^2 \hat{E}_q \ll k_{q,0}\partial_z \hat{E}_q` and neglecting the term in :math:`\hat{E}_q \partial_z k_{q,0}`,
the wave equation simplifies to an equation of the envelope field:

.. math::
    \nabla_\perp^2 \hat{E}_q + 2ik_{q,0}(z)\partial_z \hat{E}_q +2k_{q,0} \delta k_q(r,z) \hat{E}_q=-\frac{1}{\epsilon_0} \left(\frac{q\omega}{c}\right)^2 P_q e^{-i\int^z \Delta k_{q,0}(z')dz'},

where :math:`\nabla_\perp^2 = \partial_r^2  + 1/r \partial_r` and :math:`\Delta k_{q,0}(z)=k_{q,0}(z)-qk_{1,0}(z)`. Introducing the change of variable:


.. math::
    \hat{E}_q=E_q e^{-i\int^z \Delta k_{q,0}(z')dz},

we obtain the governing equation:

.. math::
    :label: EqEnvelopeQ

    \nabla_\perp^2 E_q + 2ik_{q,0}(z)\partial_z E_q +2k_{q,0} [\delta k_q(r,z)+ \Delta k_{q,0}(z)]E_q=-\frac{1}{\epsilon_0} \left(\frac{q\omega}{c}\right)^2 P_q.

Since :math:`k_{1,0}` is purely real for optical lasers, the energy calculated from :math:`E_q` is identical to that calculated from the full electric field :math:`\mathcal{E}_q`.

It is essential to distinguish between the optical period of the laser oscillations and its temporal envelope. 
Under the slowly varying envelope approximation (SVEA), these two time-scales can be treated independently. 
The rapid oscillations are treated in the frequency domain by decomposing the field into its harmonic components, while the corresponding field envelopes remain functions of time. 
Consequently, although the propagation equations are formulated in the frequency domain for each harmonic, they are solved at successive times discretising the slowly varying envelopes. 
The same separation of time scales applies to other slowly varying quantities, such as plasma density and wave vectors.


----

Dipole Moment
^^^^^^^^^^^^^

The non-linear polarisation term :math:`P_q` is defined through the atomic dipole moment :math:`d_q(r,z)` and the phase of the fundamental field :math:`\phi_1(r,z)`:

.. math::
    :label: EqPq

    P_q(r,z) = 2\rho_a(r,z)d_q(r,z)e^{\mathbf{i}q\phi_1(r,z)}.

This term, directly responsible for high-harmonic generation, contains the microscopic response of the atom interacting with the intense incident laser field. 
It can be modelled using several theoretical frameworks. Numerous codes have already been developed for calculating the non-linear polarisability within the Strong Field Approximation (SFA) [Lewenstein1994]_. 
A non-exhaustive list of such approaches includes [Priori2000]_, [Shon2000]_, [Gaarde2008]_, [Tosa2009]_, [Tosa2016]_, and [Schröder2025]_.
While computationally efficient, these codes are limited by the SFA and can suffer from inaccuracies compared to dipole moments calculated from the Time-Dependent Schrödinger Equation (TDSE).

An alternative approach was introduced by [Vábek2025]_. In that procedure, the wave equation for the fundamental field is solved in the time domain.
The microscopic response is then calculated at each time step by solving the 1D-TDSE under the Single Active Electron (SAE) approximation. 
The harmonics are finally propagated using a Hankel transform. Although less restrictive, this approach is computationally demanding.

In the :program:`HarP` code, the dipole moment is provided from precomputed external tables.
These tables map the dipole moment (absolute value and phase) as a function of the laser intensity and harmonic order (see :doc:`namelist <../practice/namelist>` for the appropriate format). 

Precomputed tables allow :program:`HarP` to incorporate sophisticated microscopic response models without increasing the computational cost. 
However, these tables are linked to a specific laser wavelength and atomic species.
Additionally, this approach relies on the assumption that the dipole moment is a function of the local fundamental field only. This being justified by the slowly varying envelope approximation.

In the current version of :program:`HarP`, two tables, computed from **3D-TDSE** within **SAE** approximation, are provided for an
incident laser wavelength of :math:`810` nm interacting with Neon and Argon. These tables were generated following the methodology introduced in [KulanderShore1990]_.
A driving field featuring a five-cycle ramp-up, a five-cycle constant-amplitude plateau, and a five-cycle ramp-down is applied to the atom.
The interaction is simulated by solving the 3D-TDSE as modelled in [Carlstrom2022]_. 
Once transient dynamics have subsided, the dipole moment is Fourier-transformed over the plateau region. Repeating this simulation across a range of peak laser intensities yields a map of the dipole response as a function of laser intensity and harmonic order for a given atomic species and driving wavelength.

Since the dipole tables are associated with a specific wavelength, the tabulated approach requires a pulse duration that is long relative to the optical period (:math:`\mathbf{T_{fwhm}\gg 1/\omega}`). 
Otherwise, spectral broadening of the fundamental field causes multiple photon energies to contribute to harmonic generation.

.. warning::
    Because of the spin degeneracy in the SAE approximation, the tabulated dipole moment and ionisation rate are a factor of 2 smaller than the ones defined in Eqs. :eq:`EqPq` and :eq:`EqDensity`. This factor is taken into account when loading the tables through the ``statistical_factor`` arguments.


----

Density
^^^^^^^

As the intense laser pulse propagates through the target medium, it ionises the gas, altering the local atomic and electronic densities. This evolution is governed by the population equations:


.. math::
    :label: EqDensity
    
    \begin{aligned}
        \rho_a(r,z,t) &= \rho_g(r,z)e^{-\int^t w(E_1(r,z,t'))dt'}, \\
        \rho_e(r,z,t) &= \rho_g(r,z)-\rho_a(r,z,t) 
    \end{aligned}

where :math:`w` denotes the ionisation rate. Because the fundamental photon energy is significantly lower than the harmonic energies, the ionisation rate is calculated solely from the local fundamental field amplitude.

The ionisation rate is evaluated by interpolating precomputed tables, mapping :math:`w` to the laser intensity (see :doc:`namelist <../practice/namelist>` for the appropriate format).
These tables can be generated using standard cycle-averaged ionisation rate models such as PPT [Perelomov1966]_ or ADK [Ammosov1986]_.
For maximum accuracy, the ionisation rates can also be extracted directly from 3D-TDSE SAE calculations, also used for the dipole calculation. 
The current version of :program:`HarP` provides 3D-TDSE SAE rate tables for an :math:`810` nm laser propagating in Argon and Neon.
