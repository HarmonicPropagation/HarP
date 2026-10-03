Far Field Propagation
---------------------

Since the HarP algorithm is limited by the longitudinal box size, to propagate the field far from the exit of the gas, it is convenient to solve the field propagation in vacuum using the Hankel transform. 
In this section, we derive the mathematical and numerical framework for this purpose.

----

Hankel transformation
^^^^^^^^^^^^^^^^^^^^^

For the sake of completeness, we remind the Hankel transform of order :math:`\nu` of an arbitrary function :math:`f(r)` [HankelTransform]_:

.. math::
    F_\nu(k) = \int_0^\infty f(r)J_\nu (k r) r dr

and its inverse transformation:

.. math::
    f(r) = \int_0^\infty F_\nu(k)J_\nu (k r) k dk,

where :math:`J_\nu` denotes the Bessel function of the first kind of order :math:`\nu\geq -1/2`. By construction, the Bessel function :math:`J_\nu` satisfies the differential equation:

.. math::
    x^2 J_\nu''(x) + J_\nu'(x) + x^2 J_\nu(x) = \nu^2 J_\nu(x),

where the prime notation corresponds to the successive derivatives of :math:`J_\nu(x)` with respect to :math:`x`.

----

Envelope equation in vacuum
^^^^^^^^^^^^^^^^^^^^^^^^^^^

In vacuum the **paraxial Helmholtz equation** for the field envelope of harmonic :math:`q` reads:

.. math::
    :label: EqParaxialHelmholtz

    r^2 \partial_r^2 E_q + r \partial_r E_q + 2i \frac{q\omega}{c}r^2 \partial_z E_q = 0. 

    
This could be directly derived from Eq. :eq:`EqEnvelopeQ` noting that there is no variation of the background densities (in vacuum :math:`\Delta k_q = 0` and :math:`\delta k_q = 0`).

Let :math:`\tilde{E}_{q}` be the Hankel transform of order :math:`0` of the field envelope :math:`E_{q}`:

.. math:: 
    E_{q}(r,z) = \int_0^\infty \tilde{E}_{q}(k,z) J_0 (k r) k dk


and injecting in the paraxial Helmholtz equation Eq. :eq:`EqParaxialHelmholtz`, we obtain the differential equation for the Hankel transform:

.. math:: 
    -k^2 \tilde{E}_{q}(k,z) + 2i \frac{q\omega}{c} \partial_z \tilde{E}_{q}(k,z) = 0,

which admits the solution:

.. math:: 
    \tilde{E}_{q}(k,z) = \tilde{E}_{q}(k,z=z_0) \exp{\left[-i \frac{c k^2}{q\omega} (z-z_0)\right]}.

These properties are used to propagate the field far from the exit of the target as described below.

----

Numerical procedure
^^^^^^^^^^^^^^^^^^^

The field envelope at the end of the simulation box :math:`E_q(r,z=Lz)` is calculated using the :doc:`HarP algorithm <../theory/harp_algorithm>`. 
The numerical solution is used as an input for the propagation at any given longitudinal position :math:`z \geq Lz`. 
To compute the field far from the exit box, we apply the following three distinct steps:

    #. First, the Hankel transform :math:`\tilde{E}_{q}` of the field envelope at the rear of the target is calculated according to:
        .. math:: 
            \tilde{E}_{q}(k,z=Lz) = \int_0^\infty E_{q}(r,z=Lz) J_0 (k r) r dr
    #. Second, it is propagated to an arbitrary position using:
        .. math:: 
            \tilde{E}_{q}(k,z) = \tilde{E}_{q}(k,z=Lz) \exp{\left[-i \frac{c k^2}{q\omega} (z-Lz)\right]}
    #. Finally, the field envelope is reconstructed from the inverse Hankel transform:
        .. math:: 
            E_{q}(r,z) = \int_0^\infty \tilde{E}_{q}(k,z) J_0 (k r) k dk


At each time slice, when the field reaches the simulation box limit, these three steps are applied to propagate the field. 
The vector :math:`k` is uniformly distributed from :math:`0` to :math:`\min(\pi / \Delta r, \alpha q\omega/c )`. 
The limit corresponds to the minimum between the Nyquist frequency and the paraxial approximation. It ensures that at least two points discretise the wavelength through the condition :math:`k<\pi / \Delta r` (Nyquist limit). 
It also ensures that the angle between the wave and the optical axis remains small (paraxial approximation :math:`k_\perp \ll q\omega/c`) through the condition :math:`k_\perp \ll \alpha q\omega/c` where :math:`\alpha` is an arbitrary number below unity. 
This discretisation is used to compute the Hankel transform integrals with a Riemann sum. 

.. note::
    As the transverse dimension is uniformly discretised in the HarP loop, the wavevector :math:`k` is also uniformly discretised.
    However, for the numerical evaluation of the Hankel transforms, a discretisation based on the zeros of the Bessel function provides a more appropriate choice to ensure the grid resolves the oscillation of :math:`J_0`.
    We therefore plan to implement a non-uniform, Bessel-zero-based discretisation in future work.


This calculation is not part of the HarP algorithm but a branch inside the Diagnostics ``DiagFarField`` (practical application described in the :doc:`namelist <../practice/namelist>`). 
It allows the efficient calculation of the field far from the exit of the gas without increasing the simulation box length.



