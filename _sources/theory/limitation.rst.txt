Limitations
-----------

The :program:`HarP` procedure relies on several approximations that are listed below.

----

Physical Constraints
^^^^^^^^^^^^^^^^^^^^

        #. The driving laser field is assumed to be **linearly polarised**.
        #. The interaction is restricted to **cylindrical geometry**.
        #. The **polarisation response** is assumed to be induced solely by the **local fundamental field**. With current laser technology, the harmonic energy is several orders of magnitude lower than that of the fundamental field. Additionally, for a long laser pulse with respect to its oscillation period, only the local field intensity should be responsible for the non-linear polarisation.
        #. The **ionic susceptibility is neglected** with respect to the electronic and atomic susceptibility.
        #. The propagation equations are solved within the **paraxial approximation**.
        #. The temporal envelope is assumed to vary slowly compared to the laser oscillation period. Consequently, the **laser duration** :math:`\mathbf{T_{fwhm}}` must be significantly **longer** than the reference **optical period** :math:`\mathbf{T_r = 1/\omega}`.
        #. Although users can provide their own tables for the microscopic response, the tables currently distributed with :program:`HarP` are obtained by solving the 3D time-dependent Schrödinger equation within the **single active electron approximation**, following the methodology of [KulanderShore1990]_.
        #. The target **density** must be much **lower than the reference density** :math:`\mathbf{\rho_r}`. This condition ensures laser propagation through the target (no reflection), and guarantees the validity of the refraction index perturbation approximation :math:`\chi \ll 1`.
        #. Within the single active electron approximation, **only the first ionisation** is considered significant.
        #. The atomic form factors are provided by the user. While their values are well established for high-energy photons [Henke1993]_, [NIST]_, significant discrepancies may exist for low-order harmonics depending on the approximations used in the calculation. Consequently, variations of the absorption cross-section may give rise to differences in harmonics properties.

----

Numerical Constraints
^^^^^^^^^^^^^^^^^^^^^

        #. The dipole polarisabilities are interpolated from **precomputed tables**. These tables should provide high-resolution dipole moments as a function of laser intensity and harmonic order for a given set of wavelengths and atomic species.
        #. Only **two tables** are currently available: **argon and neon at** :math:`\mathbf{810}` **nm**. Extension to other wavelengths associated with conventional laser gain media are in progress.
        #. The spatial and temporal **grids are uniformly discretised**. While most harmonics are generated near :math:`r=0` in most cases, an adaptive grid could be appropriate to reduce the radial discretisation resolution. Additionally, for far-field projection using a Hankel transform, a non-uniform grid based on the zeros of the Bessel function would provide better resolution.
        #. The computational complexity is :math:`\mathbf{O((N_q+1) \times N_t \times N_z \times N_r)}`. The estimated execution time for a full simulation is :math:`6 \times (N_q+1) \times (N_r/512) \times (N_z/256) \times (N_t/256)` seconds, measured on an Apple MacBook Air M5 (16 GB RAM).
