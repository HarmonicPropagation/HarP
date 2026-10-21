Namelist
--------

The input for  :program:`HarP`  is a Python script that defines the simulation geometry through a set of classes. Their dependencies and relationships are described below.

----

Main
^^^^

.. py:class:: Main

    The main class stores all the general variables related to the numerical calculation::

        Main(
            wavelength = 0.81*1e-6,
            solver = "CrankNicolson",
            boundary_condition = "absorbing",
            anne_eq = True,
            predictor_corrector = False,
            dr = None,
            dz = None,
            nr = 512,
            nz = 256,
            Lr = 10.,
            Lz = 100.,
            dt = None,
            nt = 256,
            simulation_time = 100.,
            verbose = True,
            print_every = 10,
            approx_check = False
        )

    .. py:attribute:: wavelength

        :type: Float
        :default: ``810*1E-9``

        Reference wavelength of the incident laser. The value is provided in SI (meters). It is used to link SI units with code units.

        .. note::

            In case harmonics generation is considered, make sure harmonics tables are related to this wavelength.

    .. py:attribute:: solver

        :type: String
        :default: ``"CrankNicolson"``

        Scheme used for the solution of the propagation equation. 

        .. warning::

            Only the Crank-Nicolson scheme is available in the current HarP version.
    
    .. py:attribute:: boundary_condition

        :type: String
        :default: ``"absorbing"``

        Boundary condition at the outermost radial point.
        
            * Set to ``"absorbing"`` for the Sommerfeld absorbing condition :math:`\partial_r E_q(r_{max},z)=\mathbf{i}k_q^0 E_q(r_{max},z)`.
            * Set to ``"reflective"`` for the Dirichlet boundary condition :math:`E_q(r_{max},z)=0`.
        
    .. py:attribute:: anne_eq

        :type: Boolean
        :default: ``True``

        Definition of the wavevector from [LHuillier1992]_ (Eqs. :eq:`EqWavevectors`). If false, all the spatial dependencies of the density are introduced inside the perturbation :math:`\delta k_q`.

        .. warning::

            Set ``anne_eq=False`` if the initial target density depends on the radial coordinate.

    .. py:attribute:: predictor_corrector

        :type: Boolean
        :default: ``True``

        Uses the Predictor Corrector method for the calculation of the perturbation :math:`\delta k_1 \propto \lvert E_1 \rvert^2`.

    .. py:attribute:: dr, dz, dt

        :type: Float
        :default: ``None``

        Radial, longitudinal and temporal spacing of the grid (code units).

    .. py:attribute:: nr, nz, nt

        :type: Integer
        :default: ``None``

        Number of radial, longitudinal and temporal grid points.

    .. py:attribute:: Lr, Lz

        :type: Float
        :default: ``None``

        Radial and longitudinal length of the box (code units).

    .. py:attribute:: simulation_time

        :type: Float
        :default: ``None``

        Temporal duration of the interaction (code units).

    .. note::

        Either the number of grid points or the grid spacing must be specified.
        The computational domain is uniformly discretised with steps ``dr``, ``dz``, and ``dt`` over the intervals ``[0, Lr]``, ``[0, Lz]``, and ``[0, simulation_time]``, respectively.

    .. py:attribute:: print_every

        :type: Integer
        :default: ``nt-1``

        Number of time steps between successive prints of the elapsed simulation time. By default only the first and the last iteration are printed.

    .. py:attribute:: verbose

        :type: Boolean
        :default: ``True``

        Enable or disable console output during the simulation.
    
    .. py:attribute:: approx_check

        :type: Boolean
        :default: ``False``

        Enable or disable the test of the different approximations during the calculation (see :doc:`limitation <../theory/limitation>` for the approximations).

        
                    * The low density approximation is tested at the beginning of the calculation. In the case of a cell with a density larger than :math:`n_r/5`, a warning is printed.
                    * The slowly varying envelope approximation in time is tested at the begining of the calculation comparing :math:`\lvert \frac{dE_1}{dt}\rvert` with :math:`\omega E_1`.  In the case of a ratio larger than :math:`0.2`, a warning is printed.
                    * The paraxial approximation is tested during the HarP loop comparing :math:`\lvert \frac{dE_1}{dz}\rvert` with :math:`k E_1`. In the case of a ratio larger than :math:`0.2`, a warning is printed.

----

Pulse
^^^^^

.. py:class:: Pulse

    The pulse is defined by three components: its amplitude, radial profile, and temporal envelope. At each time step, the electric field envelope at :math:`z=0` is reconstructed as the product of these three quantities::

        Pulse(
            amplitude                   = E0,
            radial_profile              = eprofile,
            time_envelope               = tenv
        )

    .. py:attribute:: amplitude

        :type: Float
        :default: ``1.``

        Amplitude factor in front of the time envelope profile (code units).

    .. py:attribute:: radial_profile

        :type: Python function of argument (:math:`r`)
        :default: ``None``

        Radial profile of the electric field envelope of the laser at :math:`z=0`.

        :Example:

            For a Gaussian beam, where ``waist`` is the beam waist, ``z_R`` the Rayleigh length, and ``z_focus`` the focal position:

            .. code-block:: python

                def radial_profile(r):   # radial profile of the pulse at z=0
                    wz     = waist * np.sqrt(1+ (zfocus/z_R)**2)
                    R      = zfocus*(1.+(z_R/zfocus)**2)
                    PhiG   = np.arctan2(zfocus,z_R)
                    return waist/wz * np.exp(-(r/wz)**2)*np.exp(-1j*(r**2/2./R-PhiG))

    .. py:attribute:: time_envelope

        :type: Python function of argument (:math:`t`)
        :default: ``None``

        Temporal profile of the electric field envelope of the laser.

        :Example:

            For a Gaussian beam, assuming ``tfwhm`` the temporal full width at half maximum in intensity and ``tfocus`` the temporal position of the pulse peak:

            .. code-block:: python

                def time_envelope(t):
                    return np.exp(-2*np.log(2)*(t-tfocus)**2/tfwhm**2)

----

Medium
^^^^^^

.. py:class:: Medium
    
    The medium properties, including the initial target profile, polarizability, second hyperpolarizability and form factors::

        Medium(
            gaz_profile      = gazprofile,
            polarizability_0 = 1.8247e-40,      # in SI from P. Soldan et al. 2001 PCCP
            gamma_kerr       = 7.2767118e-62/6.,  # in SI from P. Soldan et al. 2001 PCCP
            table_path       = "form_factor/f12_Ar_henkegov.txt",
            unit_ev          = 1.,         # photon energy in ev
            unit_f1          = 1.,        # polarizability as form factor f1
            unit_f2          = 1.,         # absorption cross section as form factor f2
            skiprows         = 1
        )

    .. py:attribute:: gaz_profile

        :type: Python function of arguments (:math:`r,z`)
        :default: ``None``

        Target initial density as a function of :math:`r` and :math:`z`.

        :Example:

            For a simple Heaviside target profile of density ``n0`` and length ``Lgaz`` :

            .. code-block:: python

                def gazprofile(r,z):
                    """if depend on r need to set anne_eq=False"""
                    return n0*np.heaviside(z-0.,1.)*np.heaviside(Lgaz-z,1.)

    .. py:attribute:: polarizability_0_au

        :type: Float
        :default: ``None``

        Value of the polarizability at optical wavelengths in atomic units. Some of them are given in [Soldán2001]_.

    .. py:attribute:: polarizability_0

        :type: Float
        :default: ``1.8247e-40``

        Value of the polarizability at optical wavelengths in SI. If ``polarizability_0_au`` is already defined, the value considered is the one in atomic units. (The default value is the Argon polarisability).

    .. py:attribute:: gamma_kerr_au

        :type: Float
        :default: ``None``

        Value of the second hyperpolarizability in atomic unit. Some of them are given in [Soldán2001]_.

        .. warning::

            Due to different definitions of the second hyperpolarizability, the values from [Soldán2001]_ have to be divided by 6 in the namelist.

    .. py:attribute:: gamma_kerr

        :type: Float
        :default: ``0.``

        Value of the second hyperpolarizability in SI. If ``gamma_kerr_au`` is already defined the value considered is the one in atomic units.

    .. py:attribute:: table_path

        :type: String
        :default: ``""``

        Path to the table containing the form factors :math:`f_1` and :math:`f_2` as a function of the photon energy. The table is organised into three columns: the first column contains the photon energy, the second contains :math:`f_1`, and the third contains :math:`f_2`. The units can be adjusted using the variables ``unit_ev``, ``unit_f1`` and ``unit_f2``.

    .. py:attribute:: unit_ev

        :type: Float
        :default: ``1.``

        Conversion factor applied to the photon's energy row to express values in eV (default unit of the Henke tables).

    .. py:attribute:: unit_f1

        :type: Float
        :default: ``1``

        Factor applied to the form factor :math:`f_1` row. For the Henke table, set the value to ``1``.

    .. py:attribute:: unit_f2

        :type: Float
        :default: ``1``

        Factor applied to the form factor :math:`f_2` row. For the Henke table, set the value to ``1``.

    .. py:attribute:: skiprows

        :type: Integer
        :default: ``3``

        Number of rows skipped before reading the tables.

----

Ionisation
^^^^^^^^^^

.. py:class:: Ionisation

    The ionisation class stores the table of ionisation rates::

        Ionisation(
            table_path          = "ionisation/ionisation_Ar_810_tdse.txt",
            skiprows            = 2,
            intensity_unit_Wcm2 = 1.,
            rate_unit_SI        = 2.   #add statistical factor of 2 on the rate due to SAE.
        )

    .. py:attribute:: table_path

        :type: String
        :default: ``""``

        Path to the table containing the ionisation rate :math:`w` as a function of the laser intensity. The table's first column corresponds to the laser intensity and the second column to the ionisation rate. The units can be adjusted using the variables defined below.

    .. py:attribute:: intensity_unit_Wcm2

        :type: Float
        :default: ``1``

        Factor applied to the intensity row to express values in W/cm\ :sup:`2`.

    .. py:attribute:: rate_unit_SI

        :type: Float
        :default: ``1``

        Factor applied to the rate row to express values in s :sup:`-1`.
    
    .. py:attribute:: statistical_factor

        :type: Float
        :default: ``2``

        Statistical factor applied to the rate values when loading the table. If the tables have been computed withing the SAE approximation set ``statistical_factor=2`` to account for the spin degeneracy.

    .. py:attribute:: skiprows

        :type: Integer
        :default: ``0``

        Number of rows skipped before reading the tables.

----

Harmonic
^^^^^^^^

.. py:class:: Harmonic

    The class for the harmonics propagation stores the tables of the dipole moment::

        Harmonic(
            order = [23,25,27,29],        #photon energy corresponding to order*hbar*omega
            table_path = tables_path+'tdse/dipoles_Ar_810.txt',
            intensity_unit_Wcm2 = 1e14,
            statistical_factor  = 2
        )

    .. py:attribute:: order

        :type: List or numpy array
        :default: ``None``

        Harmonic orders considered in the propagation. The order list should **not include 1**. It should only contain harmonics for which data are available in the table.

    .. py:attribute:: table_path

        :type: String
        :default: ``""``

        Path to the table containing the dipole moment values :math:`d_q` as a function of the laser intensity and harmonic order :math:`q`.
        The data block corresponding to the harmonic order :math:`q` should start with the following identifier ``harmonic      q``.
        The subsequent table is organised into three columns: the laser intensity in units of ``intensity_unit_Wcm2`` W/cm\ :sup:`2` (:math:`10^{14}` W/cm\ :sup:`2` for the provided tables), the square of the dipole absolute value :math:`\lvert d_q \rvert^2` in a.u., and the dipole phase :math:`arg(d_q)` (in radians), respectively.
        This structure is repeated for each harmonic order tabulated.

    .. py:attribute:: intensity_unit_Wcm2

        :type: Float
        :default: ``1e14``

        Factor applied to the intensity row to express values in W/cm\ :sup:`2`.
    
    .. py:attribute:: statistical_factor

        :type: Float
        :default: ``2``

        Statistical factor applied to the dipole value when loading the table. If the tables have been computed withing the SAE approximation set ``statistical_factor=2`` to account for the spin degeneracy.

----

Diagnostics
^^^^^^^^^^^

.. py:class:: Diagnostics

    To save data during and after the simulation, the :program:`HarP` input uses Diagnostics classes.
    Instead of directly handling three-dimensional arrays, these classes enable the continuous storage of one-dimensional arrays in **HDF5** files.
    The data are stored in a file named ``diag_path+"Diag"+deposited_quantity+add_suffix+".h5"``.
    For deposited quantities depending on the harmonic order, the first-level keys are organised as ``"order_XX"``. Each order key contains a set of datasets named ``"t_XX"``. Each dataset is a 2D array indexed ``[z, r]``, at the saved time :math:`t_{XX}` (``XX`` counts saved steps, i.e. ``it // every_t``).
    The spatial and temporal coordinates of the stored data are provided through the first-level keys ``"r"``, ``"z"``, and ``"t"``.

    The implementation is similar to and inspired by [SMILEI]_::

        Diagnostics(
            diag_path          = "",
            add_suffix         = "",
            deposited_quantity = "Field",
            order              = [1,23],
            every_r            = 1,
            every_z            = 1,
            every_t            = 32,
            save_every         = 32,
            SI                 = True
        )

    .. py:attribute:: diag_path

        :type: String
        :default: ``""``

        Path to the directory where the diagnostic is stored.

    .. py:attribute:: add_suffix

        :type: String
        :default: ``""``

        String added at the end of the filename.

    .. py:attribute:: deposited_quantity

        :type: String
        :default: ``""``

        The physical quantity is chosen in the following list:

            * ``"Field"``  modified electric field envelope :math:`E_q` ``complex128``
            * ``"RealField"`` electric field envelope :math:`\mathcal{E}_q` ``complex128``
            * ``"Phase"`` fundamental envelope phase :math:`\int^z_0 dz' k_1(z')` ``complex128``
            * ``"S"`` non-linear polarisability :math:`P_q` ``complex128``
            * ``"LocalRate"`` ionisation rate :math:`w(r,z,t)` ``float64``
            * ``"IntegratedRate"`` integrated ionisation rate  :math:`\int_0^t dt' w(r,z,t')` ``float64``
            * ``"ElectronDensity"`` electron density :math:`n_e(r,z,t)` ``float64``
            * ``"AtomicDensity"`` atomic density :math:`n_a(r,z,t)` ``float64``
            * ``"GazDensity"`` initial target density :math:`n_g(r,z)` ``float64``
            * ``"k"`` wavevector :math:`k_q(r,z,t)` ``complex128``
            * ``"dk"`` perturbation of the wavevector :math:`\delta k_q(r,z,t)` ``complex128``

   .. py:attribute:: order

        :type: List or numpy array
        :default: ``[]``

        List of integer associated with the harmonic orders. The integer :math:`1` corresponds to the diagnostic of the fundamental field. If :math:`1` is not provided by the user, the diagnostic will still save the data for the fundamental field.

    .. py:attribute:: every_r

        :type: Integer
        :default: ``1``

        Number of radial steps between two consecutive data saves. For example, if ``every_r=2``, only the deposited quantities at :math:`r=r_{2j}` are stored.

    .. py:attribute:: every_z

        :type: Integer
        :default: ``1``

        Number of longitudinal steps between two consecutive data saves. For example, if ``every_z=4``, only the deposited quantities at :math:`z=z_{4j}` are stored.

    .. py:attribute:: every_t

        :type: Integer
        :default: ``1``

        Number of times steps between two consecutive data saves. For example, if ``every_t=8``, only the deposited quantities at :math:`t=t_{8j}` are stored.

    .. py:attribute:: save_every

        :type: Integer
        :default: ``1``

        The data are saved at the final time step of the simulation. Additionally, they are periodically saved and reset whenever the number of collected time steps reaches ``save_every``.
        Opening, writing, and closing files introduces an additional overhead. This variable is therefore used to limit the saving frequency.
        Conversely, reducing the saving frequency may result in increased memory consumption and saturation of the available RAM.

    .. py:attribute:: SI

        :type: Boolean
        :default: ``True``

        If ``True``, the deposited quantities are saved in SI units. Otherwise, they are saved in code units.

----

DiagEnergy
^^^^^^^^^^

.. py:class:: DiagEnergy

    The energy contained inside the harmonics of order q can be saved using the ``DiagEnergy`` class.
    The data is stored in an **HDF5** file named ``diag_path+"DiagEnergy"+add_suffix+".h5"``.
    The first-level keys correspond to the harmonic order and are named ``"order_XX"``. Each order key contains the energy  of the corresponding harmonic as a function of the longitudinal position :math:`z`. In SI units, the harmonic energy reads:

    .. math::

        \pi \epsilon_0 c \int_0^T \int_0^{Lr} \lvert E_q(r,z,t) \rvert^2 r dr dt

    where T corresponds to the total simulation time. It is computed internally in code units using a Riemann sum.
    The longitudinal coordinate is available through the first-level key ``"z"``. The class is defined as follows::

        DiagEnergy(
            diag_path  = "",
            add_suffix = "",
            order   = [1,23],
            every_z = 1,
            SI      = True
        )

    .. py:attribute:: diag_path

        :type: String
        :default: ``""``

        Path to the directory where the diagnostic is stored.

    .. py:attribute:: add_suffix

        :type: String
        :default: ``""``

        String added at the end of the filename.

    .. py:attribute:: order

        :type: List or numpy array
        :default: ``[]``

        List of integer associated with the harmonic orders. The integer :math:`1` corresponds to the diagnostic of the fundamental field. If :math:`1` is not provided by the user, the diagnostic will still save the data for the fundamental field.

    .. py:attribute:: every_z

        :type: Integer
        :default: ``1``

        Number of longitudinal steps between two consecutive data saves. For example, if ``every_z=4``, only the deposited quantities at :math:`z=z_{4j}` are stored.

    .. py:attribute:: SI

        :type: Boolean
        :default: ``True``

        If ``True``, the energy is saved in Joule. Otherwise, it is saved in code units.

----

DiagFarField
^^^^^^^^^^^^

.. py:class:: DiagFarField

    To obtain the field far from the end of the simulation box, the DiagFarField class uses the Hankel transformation to propagate the field into vacuum at an arbitrary position :math:`z>Lz`. The theoretical framework is explained in the :doc:`Far field propagation <../theory/vacuum_propagation>` section.
    
    The data is stored in an **HDF5** file named ``diag_path+"DiagFarField"+add_suffix+".h5"``.
    The class parameters are defined as follows::

        DiagFarField(
            diag_path  = "",
            add_suffix = "",
            positions  = [],
            order   = [1,23],
            paraxial_max_factor = 0.5,
            every_r = 1,
            every_t = 1,
            SI      = True
        )

    .. py:attribute:: diag_path

        :type: String
        :default: ``""``

        Path to the directory where the diagnostic is stored.

    .. py:attribute:: add_suffix

        :type: String
        :default: ``""``

        String added at the end of the filename.
    
    .. py:attribute:: positions

        :type: List or numpy array
        :default: ``[]``

        List of longitudinal positions (in code units) at which the field is evaluated. The positions are defined with respect to the end of the simulation box, meaning that the origin of these positions is :math:`Lz`.

        :example:

            To get the field :math:`1` mm and :math:`1000` mm after the box, set

            .. code-block:: python

                positions=[1*mm,1000*mm],
            
            where mm :math:`=10^{-3}/L_r` is the conversion of millimeter distance in code units.


    .. py:attribute:: order

        :type: List or numpy array
        :default: ``[]``

        List of integer associated with the harmonic orders. The integer :math:`1` corresponds to the diagnostic of the fundamental field. If :math:`1` is not provided by the user, the diagnostic will still save the data for the fundamental field.


    .. py:attribute:: paraxial_max_factor

        :type: Float
        :default: ``0.5``

        The factor :math:`\alpha` is used to ensure the validity of the paraxial approximation. It defines the vector :math:`k` in the Hankel space such as :math:`k \in [0,\min(\pi/dr, \alpha q\omega/c )]`. See the :doc:`Far field propagation <../theory/vacuum_propagation>` section for more details.

    .. py:attribute:: every_r

        :type: Integer
        :default: ``1``

        Number of radial steps between two consecutive data saves. For example, if ``every_r=2``, only the deposited quantities at :math:`r=r_{2j}` are stored.

    .. py:attribute:: every_t

        :type: Integer
        :default: ``1``

        Number of time steps between two consecutive data saves. For example, if ``every_t=8``, only the deposited quantities at :math:`t=t_{8j}` are stored.

    .. py:attribute:: SI

        :type: Boolean
        :default: ``True``

        If ``True``, the field is saved in SI. Otherwise, it is saved in code units.
