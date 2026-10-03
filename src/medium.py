from typing import Callable, Any
from components import HarpComponent
import numpy as np
from numpy.typing import NDArray
from scipy.interpolate import interp1d
import time

type SpatialGrid = NDArray[np.float64]     # e.g., main.r
type DensityField = NDArray[np.float64]    # e.g., n_gaz, n_atomic, n_electron
type FormFactor = NDArray[np.float64]      # e.g., data_f1, data_f2
type PhotonEnergy = NDArray[np.float64]    # e.g., data_qf1, data_qf2

type GazProfileFunc = Callable[[SpatialGrid, float], DensityField]


LOW_DENSITY_WARN_RATIO = 0.2

class Medium(HarpComponent): # type: ignore[misc]
    gaz_profile: GazProfileFunc | None = None
    polarizability_0_au: float | None = None  # value in a.u.
    polarizability_0: float = 1.8247e-40      # value in SI
    gamma_kerr_au: float | None = None        # value in a.u.
    gamma_kerr: float = 0.                    # value in SI
    table_path: str = ""
    skiprows: int = 3
    unit_ev: float = 1.
    unit_f1: float = 1.
    unit_f2: float = 1.

    def __init__(self, **kwargs: Any) -> None:
        super(Medium, self).__init__(**kwargs)
        
        self.n_atomic: DensityField | None = None
        self.n_gaz: DensityField | None = None
        self.n_electron: DensityField | None = None

        self.data_f1: FormFactor | None = None
        self.data_f2: FormFactor | None = None
        self.data_qf1: PhotonEnergy | None = None
        self.data_qf2: PhotonEnergy | None = None

        self.compute_f1: interp1d | None = None
        self.compute_f2: interp1d | None = None

        self.timer: float | None = None

    def initialize(self, main: Any) -> None:

        required_attrs = ['polarizability_0', 'gamma_kerr', 'gaz_profile']
        for attr in required_attrs:
            if getattr(self, attr, None) is None:
                raise ValueError(f"Medium initialization failed: '{attr}' must be defined.")
        
        assert self.gaz_profile is not None
        
        if main.approx_check:
            self.check_low_density(main)

        self.n_gaz = self.gaz_profile(main.r, main.z[0])

        self.n_atomic = np.ones(main.nr)
        self.n_atomic[:] = self.n_gaz
        self.n_electron = np.zeros(main.nr)

        self.polarizability_0 = self.polarizability_0 * main.Pol_r
        if self.polarizability_0_au is not None:
            self.polarizability_0 = self.polarizability_0_au * main.Polau_r

        self.gamma_kerr = self.gamma_kerr * main.Ker_r
        if self.gamma_kerr_au is not None:
            self.gamma_kerr = self.gamma_kerr_au * main.Kerau_r

        self.timer = 0.
        self.load_table(main)


    def check_low_density(self, main: Any) -> None:
        assert self.gaz_profile is not None
        ratio = max(np.max(self.gaz_profile(main.r, z)) for z in main.z)
        if ratio > LOW_DENSITY_WARN_RATIO:
            main.printer(
                f"       WARNING: low-density assumption is strained - max(n_gaz)/N_r = {ratio:.3g} "
                f"(recommended << 1, see doc/source/theory/limitation.rst, Physical Constraint #8)."
            )

    def load_table(self, main: Any) -> None:
        try:
            data = np.loadtxt(self.table_path, skiprows=self.skiprows, usecols=[0, 1, 2])
        except Exception:
            raise Exception("   ERROR in the harmonic table path: file not found    ")

        self.data_qf1  = data[:, 0] * self.unit_ev / main.ev_r # omega in omega_r (equivalently q order)
        self.data_qf2  = data[:, 0] * self.unit_ev / main.ev_r
        self.data_f1   = data[:, 1] * self.unit_f1
        self.data_f2   = data[:, 2] * self.unit_f2

        change = np.abs(self.data_f1) > 100
        self.data_f1[change] = -self.polarizability_0 * self.data_qf1[change]**2

        self.compute_f1 = interp1d(self.data_qf1[:], self.data_f1[:], kind='linear', bounds_error=False, fill_value=0.)
        self.compute_f2 = interp1d(self.data_qf2[:], self.data_f2[:], kind='linear', bounds_error=False, fill_value=0.)

    def reinitialize(self, sim: Any) -> None:
        assert self.gaz_profile is not None
        self.n_gaz      = self.gaz_profile(sim.main.r, sim.main.z[0])
        self.n_atomic   = self.n_gaz * np.exp(-sim.ionisation.integrated_rate[0])
        self.n_electron = self.n_gaz - self.n_atomic


    def update(self, sim: Any, iz: int) -> None:
        self.timer      = time.time()

        assert self.gaz_profile is not None
        self.n_gaz      = self.gaz_profile(sim.main.r, sim.main.z[iz])
        self.n_atomic   = self.n_gaz * np.exp(-sim.ionisation.integrated_rate[iz])
        self.n_electron = self.n_gaz - self.n_atomic

        sim.main.medium_timer += time.time() - self.timer
