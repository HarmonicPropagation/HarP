from typing import Callable, Any, TYPE_CHECKING, cast
import numpy as np
from numpy.typing import NDArray
import time
from components import HarpComponent

type SpatialGrid  = NDArray[np.float64]
type TimeGrid     = NDArray[np.float64]
type ComplexField = NDArray[np.complex128]

SVEA_WARN_RATIO  = .2 #a bit less than 5 cylces
#SVEA_ABORT_RATIO = .5 #approx 2 cycles

if TYPE_CHECKING:
    from simulation import Simulation
    from medium import Medium
    from components import Main

class Pulse(HarpComponent):  # type: ignore[misc]

    amplitude: float
    radial_profile: Callable[[SpatialGrid], SpatialGrid]
    time_envelope: Callable[[float | TimeGrid], float | TimeGrid]

    def __init__(self, **kwargs: Any) -> None:

        self.amplitude      = 1.0
        self.radial_profile = lambda r: np.ones_like(r)
        self.time_envelope  = lambda t: 1.0

        super(Pulse, self).__init__(**kwargs)

        self.E: ComplexField | None = None
        self.chi_at: complex | None = None
        self.k: Any = None
        self.dk: ComplexField | None = None
        self.kold: Any = None
        self.dkold: ComplexField | None = None
        self.compute_k: Callable[[], Any] | None = None
        self.compute_dk: Callable[[], ComplexField] | None = None
        self.medium: "Medium" | None = None
        self.chi_kerr: ComplexField | None = None
        self.timer: float = 0.
        #self.paraxial_check: bool = False

        self.Phase: ComplexField | None = None
        self.intPhase: Any = 0. + 0j
        self.Etilde: ComplexField | None = None
        self.Eold: ComplexField | None = None
    
    def check_time_svea(self, main: "Main") -> None:
        env = np.array([abs(self.time_envelope(t)) for t in main.t])
        denvdt     = np.gradient(env, main.t)
        keep = env>1e-5  #avoid numerical artefact
        if len(denvdt[keep])>0:
            svea_ratio = np.max(np.abs(denvdt[keep])/np.abs(env[keep]))
            if svea_ratio > SVEA_WARN_RATIO:
                main.printer(
                    f"WARNING: time SVEA is marginal \n"
                    f"              |dE/dt| / (omega * |E|) = {svea_ratio:.3f} (Warning for |dE/dt| / (omega * |E|) > {SVEA_WARN_RATIO:.1f}). \n"
                    f"  The separation of time scales the tabulated-dipole approach relies on breaks down at this pulse duration - see doc/source/theory/limitation.rst. \n"
                )
    
    def check_space_svea(self,main: "Main") -> None:
        assert self.E is not None
        assert self.Eold is not None
        dEdz=np.abs(self.E-self.Eold)/main.dz
        keep= np.abs(self.E)>1e-6 #avoid numerical artefact
        if len(dEdz[keep])>0:
            svea_ratio = np.max(dEdz[keep] / np.abs(self.E[keep]))
            if svea_ratio > SVEA_WARN_RATIO:
                #raise ValueError(
                main.printer(
                    f"WARNING: SVEA (Physical Constraint #5) violated: the pulse space profile should satisfy: \n"
                    f"            |dE/dz| / (k * |E|) = {svea_ratio:.3f}  << 1. (Warning for |dE/dz| / (k * |E|) > {SVEA_WARN_RATIO:.1f}.) \n"
                    f"  The equations relies on the paraxial approximation - see doc/source/theory/limitation.rst. \n"
                )
        
    def initialize(self, main: "Main", medium: "Medium") -> None:
        if main.approx_check:
            self.check_time_svea(main)
        self.medium = medium
        if main.anne_eq:
            self.compute_k  = self.compute_k_Anne
            self.compute_dk = self.compute_dk_Anne
        else:
            self.compute_k  = self.compute_k_AnneBar
            self.compute_dk = self.compute_dk_AnneBar

        assert self.medium is not None
        self.chi_at = self.medium.polarizability_0 + 0.j
        self.chi_kerr = np.zeros(main.nr, dtype=complex)

        self.k = self.compute_k()
        self.dk = np.zeros(main.nr, dtype=complex)

        self.kold = 1. + 0j
        self.dkold = np.zeros(main.nr, dtype=complex)

        self.Phase = np.zeros(main.nz, dtype=complex)
        self.intPhase = 0. + 0j

        self.E      = np.zeros(main.nr, dtype=complex)
        self.E[:]   = self.radial_profile(main.r) * self.time_envelope(main.t[0]) * self.amplitude

        if main.boundary_condition=="reflective":
            self.E[-1]=0.+0j

        if main.predictor_corrector:
            self.Etilde = np.zeros(main.nr, dtype=complex)
        
        if main.approx_check:
            self.paraxial_check = True
            self.Eold = np.zeros(main.nr, dtype=complex)

    def compute_chi_at(self) -> None:
        assert self.medium is not None
        self.chi_at   = self.medium.polarizability_0 + 0.j

    def compute_chi_kerr(self, Efield: ComplexField) -> None:
        assert self.medium is not None
        self.chi_kerr = 3/4. * self.medium.gamma_kerr * np.abs(Efield)**2 + 0.j

    def compute_k_Anne(self) -> Any:
        assert self.medium is not None
        assert self.chi_at is not None
        return 1. + self.chi_at * self.medium.n_gaz / 2. + 0.j

    def compute_dk_Anne(self) -> ComplexField:
        assert self.medium is not None
        assert self.medium.n_electron is not None
        assert self.medium.n_atomic is not None
        assert self.chi_at is not None
        assert self.chi_kerr is not None

        result = - self.chi_at * self.medium.n_electron / 2. - self.medium.n_electron / 2. + self.chi_kerr * self.medium.n_atomic / 2. + 0.j
        
        return cast(ComplexField, result)

    def compute_k_AnneBar(self) -> Any:
        return 1. + 0.j

    def compute_dk_AnneBar(self) -> ComplexField:
        assert self.medium is not None
        assert self.medium.n_electron is not None
        assert self.medium.n_atomic is not None
        assert self.chi_at is not None
        assert self.chi_kerr is not None

        result = (self.chi_at + self.chi_kerr) * self.medium.n_atomic / 2. - self.medium.n_electron / 2. + 0.j

        return cast(ComplexField, result)

    def reinitialize(self, sim: "Simulation", it: int) -> None:
        assert self.E is not None
        assert self.dk is not None
        assert self.Phase is not None
        assert self.compute_k is not None
        assert self.compute_dk is not None

        self.E[:] = self.radial_profile(sim.main.r) * self.time_envelope(sim.main.t[it]) * self.amplitude
        if sim.main.boundary_condition== "reflective":
            self.E[-1]=0.+0j
        
        self.compute_chi_at()
        self.compute_chi_kerr(self.E)
        self.k      = self.compute_k()
        self.dk[:]  = self.compute_dk()
        self.intPhase = 0. + 0j
        self.Phase[0] = 0. + 0j

    def update(self, sim: "Simulation", iz: int) -> None:
        self.timer = time.time()
        assert self.E is not None
        assert self.dk is not None
        assert self.dkold is not None
        assert self.Phase is not None
        assert self.compute_k is not None
        assert self.compute_dk is not None

        #self.compute_chi_at()
        self.compute_chi_kerr(self.E)
        self.kold     = self.k
        self.dkold[:] = self.dk[:]
        self.k        = self.compute_k()
        self.dk[:]    = self.compute_dk()

        if sim.main.predictor_corrector:
            self.compute_predictor_corrector(sim)
        
        if sim.main.approx_check:
            assert self.Eold is not None
            self.Eold[:]=self.E[:]
            self.E[:] = sim.solver.solve(self)
            self.check_space_svea(sim.main)
        else:
            self.E[:] = sim.solver.solve(self)
        
        sim.main.pulse_timer += time.time() - self.timer

        self.intPhase += self.k * sim.main.dz
        self.Phase[iz+1] = self.intPhase

    def compute_predictor_corrector(self, sim: "Simulation") -> None:
        assert self.Etilde is not None
        assert self.E is not None
        assert self.dk is not None
        assert self.compute_dk is not None

        self.Etilde[:] = sim.solver.solve(self)
        self.compute_chi_kerr((self.Etilde + self.E) / 2.)
        self.dk[:] = self.compute_dk()
